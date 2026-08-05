/*
 * This class finds the single molecule spots using the mapping.
 */

import com.fasterxml.jackson.databind.ObjectMapper;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.Concatenator;
import ij.plugin.filter.MaximumFinder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, headless = true,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Spot Finder", weight = 2.0)})
public class smFRETSpotFinder implements Command {
    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter(description = "image stack to analyze", label = "Image", style = "open")
    File inputImageName;

    @Parameter(description = "channel to channel mapping JSON file", label = "Mapping JSON file", style = "open")
    File mappingFile;

    @Parameter (description = "first slice for averaging", min = "1")
    Integer startSlice = 1;

    @Parameter (description = "last slice for averaging", min = "1")
    Integer endSlice = 30;

    @Parameter (description = "spot SNR detection threshold", min = "1.0")
    Double spotThreshold = 6.0;

    @Parameter (description = "spot tolerance threshold (for MaximaFinder plugin)", min = "1.0")
    Double spotTolerance = 10.0;

    @Parameter (description = "spot prominence compared to pixels on 2 x spot Sigma radius", min = "1.0")
    Double spotProminence = 1.8;

    @Parameter (description = "spot size (sigma, pixels)", min = "0.2")
    Double spotSigma = 2.0;

    @Parameter (description = "camera offset / black level", min = "0")
    Integer cameraBlackLevel = 0;

    @Parameter (description = "camera gain (e-/ADU)", min = "0.1")
    Double cameraGain = 1.0;

    @Parameter (description = "minimum allowed distance between spots (pixels)", min = "1")
    Integer spotSpacing = 3;

    @Parameter (description = "margin around the edge of the channels (pixels)", min = "1")
    Integer edgeMargin = 5;

    // Radius masked as foreground around each spot, in pixels. Not adjustable,
    // because sweeping it against a known background over simulated PSFs from
    // sigma 1 to 3 found the choice barely matters: anything from 2 to 6 costs
    // at most 11% more error at any spot size, and 4 is within 2% of the best
    // everywhere. That is a consequence of how the background is estimated -
    // the clip finds contaminated pixels for itself, so the mask is no longer
    // the thing standing between spot light and the estimate, and there is
    // little left for a careful choice to buy.
    private static final int spotMargin = 4;

    @Parameter (description = "background clipping threshold, 0 to derive it from the spot size", min = "0.0")
    Double backgroundKappa = 0.0;

    // Rounds of clipping in backgroundEstimate. It settles in three or four -
    // each round removes the brightest leftovers and the ones after that find
    // nothing new to remove.
    private static final int backgroundClipRounds = 4;

    // kappa = backgroundKappaIntercept + backgroundKappaSlope * spotSigma, when
    // backgroundKappa is left at zero. Measured by sweeping radius, kappa and
    // smoothing against a known background over simulated aberrated Airy PSFs
    // from sigma 1 to 3: the best kappa runs 1.25, 1.00, 1.00, 0.80, 0.60 and
    // fits this line to an rms of 0.05.
    //
    // It *falls* with spot size, which is not the obvious direction. A wider
    // PSF spreads the same wing photons over more pixels, so the contamination
    // is a smaller excursion above the background in each one, and catching it
    // needs a tighter clip rather than a looser one.
    private static final double backgroundKappaIntercept = 1.5;
    private static final double backgroundKappaSlope = -0.3;

    // Smallest derived kappa. The fit reaches 0.6 at sigma 3 and would keep
    // going; below about half a robust sigma the clip starts removing ordinary
    // background noise along with the spot light.
    private static final double backgroundKappaFloor = 0.5;

    // Smoothing scale for the background estimate, in pixels. A constant, and
    // measured to be one: the best value sat at 12 to 16 pixels at every spot
    // size tried, with no trend against it (the fit is 14.4 - 0.00 * sigma).
    // It is a property of how fast the illumination varies, which is a fact
    // about the microscope rather than about the spots.
    //
    // This used to be 2 * spotMargin, which with the default margin gave 8 and
    // cost 40 to 60% more error than 14 at every sigma above 1.
    //
    // If this is ever changed, backgroundDecimation has to be revisited with it -
    // the two are tied, and the tie is not visible from either one alone.
    private static final double backgroundSmoothing = 14.0;

    // Grid the background estimate is smoothed on, as a divisor of the full
    // resolution: the estimate is binned by this, smoothed at
    // backgroundSmoothing / backgroundDecimation, and interpolated back. 1 would
    // smooth at full resolution.
    //
    // Worth 2.0x on stage 3 of a 1295 frame movie - 209 s to 103 s - for a
    // background estimate that moves by an rms of 0.02 ADU on a level of 16,
    // no spot moved on either test movie, and traces shifted by a median of
    // 0.16%. For scale, removing the 8 bit quantization moved traces by 2.9%.
    //
    // 8 was measured and rejected: it saves a further 9 s and changes which
    // molecules get measured (496 spots to 495 on the long movie, 2 lost and 1
    // gained). Time is the cheaper thing to spend.
    //
    // WARNING: this is only safe because backgroundSmoothing is 14, which leaves
    // 3.5 on the coarse grid. Binning is a box prefilter, so it adds
    // backgroundDecimation^2 / 12 to the variance of the effective kernel - 1.3
    // square pixels against sigma 14's 196, which is why the estimate barely
    // moves. Lower backgroundSmoothing without lowering this and that ratio
    // degrades as the square: at sigma 7 it is already 2.7% of the kernel
    // variance and at sigma 3.5 it is 11%. Keep backgroundSmoothing /
    // backgroundDecimation at 3.5 or above. There is nothing to lose by doing so
    // - below about that sigma imglib2's convolution is dominated by its own
    // per-call overhead rather than by the kernel, so a smaller coarse-grid sigma
    // buys no speed either. This does not touch spotSigma: the decimation lives
    // in maskedSmooth, which only ever runs at backgroundSmoothing.
    private static final int backgroundDecimation = 4;

    // Member variables.
    public ImagePlus backgroundMask;
    public java.util.List<String> columnHeaders = Arrays.asList("x", "y", "snr", "prominence"); // the first two fields should always be "x","y".
    private final boolean diagnostic_mode = true;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    public ImagePlus overlapMask;
    private String saveRootName;
    private final smFRETChannelMapper smfcm = new smFRETChannelMapper();  // smFREChannelMapper object.

    /**
     * Estimate background of an image.
     *
     * The background is a Gaussian weighted mean of the pixels we trust - inside the overlap
     * region and away from a spot - with the rest filled in by the same weighting. Then the
     * pixels sitting more than backgroundKappa robust standard deviations *above* that estimate
     * are dropped and it is computed again, a few times over.
     *
     * That second part is the point. Telling the estimator where the spots are is not enough,
     * because a real PSF has wings that reach well past any masking radius a crowded field can
     * afford, and the light in them lands in the pixels the estimator was told to trust. It
     * then reads the background high exactly where the background is about to be subtracted,
     * so it is subtracted twice. Letting the estimator find the contaminated pixels for itself
     * is what fixes that, and it makes the masking radius matter much less.
     *
     * The clipping is deliberately one-sided. Contamination only ever makes a pixel brighter,
     * so rejecting symmetrically would throw away good dark pixels and pull the estimate down.
     *
     * One-sidedness has a price, and it is paid back explicitly - see truncationConstants. Keeping
     * only the pixels below kappa spread leaves a truncated sample whose mean sits below the true
     * background, so the estimate reads low even where there is nothing to clip but noise. Measured
     * at about 1 ADU, roughly 2.8 times the single-round prediction because four rounds compound it.
     * That matters out of all proportion to its size: a trace is 4*pi*sigma^2 times frame minus
     * background, so 1 ADU is 50 units of trace, and spotFilterSNR scales the same error by its own
     * 2*pi*sigma^2 while the signal it is compared against grows only as sigma. The result was
     * an SNR that climbed with spot size - +15% from sigma 1 to 3 at fixed true significance, and
     * +46% at fixed molecular brightness - so spotThreshold silently meant something different on
     * every movie. Correcting the level each round rather than once at the end is what keeps the
     * clip in the right place and stops the bias compounding; it takes the drift to 0.5%.
     *
     * The estimate is returned as float whatever the input was. It used to be rounded back into
     * the input's type, which on an 8 bit movie meant the background was quantized to whole ADU -
     * and since a trace is 4*pi*sigma^2 times the difference between the frame and the background,
     * one ADU there is 50 units of trace. That quantization was carrying about 0.3 ADU of rounding
     * noise into every measurement, and it made any small disagreement in the estimate land as a
     * full step rather than a small one.
     *
     * Neither the clipping threshold nor the smoothing scale is asked for. The threshold comes
     * from spotSigma - see clippingThreshold - and the smoothing is a constant. Both were swept
     * against a known background over simulated PSFs from sigma 1 to 3, and what came out is
     * that only one of the three settings genuinely follows the spot size:
     *
     *   - the masking radius hardly matters. Anything from 2 to 6 pixels costs at most 11% more
     *     error at any sigma, and 4 is within 2% everywhere. That is this estimator working as
     *     intended: once the clip finds the contaminated pixels itself, there is little left for
     *     the mask to do, so spotMargin no longer needs to be chosen with care;
     *   - the threshold is sharp, and one grid step away costs 15 to 40%. It is worth deriving;
     *   - the smoothing scale showed no trend against sigma at all. It is set by how fast the
     *     illumination varies, not by the spots.
     */
    public ImagePlus backgroundEstimate(ImagePlus image) {
        return backgroundEstimate(image, null).image;
    }

    /**
     * As backgroundEstimate(image), but starting from a level a previous frame settled on.
     *
     * Consecutive frames of a movie share all but one frame of the temporal average behind
     * them, and their backgrounds come out nearly identical - measured at 0.035 ADU mean
     * change, with 3.5% of pixels moving at all. Recomputing the level from nothing each time
     * therefore spends four or five sigma 14 smoothings rediscovering what the last frame
     * already knew.
     *
     * It is the level that carries over, not the mask. The clip only ever removes pixels, so
     * a mask handed from frame to frame would erode all the way down a long movie with
     * nothing to restore it. The mask is rebuilt from the trusted pixels every frame; only
     * the starting estimate is inherited, and one round of clipping against it lands where
     * five rounds from scratch would.
     */
    public Background backgroundEstimate(ImagePlus image, float[] seed) {
        ImageProcessor original = image.getProcessor();

        // duplicate, because convertToFloatProcessor hands back the same object when the input is
        // already float, and Shared would then be aliasing the caller's pixels.
        Shared source = new Shared((FloatProcessor) original.convertToFloatProcessor().duplicate());
        float[] values = source.pixels;

        boolean[] keep = trustedPixels(source.width, source.height);
        float[] scratch = new float[values.length];
        double smoothing = backgroundSmoothing;
        double kappa = clippingThreshold();

        Shared estimate = null;
        float[] level = ((seed != null) && (seed.length == values.length)) ? seed : null;

        // Whether keep has been through the clip, and so whether the statistics read off it are
        // the truncated ones. The first round clips nothing, so its level and spread are honest.
        boolean clipped = false;
        truncationConstants(kappa);

        for (int round = 0; round < backgroundClipRounds; round++) {

            // Without a seed the first round has to smooth before it can clip. With one it clips
            // straight away, against a level that is already close to this frame's answer.
            if (level == null) {
                estimate = maskedSmooth(source, keep, smoothing, scratch);
                level = estimate.pixels;
            }

            double spread = robustSpread(values, level, keep, scratch);
            if (spread <= 0.0) {
                break;
            }
            if (clipped) {
                spread /= cachedMadFactor;
            }

            boolean[] fresh = new boolean[keep.length];
            boolean changed = false;
            boolean any = false;
            for (int i = 0; i < keep.length; i++) {
                fresh[i] = keep[i] && (values[i] - level[i]) < kappa * spread;
                changed |= (fresh[i] != keep[i]);
                any |= fresh[i];
            }
            if (!any || !changed) {
                break;
            }

            keep = fresh;
            clipped = true;
            estimate = maskedSmooth(source, keep, smoothing, scratch);
            level = estimate.pixels;

            // The level is now the mean of a sample truncated at kappa spread above it, which sits
            // low by a known amount. Putting it back here rather than at the end keeps the next
            // round's clip in the right place, which is what stopped the bias compounding.
            float bump = (float) (cachedMeanOffset * spread);
            for (int i = 0; i < level.length; i++) {
                level[i] += bump;
            }
        }

        // Reachable when a seeded first round finds nothing to clip: the level then still belongs
        // to the previous frame and has to be replaced, or the estimate would never again be
        // derived from the frame it is subtracted from.
        if (estimate == null) {
            estimate = maskedSmooth(source, keep, smoothing, scratch);
        }

        ImagePlus backgroundImage = new ImagePlus("background_image", estimate.processor);
        return new Background(backgroundImage, estimate.pixels);
    }

    /**
     * A background estimate and the level it settled on, so the next frame can start from it.
     */
    public static final class Background {
        public final ImagePlus image;
        public final float[] level;

        Background(ImagePlus image, float[] level) {
            this.image = image;
            this.level = level;
        }
    }

    /**
     * An imglib2 image and an ImageJ1 processor over one array of pixels.
     *
     * Two things in this class have to stay ImageJ1, because no imglib2 equivalent produces the
     * same numbers: GaussianBlur, which differs from Gauss3 by 0.04 ADU at the spot scale and
     * 0.8 ADU - 3.9% - at the background scale, where the kernel is downscaled; and the ROI
     * rasterizers, where fillOval covers 69 pixels at radius 4 against 49 for the analytic disc.
     * Sharing the pixel array lets those two run in place while everything else is imglib2, with
     * no conversion between them.
     */
    /**
     * Gaussian smoothing, imglib2's.
     *
     * This was ImageJ1's GaussianBlur. The two do not agree exactly - measured at 0.04 ADU at the
     * spot scale and 0.8 ADU, 3.9%, at the sigma 14 background scale, where ImageJ1 downscales its
     * kernel - so numbers derived under the old one shifted slightly when this changed.
     *
     * The out of bounds strategy is the caller's, because it is not the same question everywhere:
     * a normalized convolution wants zero outside, so that absent pixels carry no weight in either
     * the numerator or the denominator, while plain smoothing wants mirroring, so that the edge is
     * not dragged toward zero.
     */
    static void gauss(double sigma, RandomAccessible<FloatType> source, RandomAccessibleInterval<FloatType> target) {
        try {
            Gauss3.gauss(sigma, source, target);
        } catch (IncompatibleTypeException e) {
            // Both sides are FloatType, so this cannot happen.
            throw new smFRETAnalysisException("Error: gaussian smoothing failed", e);
        }
    }

    static final class Shared {
        final int height;
        final ArrayImg<FloatType, FloatArray> img;
        final float[] pixels;
        final FloatProcessor processor;
        final int width;

        Shared(int width, int height) {
            this(new FloatProcessor(width, height));
        }

        Shared(FloatProcessor source) {
            width = source.getWidth();
            height = source.getHeight();
            pixels = (float[]) source.getPixels();
            img = ArrayImgs.floats(pixels, width, height);
            processor = source;
        }
    }

    private static double cachedKappa = Double.NaN;
    private static double cachedMeanOffset;
    private static double cachedMadFactor;

    /**
     * Normal CDF, via Abramowitz and Stegun 7.1.26. Good to about 1e-7, which is far past what
     * a clipping correction needs.
     */
    private static double normalCdf(double x) {
        double z = x / Math.sqrt(2.0);
        double sign = (z < 0.0) ? -1.0 : 1.0;
        double a = Math.abs(z);
        double t = 1.0 / (1.0 + 0.3275911 * a);
        double poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
                + t * (-1.453152027 + t * 1.061405429))));
        double erf = sign * (1.0 - poly * Math.exp(-a * a));
        return 0.5 * (1.0 + erf);
    }

    private static double normalPdf(double x) {
        return Math.exp(-0.5 * x * x) / Math.sqrt(2.0 * Math.PI);
    }

    /** Smallest x with normalCdf(x) >= p, by bisection. */
    private static double normalQuantile(double p) {
        double lo = -10.0;
        double hi = 10.0;
        for (int i = 0; i < 80; i++) {
            double mid = 0.5 * (lo + hi);
            if (normalCdf(mid) < p) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /**
     * The two things one-sided clipping does to the level, as multiples of the true spread.
     *
     * Keeping only residuals below kappa sigma leaves a truncated Gaussian, and both statistics
     * the estimator reads off it are biased:
     *
     *   - its *mean* - which is what maskedSmooth computes - sits phi(kappa)/Phi(kappa) below the
     *     true mean. That is the bias being corrected;
     *   - its *MAD* - which is what robustSpread reads - is smaller than an untruncated Gaussian's,
     *     so the spread comes out low and a correction scaled by it would under-correct. The MAD is
     *     taken about the truncated median, matching what robustSpread does, and solved numerically
     *     because it has no closed form.
     *
     * Both reduce to their untruncated values as kappa grows: the offset to zero, the MAD factor to
     * 1.4826 * 0.6745 = 1.
     */
    private static void truncationConstants(double kappa) {
        if (kappa == cachedKappa) {
            return;
        }
        double mass = normalCdf(kappa);
        cachedMeanOffset = normalPdf(kappa) / mass;

        // Median of the truncated distribution, which is what robustSpread centres on.
        double median = normalQuantile(0.5 * mass);

        // Half-width containing half the truncated mass about that median.
        double lo = 0.0;
        double hi = 10.0;
        for (int i = 0; i < 80; i++) {
            double d = 0.5 * (lo + hi);
            double covered = normalCdf(Math.min(median + d, kappa)) - normalCdf(median - d);
            if (covered < 0.5 * mass) {
                lo = d;
            } else {
                hi = d;
            }
        }
        cachedMadFactor = 1.4826 * 0.5 * (lo + hi);
        cachedKappa = kappa;
    }

    /**
     * How far above the estimate a pixel may sit before it is called spot light.
     *
     * Derived from spotSigma unless backgroundKappa says otherwise, because it is not an
     * independent property of the experiment - it is set by how much of a spot's light lands
     * in the pixels the estimator is trusting, and that follows from the spot's size.
     *
     * Setting backgroundKappa to anything above zero overrides this. That escape hatch is
     * there because the relationship was measured on simulated data with a smooth Gaussian
     * illumination profile, and the one real movie it has been checked against wanted a looser
     * clip than the line predicts - 1.8 against 0.9. Two things could explain that and they
     * have not been separated: that movie is analysed with a spotsigma of 2.0 while its PSF
     * actually fits at 1.36, so the input to the formula is wrong there; and a structured
     * illumination profile has real background variation that hard clipping would eat.
     */
    public double clippingThreshold() {
        if (backgroundKappa > 0.0) {
            return backgroundKappa;
        }
        return Math.max(
            backgroundKappaFloor,
            backgroundKappaIntercept + backgroundKappaSlope * spotSigma);
    }

    /**
     * Pixels the background estimate can be built from: inside the overlap and off a spot.
     */
    private boolean[] trustedPixels(int width, int height) {
        RandomAccessibleInterval<FloatType> overlap = ImageJFunctions.convertFloat(overlapMask);
        RandomAccessibleInterval<FloatType> background = ImageJFunctions.convertFloat(backgroundMask);

        boolean[] keep = new boolean[width * height];
        int[] index = {0};
        LoopBuilder.setImages(Views.flatIterable(overlap), Views.flatIterable(background))
                .forEachPixel((o, b) -> keep[index[0]++] = (o.get() > 0.0f) && (b.get() > 0.0f));
        return keep;
    }

    /**
     * Gaussian weighted mean of the kept pixels, ignoring the rest.
     *
     * Smoothing the image and the mask with the same kernel and dividing one by the other is
     * a normalized convolution: every pixel gets the average of its kept neighbours, weighted
     * by distance, whether or not it was kept itself. That fills the spots and the edges in
     * one pass, which is what the two separate inpainting passes used to do.
     */
    private Shared maskedSmooth(Shared image, boolean[] keep, double sigma, float[] scratch) {
        Shared weighted = new Shared(image.width, image.height);
        Shared total = new Shared(image.width, image.height);

        for (int i = 0; i < image.pixels.length; i++) {
            if (keep[i]) {
                weighted.pixels[i] = image.pixels[i];
                total.pixels[i] = 1.0f;
            }
        }

        // Zero outside, which is what makes this a normalized convolution rather than a smoothing:
        // a pixel off the edge contributes nothing to the weighted sum and nothing to the weight,
        // so the ratio stays the mean of the kept pixels that actually exist.
        Shared weightedSmooth = new Shared(image.width, image.height);
        Shared totalSmooth = new Shared(image.width, image.height);
        if (backgroundDecimation > 1) {
            smoothDecimated(weighted, sigma, backgroundDecimation, weightedSmooth);
            smoothDecimated(total, sigma, backgroundDecimation, totalSmooth);
        } else {
            gauss(sigma, Views.extendZero(weighted.img), weightedSmooth.img);
            gauss(sigma, Views.extendZero(total.img), totalSmooth.img);
        }

        // Only reached by a pixel with no kept neighbour anywhere in the kernel, which needs a
        // masked patch several times the smoothing scale across. It is a floor, not a fill.
        final float floor = (float) subsetMedian(image.pixels, keep, scratch);

        Shared smoothed = new Shared(image.width, image.height);
        LoopBuilder.setImages(smoothed.img, weightedSmooth.img, totalSmooth.img)
                .forEachPixel((out, sum, count) ->
                        out.set((count.get() > 1.0e-6f) ? (sum.get() / count.get()) : floor));
        return smoothed;
    }

    /**
     * Gaussian smoothing done on a coarser grid and interpolated back.
     *
     * Gauss3 is a direct separable convolution - it truncates the kernel at three sigma and
     * convolves, with no FFT anywhere - so it costs pixels times sigma, measured flat at 0.10 ns
     * per tap-pixel from sigma 7 to sigma 28. Decimating by d cuts both factors: d squared fewer
     * pixels and a kernel d times narrower. The blur alone therefore drops as d cubed, 3.34 ms to
     * 0.26 ms at d of 4.
     *
     * The blur is not the whole cost, which is why d of 4 is the knee rather than a waypoint.
     * Binning and interpolating back are both full resolution passes and do not shrink with d, so
     * by d of 4 the three phases cost about the same and the total has fallen 3.34 ms to 0.75 ms -
     * a factor of 4.5, not 64. Past that there is little of the blur left to remove.
     *
     * It is an approximation, but a mild and quantifiable one. Summing d by d blocks is a box
     * prefilter, which adds d squared over twelve to the variance of the effective kernel - 1.3
     * square pixels at d of 4 against sigma 14's 196, so the kernel this actually applies is
     * sigma 14.05. The estimate is a surface that varies on a 14 pixel scale being sampled every
     * d pixels, which is far above what it needs. See backgroundDecimation for what that ratio
     * costs if sigma is lowered without lowering d.
     *
     * Both the weighted sum and the weight are summed rather than averaged, and each is
     * decimated the same way, so the d squared factors cancel when the caller divides one by
     * the other. Blocks that hang off the right or bottom edge are simply short, which the
     * weight image already accounts for.
     */
    private static void smoothDecimated(Shared image, double sigma, int decimation, Shared target) {
        Shared binned = bin(image, decimation);
        Shared smoothed = new Shared(binned.width, binned.height);
        gauss(sigma / decimation, Views.extendZero(binned.img), smoothed.img);
        unbin(smoothed, decimation, target);
    }

    /**
     * Sum each decimation by decimation block of pixels into one.
     */
    private static Shared bin(Shared image, int decimation) {
        Shared binned = new Shared(
                (image.width + decimation - 1) / decimation,
                (image.height + decimation - 1) / decimation);
        for (int y = 0; y < image.height; y++) {
            int source = y * image.width;
            int row = (y / decimation) * binned.width;
            for (int x = 0; x < image.width; x++) {
                binned.pixels[row + (x / decimation)] += image.pixels[source + x];
            }
        }
        return binned;
    }

    /**
     * Interpolate a binned image back to full size.
     *
     * The block starting at full resolution x covers x to x + decimation - 1, so its centre - the
     * point the binned pixel actually stands for - is half a block in from its corner. Ignoring
     * that offset would shift the whole background estimate by (decimation - 1) / 2 pixels.
     */
    private static void unbin(Shared binned, int decimation, Shared target) {
        int width = target.width;
        int height = target.height;

        int[] left = new int[width];
        int[] right = new int[width];
        float[] acrossWeight = new float[width];
        interpolationWeights(width, binned.width, decimation, left, right, acrossWeight);

        int[] above = new int[height];
        int[] below = new int[height];
        float[] downWeight = new float[height];
        interpolationWeights(height, binned.height, decimation, above, below, downWeight);

        float[] source = binned.pixels;
        for (int y = 0; y < height; y++) {
            int upper = above[y] * binned.width;
            int lower = below[y] * binned.width;
            float fy = downWeight[y];
            int row = y * width;
            for (int x = 0; x < width; x++) {
                int a = left[x];
                int b = right[x];
                float fx = acrossWeight[x];
                float top = source[upper + a] + fx * (source[upper + b] - source[upper + a]);
                float bottom = source[lower + a] + fx * (source[lower + b] - source[lower + a]);
                target.pixels[row + x] = top + fy * (bottom - top);
            }
        }
    }

    /**
     * The two binned samples each full resolution pixel sits between, and how far between them.
     *
     * Separable, so one pass along each axis is enough and the weights are reused down every row.
     * This is deliberately arithmetic rather than an imglib2 NLinearInterpolator: interpolating
     * through a RealRandomAccess costs a setPosition and a get per pixel, which measured at 4.8 ms
     * on a 256 by 512 frame - more than the full resolution sigma 14 blur it exists to avoid.
     *
     * Positions off either end clamp to the last sample, matching extendBorder.
     */
    private static void interpolationWeights(int size, int binnedSize, int decimation,
                                             int[] lower, int[] upper, float[] weight) {
        double centre = 0.5 * (decimation - 1);
        for (int p = 0; p < size; p++) {
            double position = (p - centre) / decimation;
            int index = (int) Math.floor(position);
            double fraction = position - index;
            if (index < 0) {
                index = 0;
                fraction = 0.0;
            } else if (index >= (binnedSize - 1)) {
                index = Math.max(binnedSize - 2, 0);
                fraction = (binnedSize > 1) ? 1.0 : 0.0;
            }
            lower[p] = index;
            upper[p] = Math.min(index + 1, binnedSize - 1);
            weight[p] = (float) fraction;
        }
    }

    /**
     * Standard deviation of the residual over the kept pixels, from the median absolute deviation.
     *
     * The ordinary standard deviation of a contaminated sample is inflated by exactly the pixels
     * being looked for, which would make the clip too loose to catch them.
     */
    private double robustSpread(float[] values, float[] level, boolean[] keep, float[] scratch) {
        int count = 0;
        for (int i = 0; i < keep.length; i++) {
            if (keep[i]) {
                scratch[count++] = values[i] - level[i];
            }
        }
        if (count == 0) {
            return 0.0;
        }

        double median = median(scratch, count);
        for (int i = 0; i < count; i++) {
            scratch[i] = (float) Math.abs(scratch[i] - median);
        }
        return 1.4826 * median(scratch, count);
    }

    /**
     * Median of the kept entries of values. Uses scratch as working space.
     */
    private double subsetMedian(float[] values, boolean[] keep, float[] scratch) {
        int count = 0;
        for (int i = 0; i < keep.length; i++) {
            if (keep[i]) {
                scratch[count++] = values[i];
            }
        }
        return (count == 0) ? 0.0 : median(scratch, count);
    }

    /**
     * Median of the first count entries, which are partially reordered in place.
     *
     * Selection rather than a sort. This is called twice per clipping round, on up to every pixel
     * in the frame, so on a long movie it ran to a few thousand full sorts of 131072 floats -
     * 14% of the profile. Finding the middle element does not require ordering the rest, and both
     * routines leave the same multiset behind, which is what robustSpread relies on when it reuses
     * the scratch array for the absolute deviations.
     */
    private static double median(float[] values, int count) {
        int middle = count / 2;
        double high = select(values, count, middle);
        if ((count % 2) != 0) {
            return high;
        }

        // Everything below the selected position is no greater than it, so the element below the
        // median is the largest of them - no second selection needed.
        float low = values[0];
        for (int i = 1; i < middle; i++) {
            if (values[i] > low) { low = values[i]; }
        }
        return 0.5 * ((double) low + high);
    }

    /**
     * The k-th smallest of the first count entries, partitioning in place around it.
     *
     * Quickselect with a median-of-three pivot, which keeps it deterministic - the same input
     * gives the same work, so a run is reproducible - and keeps sorted or nearly sorted input,
     * which is the pathological case for a naive pivot, away from quadratic.
     */
    private static float select(float[] values, int count, int k) {
        int low = 0;
        int high = count - 1;

        while (low < high) {
            int middle = low + ((high - low) >> 1);
            if (values[middle] < values[low]) { swap(values, middle, low); }
            if (values[high] < values[low]) { swap(values, high, low); }
            if (values[high] < values[middle]) { swap(values, high, middle); }
            float pivot = values[middle];

            int i = low;
            int j = high;
            while (i <= j) {
                while (values[i] < pivot) { i++; }
                while (values[j] > pivot) { j--; }
                if (i <= j) {
                    swap(values, i, j);
                    i++;
                    j--;
                }
            }

            if (k <= j) { high = j; }
            else if (k >= i) { low = i; }
            else { return values[k]; }
        }
        return values[k];
    }

    /**
     * select helper.
     */
    private static void swap(float[] values, int a, int b) {
        float temp = values[a];
        values[a] = values[b];
        values[b] = temp;
    }

    /**
     * Count the number of good spots.
     */
    private int countGoodSpots(double[][] spots){
        int numGoodSpots = 0;

        for (int i=0; i < spots.length; i++) {
            if (spots[i][0] > 0.5){
                numGoodSpots += 1;
            }
        }
        return numGoodSpots;
    }

    /**
     * Creates the mask that is used to identify spots that are not in the overlap
     * region of both channels.
     *
     * Spots on the pixels where this image is <= overlapThreshold are filtered out.
     */
    private ImagePlus createOverlapMask (int imageWidth, int imageHeight){
        // If you use 16 bit images they don't get auto-scaled when converted back from float32, but 8
        // bit images do, why IDK. TurboReg, as used by smFRETChannelMapper splitImagePlut() method,
        // converts images to float32 when it transforms them.
        ImagePlus overlapMask = IJ.createImage("overlap_mask", "16-bit black", imageWidth, imageHeight, 1);

        // Draw filled rectangles in allowed region.
        ImageProcessor ip = overlapMask.getProcessor();
        int hw = overlapMask.getWidth()/2;
        ip.setColor(100);
        ip.setLineWidth(0);
        ip.fillRect(edgeMargin, edgeMargin, hw - 2*edgeMargin, overlapMask.getHeight() - 2*edgeMargin);
        ip.fillRect(hw + edgeMargin, edgeMargin, hw - 2*edgeMargin, overlapMask.getHeight() - 2*edgeMargin);
        overlapMask.updateAndDraw();

        // Transform and overlap allowed regions. The two halves are added and thresholded, so a
        // pixel survives only where both channels contributed their full 100.
        java.util.List<ImagePlus> overlapImages = smfcm.splitImagePlus(overlapMask, true);
        RandomAccessibleInterval<FloatType> targetHalf = ImageJFunctions.convertFloat(overlapImages.get(0));
        RandomAccessibleInterval<FloatType> sourceHalf = ImageJFunctions.convertFloat(overlapImages.get(1));

        Shared overlapSum = new Shared((int) targetHalf.dimension(0), (int) targetHalf.dimension(1));
        LoopBuilder.setImages(overlapSum.img, targetHalf, sourceHalf)
                .forEachPixel((out, a, b) -> out.set(((a.get() + b.get()) > 190.0f) ? 1.0f : 0.0f));

        return new ImagePlus("overlap_mask", overlapSum.processor);
    }

    /**
     * Creates mask for (circular) neighborhood around spots.
     *
     * Spots on the pixels where this image is > 1 are filtered out.
     */
    private ImagePlus createSpotsNeighborhoodMask(double[][] spots, int imageWidth, int imageHeight, int radius){
        Shared neighborhoodMask = new Shared(imageWidth, imageHeight);

        // One spot at a time into a scratch image, then accumulated - the count of overlapping
        // neighbourhoods is what the two readers of this mask key off, so they cannot simply be
        // drawn on top of each other. The drawing is ImageJ1 because its oval rasterization is
        // not the analytic disc and the mask footprint has to stay what it was; see Shared.
        Shared spot = new Shared(imageWidth, imageHeight);
        for (int i = 0; i < spots.length; i++) {
            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            Arrays.fill(spot.pixels, 0.0f);
            spot.processor.setColor(1);
            spot.processor.setLineWidth(0);
            spot.processor.fillOval(x-radius, y-radius, 2*radius+1, 2*radius+1);

            LoopBuilder.setImages(neighborhoodMask.img, spot.img)
                    .forEachPixel((total, one) -> total.set(total.get() + one.get()));
        }

        return new ImagePlus("neighborhood_mask", neighborhoodMask.processor);
    }

    /**
     * Return maxima in the image as double array.
     *
     * The first element for each spot is whether it passes the current filters.
     */
    private double[][] getMaxima(ImagePlus image){
        ImageProcessor imageProc = image.getProcessor();
        MaximumFinder mf = new MaximumFinder();
        Polygon spotsPoly = mf.getMaxima(imageProc, spotTolerance, true);

        double[][] spotsArr = new double[spotsPoly.npoints][3];
        for (int i=0; i < spotsPoly.npoints; i++){
            spotsArr[i][0] = 1.0;
            spotsArr[i][1] = spotsPoly.xpoints[i];
            spotsArr[i][2] = spotsPoly.ypoints[i];
        }
        return spotsArr;
    }

    /**
     * Copied from spotIntensityAnalysis plugin.
     *
     *  // AUTHOR:       Nico Stuurman
     *  //
     *  // COPYRIGHT:    University of California, San Francisco 2015
     *  //
     *  // LICENSE:      This file is distributed under the BSD license.
     *  //               License text is included with the source distribution.
     *  //
     *  //               This file is distributed in the hope that it will be useful,
     *  //               but WITHOUT ANY WARRANTY; without even the implied warranty
     *  //               of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
     *  //
     *  //               IN NO EVENT SHALL THE COPYRIGHT OWNER OR
     *  //               CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
     *  //               INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.
     */
    public static Overlay getSpotOverlay (double[][] spots, int radius, Color symbolColor) {
        Overlay ov = new Overlay();
        double diameter = 2.0 * (double)radius;
        for (int i = 0; i < spots.length; i++) {
            if (spots[i][0] > 0.5) {
                double x = spots[i][1] + 0.5;
                double y = spots[i][2] + 0.5;
                Roi roi = new Roi(x - radius, y - radius,
                        diameter, diameter, (int) diameter);
                roi.setStrokeColor(symbolColor);
                ov.add(roi);
            }
        }
        return ov;
    }

    /**
     * Load an existing mapping JSON file to initialize smFRETChannelMapper.
     */
    public void loadMappingJSON(String mappingFileName){
        // Set the log before loading, otherwise a load failure NPEs in smfcm's
        // exception handler instead of reporting the actual problem.
        smfcm.log = log;
        smfcm.loadMappingJSON(mappingFileName);
    }

    /**
     * Load overlay and background masks.
     */
    public void loadMasks(String masksFileName){
        ImagePlus masksImage = new ImagePlus(masksFileName);

        // Overlay mask.
        ImageProcessor ip = masksImage.getStack().getProcessor(1);
        overlapMask = new ImagePlus("overlap mask", ip);

        // Background mask.
        ip = masksImage.getStack().getProcessor(2);
        backgroundMask = new ImagePlus("background mask", ip);
    }

    /**
     * Load spot locations as a double[][].
     */
    public double[][] loadSpotLocations(String spotsFileName){
        try {
            ResultsTable rt = ResultsTable.open2(spotsFileName);
            double[][] spots = new double[rt.getCounter()][columnHeaders.size()];
            log.info("column size " + columnHeaders.size());

            for (int i = 0; i < rt.getCounter(); i++) {
                for (int j = 0; j < columnHeaders.size(); j++){
                    spots[i][j] = rt.getValue(columnHeaders.get(j), i);
                }
            }
            return spots;
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
        return null;
    }

    /**
     * converts neighborhood mask for use as a background mask.
     */
    private ImagePlus neighborhoodMaskToBackgroundMask(ImagePlus neighborhoodMask) {
        return invertNeighborhoodMask(neighborhoodMask, "foreground_mask", 0.0f);
    }

    /**
     * converts neighborhood mask for use to detect spots that are too close to each other.
     */
    private ImagePlus neighborhoodMaskToProximityMask(ImagePlus neighborhoodMask) {
        return invertNeighborhoodMask(neighborhoodMask, "overlap_mask", 1.0f);
    }

    /**
     * Mark the pixels where the neighbourhood count is at or below limit, and clear the rest.
     *
     * The two callers differ only in that limit: a pixel is background if no spot's neighbourhood
     * reached it at all, and is far enough from its neighbours if at most one did.
     */
    private ImagePlus invertNeighborhoodMask(ImagePlus neighborhoodMask, String title, float limit) {
        RandomAccessibleInterval<FloatType> counts = ImageJFunctions.convertFloat(neighborhoodMask);
        Shared inverted = new Shared((int) counts.dimension(0), (int) counts.dimension(1));

        LoopBuilder.setImages(inverted.img, counts)
                .forEachPixel((out, count) -> out.set((count.get() > limit) ? 0.0f : 1.0f));

        ImagePlus mask = new ImagePlus(title, inverted.processor);
        return mask;
    }

    /**
     * Save spot locations as a table.
     */
    private void saveSpotLocations(String spotsFileName, double[][] spots){
        try {
            ResultsTable rt = new ResultsTable();

            for (int i = 0; i < spots.length; i++) {
                // Only save good spots.
                if (spots[i][0] > 0.5) {
                    rt.incrementCounter();
                    for (int j=0; j < columnHeaders.size(); j++) {
                        rt.addValue(columnHeaders.get(j), spots[i][j+1]);
                    }
                }
            }
            rt.saveAs(spotsFileName);
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Thin wrapper around smfcm function.
     */
    public java.util.List<ImagePlus> splitImagePlus(ImagePlus image){
        return smfcm.splitImagePlus(image, true);
    }

    /**
     * Mark spots with prominence less than threshold.
     */
    private double[][] spotFilterProminence(double[][] spots, ImagePlus sumImage, ImagePlus backgroundImage) {
        double[][] filteredSpots = new double[spots.length][spots[0].length+1];
        int last_col = filteredSpots[0].length-1;

        int srad = (int)(Math.round(2.0*spotSigma));
        Shared fgImage = foreground(sumImage, backgroundImage);
        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            // Calculate spot lowest prominence over pixels in circular neighborhood.
            double spotHeight = valueAt(fgImage, x, y);
            double lowestProminence = spotHeight;
            for (int rx = -srad; rx <= srad; rx += 1){
                for (int ry = -srad; ry <= srad; ry += 1){
                    if (rx*rx+ry*ry >= (srad-1)*(srad-1) && rx*rx+ry*ry <= (srad+1)*(srad+1)){
                        double pixelHeight = Math.max(1, valueAt(fgImage, x+rx, y+ry));
                        double prominence = spotHeight/pixelHeight;
                        if (prominence < lowestProminence){
                            lowestProminence = prominence;
                        }
                    }
                }
            }
            filteredSpots[i][last_col] = lowestProminence;
            if (lowestProminence < spotProminence){
                filteredSpots[i][0] = 0.0;
            }
        }

        return filteredSpots;
    }

    /**
     * Mark spots with estimated SNR less than threshold.
     *
     * The factor of two for the camera black level is because we added the two channels together.
     */
    private double[][] spotFilterSNR(double[][] spots, ImagePlus sumImage, ImagePlus backgroundImage){
        double[][] filteredSpots = new double[spots.length][spots[0].length+1];
        int last_col = filteredSpots[0].length-1;

        Shared foreground = foreground(sumImage, backgroundImage);
        Shared fgSmooth = new Shared(foreground.width, foreground.height);
        gauss(spotSigma, Views.extendMirrorSingle(foreground.img), fgSmooth.img);

        Shared background = new Shared(
                (FloatProcessor) backgroundImage.getProcessor().convertToFloatProcessor().duplicate());
        Shared bgSmooth = new Shared(background.width, background.height);
        gauss(spotSigma, Views.extendMirrorSingle(background.img), bgSmooth.img);

        // The integrated magnitude of the spot, recovered from the peak of the smoothed image.
        // gauss() convolves with a normalized Gaussian, so a spot carrying N photons at the size
        // we are measuring at leaves N / (4 pi spotSigma^2) at its centre - the spot's width and
        // the kernel's add in quadrature - and this is the factor that gets N back.
        //
        // This was 2 pi spotSigma^2, which converts a normalized Gaussian to a unit height one but
        // leaves the width of the spot itself out of the quadrature sum, so it recovered N/2 and
        // took the noise term from an area half the right size. The ratio of the two was short of
        // the true matched filter significance by exactly sqrt(2); with this factor the number
        // reported here *is* that significance, so spotThreshold now reads in real sigma.
        //
        // Consequence: every SNR is sqrt(2) larger than it used to be, so a given spotThreshold
        // admits more spots than the same number did before. 6 now means 6 sigma, where it used
        // to mean 8.49.
        double norm = 4.0*Math.PI*spotSigma*spotSigma;
        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            double fg = norm*cameraGain*valueAt(fgSmooth, x, y);
            double bg = norm*cameraGain*(valueAt(bgSmooth, x, y) - 2*cameraBlackLevel);

            if (bg > 1){
                bg = Math.sqrt(bg);
            }
            else{
                bg = 1.0;
            }
            double snr = fg/bg;
            if (diagnostic_mode) {
                log.info(x + " " + y + " " + fg + " " + " " + bg + " " + snr);
            }
            filteredSpots[i][last_col] = snr;
            if (snr < spotThreshold) {
                filteredSpots[i][0] = 0.0;
            }
        }

        if (diagnostic_mode) {
            FileSaver fgSmoothImageSaver = new FileSaver(new ImagePlus("foreground_smooth", fgSmooth.processor));
            fgSmoothImageSaver.saveAsTiff(saveRootName + "_spotf_fg_smooth.tif");

            FileSaver bgSmoothImageSaver = new FileSaver(new ImagePlus("background_smooth", bgSmooth.processor));
            bgSmoothImageSaver.saveAsTiff(saveRootName + "_spotf_bg_smooth.tif");
        }

        return filteredSpots;
    }

    /**
     * The spot signal: the channel sum with the background estimate taken off, in float.
     *
     * This was ImageCalculator's "subtract create", which took its type from the first operand.
     * With the sum now float rather than 8 bit, negatives survive instead of clamping at zero and
     * bright pixels are no longer pinned at 255.
     */
    private Shared foreground(ImagePlus sumImage, ImagePlus backgroundImage) {
        RandomAccessibleInterval<FloatType> sum = ImageJFunctions.convertFloat(sumImage);
        RandomAccessibleInterval<FloatType> background = ImageJFunctions.convertFloat(backgroundImage);

        Shared foreground = new Shared((int) sum.dimension(0), (int) sum.dimension(1));
        LoopBuilder.setImages(foreground.img, sum, background)
                .forEachPixel((out, s, b) -> out.set(s.get() - b.get()));
        return foreground;
    }

    /**
     * Value at a pixel, or zero off the edge of the image.
     *
     * Reads used to go through ImagePlus.getPixel, which returned zero outside the image - the
     * prominence filter walks a ring around each spot and relies on that. It cannot be used here
     * any more regardless: on a float image getPixel returns the raw bits of the float, not the
     * value.
     */
    private static float valueAt(Shared image, int x, int y) {
        if ((x < 0) || (y < 0) || (x >= image.width) || (y >= image.height)) {
            return 0.0f;
        }
        return image.pixels[y * image.width + x];
    }

    /**
     * Mark spots where mask is 0.
     */
    private double[][] spotFilterWithMask(double[][]  spots, ImagePlus mask){
        double[][] filteredSpots = new double[spots.length][spots[0].length];
        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];
            // getf, not getPixel: on a float image getPixel hands back the raw bits of the float.
            if (mask.getProcessor().getf(x, y) == 0.0f){
                filteredSpots[i][0] = 0.0;
            }
        }

        return filteredSpots;
    }

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
	        log.info("starting spot finding on image " + inputImageName);

            // Root name to use for saving output, this is just the file name
            // without the extension.
            saveRootName = inputImageName.toString();
            int dotIndex = saveRootName.lastIndexOf('.');
            if (dotIndex > 0) {
                saveRootName = saveRootName.substring(0, dotIndex);
            }
            log.info("save root " + saveRootName);

            // Load the image to process.
            ImagePlus inputImage = new ImagePlus(inputImageName.toString());
            smFRETChannelMapper.requireGrayscale(inputImage, "the image " + inputImageName);
            smFRETChannelMapper.requireTimeStack(inputImage, "the image " + inputImageName);
            inputImage = smFRETChannelMapper.toFloat(inputImage);

            // Load the channel to channel mapping file.
            loadMappingJSON(mappingFile.toString());

            // Average image.
            log.info("average image - " + smFRETChannelMapper.frameCount(inputImage) + " frames");
            ImagePlus averageImage = smfcm.averageImagePlus(inputImage, startSlice, endSlice);

            // split, transform and add the two channels together.
            log.info("split and transform");
            java.util.List<ImagePlus> images = smfcm.splitImagePlus(averageImage, true);

            // Added in float. This was ImageCalculator's "add create", which takes the result type
            // from its first argument - the target half, which for an 8 bit movie is 8 bit, so the
            // sum of the two channels saturated at 255. On the example data that pinned 20 pixels,
            // and they were spot centres, which is exactly where the signal is. Everything
            // downstream of here - the maxima, the SNR, the prominence - was reading a clipped
            // image.
            RandomAccessibleInterval<FloatType> targetHalf = ImageJFunctions.convertFloat(images.get(0));
            RandomAccessibleInterval<FloatType> sourceHalf = ImageJFunctions.convertFloat(images.get(1));
            Shared sum = new Shared((int) targetHalf.dimension(0), (int) targetHalf.dimension(1));
            LoopBuilder.setImages(sum.img, targetHalf, sourceHalf)
                    .forEachPixel((out, t, s) -> out.set(t.get() + s.get()));

            ImagePlus sumImage = new ImagePlus("spot_qc_image", sum.processor);

            // find all spots in tne sum image.
            double[][] allSpots = getMaxima(sumImage);
            log.info("initial spot number " + allSpots.length);

            // filter spots that are near the edges of either channel.
            // overlapMask because this is where the channels overlap.
            overlapMask = createOverlapMask(averageImage.getWidth(), averageImage.getHeight());
            double[][] filteredSpots = spotFilterWithMask(allSpots, overlapMask);
            log.info("after edge proximity filter " + countGoodSpots(filteredSpots));

            // filter spots that are too close to each other.
            ImagePlus proximityMask = createSpotsNeighborhoodMask(allSpots, averageImage.getWidth()/2, averageImage.getHeight(), 2*spotSpacing);
            proximityMask = neighborhoodMaskToProximityMask(proximityMask);
            filteredSpots = spotFilterWithMask(filteredSpots, proximityMask);
            log.info("after spot proximity filter " + countGoodSpots(filteredSpots));

            // filter low SNR spots.
            backgroundMask = createSpotsNeighborhoodMask(allSpots, averageImage.getWidth()/2, averageImage.getHeight(), spotMargin);
            backgroundMask = neighborhoodMaskToBackgroundMask(backgroundMask);
            ImagePlus backgroundImage = backgroundEstimate(sumImage);
            filteredSpots = spotFilterSNR(filteredSpots, sumImage, backgroundImage);
            log.info("after SNR filter " + countGoodSpots(filteredSpots));

            // filter low prominence spots.
            filteredSpots = spotFilterProminence(filteredSpots, sumImage, backgroundImage);
            log.info("after prominence filter " + countGoodSpots(filteredSpots));

            // display as overlay on sum image.
            Overlay ov = getSpotOverlay(filteredSpots, spotMargin, Color.GREEN);
            sumImage.setOverlay(ov);

            if (!isHeadless) {
                ui.show(sumImage);
            }

            // save analysis results.
            String masksFileName = saveRootName + "_spotf_masks.tif";
            String spotsFileName = saveRootName +  "_spotf_spots.csv";

            // JSON file w/ analysis parameters, etc.
            Map<String, Object> mapping = new HashMap<>();
            mapping.put("camera black", cameraBlackLevel);
            mapping.put("camera gain", cameraGain);
            mapping.put("background kappa", clippingThreshold());
            mapping.put("edge margin", edgeMargin);
            mapping.put("end slice", endSlice);
            mapping.put("image name", inputImageName);
            mapping.put("mapping file", mappingFile);
            mapping.put("masks file", masksFileName);
            mapping.put("root name", saveRootName);
            mapping.put("spots file", spotsFileName);
            mapping.put("spot margin", spotMargin);
            mapping.put("spot prominence", spotProminence);
            mapping.put("spot sigma", spotSigma);
            mapping.put("spot spacing", spotSpacing);
            mapping.put("spot threshold", spotThreshold);
            mapping.put("spot tolerance", spotTolerance);
            mapping.put("start slice", startSlice);

            ObjectMapper mapper = new ObjectMapper();
            File saveFile = new File(saveRootName + "_spotf_finding.json");
            mapper.writeValue(saveFile, mapping);

            // Table w/ spot locations.
            saveSpotLocations(spotsFileName, filteredSpots);

            // QC image w/ identified spots.
            FileSaver qcImageSaver = new FileSaver(sumImage);
            qcImageSaver.saveAsTiff(saveRootName + "_spotf_qc_image.tif");

            // Masks that will be needed for extracting time traces.
            Concatenator cctr = new Concatenator();
            ImagePlus maskImages = cctr.concatenate(overlapMask, backgroundMask, false);
            FileSaver masksImageSaver = new FileSaver(maskImages);
            masksImageSaver.saveAsTiff(masksFileName);

	        log.info("finishing spot finding");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
