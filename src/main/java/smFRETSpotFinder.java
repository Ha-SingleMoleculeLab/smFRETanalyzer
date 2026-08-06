/*
 * This class finds the single molecule spots using the mapping.
 */

import com.fasterxml.jackson.databind.ObjectMapper;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
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
import org.scijava.command.Interactive;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;
import org.scijava.widget.Button;

import java.awt.*;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


/**
 * Interactive, so the parameter dialog stays open between runs.
 *
 * Spot finding is a thing you tune rather than a thing you do once - the threshold, the tolerance
 * and the prominence all have to be looked at against the image before they are right - and with a
 * plain Command that meant going back to the Plugins menu for every adjustment. Implementing
 * Interactive gets a non-modal dialog out of the same @Parameter declarations, so nothing about
 * the widgets, the labels or the min validation changes; only its lifetime does. The findSpots
 * button is what actually runs an analysis, because the alternative - running on every keystroke -
 * would re-analyse and rewrite four output files while a number is being typed. Opening the
 * dialog does nothing either - see run() for how that is kept out of the way of the macros.
 */
@Plugin(type = Command.class, headless = true,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Spot Finder", weight = 2.0)})
public class smFRETSpotFinder implements Command, Interactive {
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
    Double spotTolerance = 5.0;

    @Parameter (description = "spot prominence, as a fraction of an isolated spot", min = "0.0")
    Double spotProminence = 0.4;

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

    // Half width of the prominence annulus, pixels. Wide enough that the ring is
    // never empty at the smallest spot size - it holds 16 pixels at spotSigma 1 -
    // and narrow enough that it stays out at two sigma rather than reaching back
    // in toward the spot, which is what the old ring did.
    private static final double promRingHalfWidth = 0.6;

    // The 90th percentile of a standard normal. The ring is summarised by its 90th
    // percentile and its pixels carry zero mean noise, so the summary sits this
    // many noise sigmas high and the bias is taken back off - see
    // spotFilterProminence for what that is worth.
    private static final double promNoiseBias = 1.2816;

    // Reported prominence is capped here rather than running away when the ring is
    // empty, which would otherwise put an infinity in the spots CSV.
    private static final double promCap = 10.0;

    // Radius masked as foreground around each spot, in pixels: it follows the
    // spot size, as max(2, round(2.5 * spotSigma)).
    //
    // It was a fixed 4, justified by a sweep that found anything from 2 to 6
    // within 11% at any spot size. That sweep ran against the estimator before
    // the one-sided clip's downward bias was corrected - see truncationConstants
    // - and that bias was partly cancelling the contamination this mask exists
    // to stop, so the mask looked less important than it is. It also ran on a
    // sparser field than real data.
    //
    // Re-swept from 2 to 10 pixels at spot sizes 1.0 to 3.0, at two densities,
    // scored on the bias of the background estimate at the spot locations since
    // that is what gets subtracted from a trace. The best margin tracks the spot
    // size cleanly, unlike backgroundKappa: at 400 spots in a channel it runs
    // 2, 3, 5, 8 and beyond 10 as sigma goes 1.0 to 3.0.
    //
    // 2.5 sigma rather than the 3 sigma a PSF actually reaches, because masking
    // is not free - every masked pixel is one the estimator cannot use, and at
    // 400 spots a radius of 8 already covers most of the field. That shows up as
    // a genuine turning point: at sigma 3 the error bottoms at a margin of 8
    // (rms 0.278) and gets worse by 10 (0.288) even though the bias is still
    // falling. Against the fixed 4, this cuts the bias over the sigma range from
    // 0.225 to 0.098 ADU at 400 spots and from 0.066 to 0.030 at 150 - worth
    // roughly 15 units of trace at sigma 3, systematically.
    //
    // A scale of 2.0 was the alternative and is marginally better on the sparse
    // field, but half as good on the crowded one, and crowded is what the real
    // data is. The floor of 2 only matters below spotSigma 0.6.
    private static final double spotMarginScale = 2.5;
    private static final int spotMarginFloor = 2;

    private int spotMargin() {
        return Math.max(spotMarginFloor, (int) Math.round(spotMarginScale * spotSigma));
    }

    @Parameter (description = "background clipping threshold, 0 to derive it from the spot size", min = "0.0")
    Double backgroundKappa = 0.0;

    // Which image spots are found in. "sum" is what this always did, and stays the default.
    //
    // Short values rather than "donor/target", so that a macro can say spotchannel=donor without
    // having to bracket-quote it.
    @Parameter (description = "channel to find spots in: the donor/target half, the acceptor/source"
                              + " half, or the two added together",
                label = "Spot channel", choices = {"sum", "donor", "acceptor"})
    String spotChannel = "sum";

    @Parameter (label = "Find spots", description = "run spot finding with the settings above",
                callback = "findSpots")
    private Button findSpotsButton;

    // Rounds of clipping in backgroundEstimate. It settles in three or four -
    // each round removes the brightest leftovers and the ones after that find
    // nothing new to remove.
    private static final int backgroundClipRounds = 4;

    // Default clipping threshold, used when backgroundKappa is left at zero. A
    // constant, and re-derived against the current estimator by sweeping kappa
    // from 0.4 to 10 at spot sizes 1.0 to 3.0, against a known background under
    // a Gaussian beam, at the density of the real data - 400 spots in a channel,
    // randomly placed, log-normal brightness. Scored on the bias of the estimate
    // at the spot locations, since that is what gets subtracted from a trace.
    //
    // This used to be 1.5 - 0.3 * spotSigma, floored at 0.5. That line was fitted
    // when the one-sided clip still biased the level downwards - see
    // truncationConstants - so it was partly compensating for the estimator
    // rather than describing the spots. With the bias corrected, the best kappa
    // is *not* monotone in spot size: it runs loose at sigma 1 to 1.5, tightest
    // near sigma 2, and back to about 1.0 to 1.3 by sigma 3, which no line fits.
    // Aggregated over spot size the optimum is flat from about 1.0 to 1.6, and
    // 1.3 is the best single value with or without sigma 3 included.
    //
    // Do not read the residual error at large spot sizes as a kappa problem. At
    // sigma 3 the estimate read 0.44 ADU high at every kappa, and thinning the
    // field from 400 spots to 50 took that to 0.03 - the clip was not what was
    // failing there, spotMargin was, being a fixed 4 pixels where the PSF
    // reaches three times that. spotMargin follows spotSigma now, which takes
    // that residual down but does not abolish it.
    private static final double backgroundKappaDefault = 1.3;

    // Smoothing scale for the background estimate, in pixels. A constant, and
    // measured to be one: the best value sat at 12 to 16 pixels at every spot
    // size tried, with no trend against it (the fit is 14.4 - 0.00 * sigma).
    // It is a property of how fast the illumination varies, which is a fact
    // about the microscope rather than about the spots.
    //
    // This used to be 2 * spotMargin, which at the default spot size gave 8 and
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

    // Everything up to the summed image - loading, the temporal average, the split and the warp -
    // depends only on these four inputs and on none of the parameters being tuned, so it is kept
    // between runs and rebuilt only when one of them changes. It is most of the cost of a run,
    // which is the difference between adjusting a threshold and waiting for one.
    //
    // startSlice and endSlice are in the key because they choose which frames are averaged: change
    // either and the summed image is a different image, so the whole prefix has to be redone.
    private String cachedInputs;
    private Shared cachedTarget;
    private Shared cachedSource;
    private int cachedWidth;
    private int cachedHeight;

    // How many frames actually went into the average, after averageImagePlus has clamped the
    // requested range to the movie. The prominence filter needs it to know how far the temporal
    // average has taken the noise down.
    private int cachedFrames;

    // The displayed QC image, kept so that repeated runs update one window instead of opening one
    // window per run.
    private ImagePlus qcWindow;

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
     *     the mask to do. That held while the clip was biased; with the bias corrected the
     *     mask matters again at large spot sizes, which is why spotMargin follows spotSigma;
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
     * A constant now, not a function of spotSigma - see backgroundKappaDefault for the sweep
     * that says so. It used to be derived from the spot size on the reasoning that the clip has
     * to catch the wings of whatever PSF it is given; that is true, but it is not what the old
     * line was measuring, because it was fitted while the clip still biased the level downwards
     * and so was trading contamination against a bias that no longer exists.
     *
     * Setting backgroundKappa above zero still overrides it. The escape hatch stays because the
     * sweep behind the default assumes a smooth Gaussian beam, and a real illumination profile
     * can have structure at scales a hard clip would eat. The example long movie wanted 1.8,
     * which is much closer to this constant than to the 0.93 the old line gave it.
     */
    public double clippingThreshold() {
        if (backgroundKappa > 0.0) {
            return backgroundKappa;
        }
        return backgroundKappaDefault;
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
     * Mark spots whose surroundings are too bright for an isolated molecule.
     *
     * The spot's height is compared against the level on a ring at two sigma, and the ratio is
     * divided by what a perfect isolated spot of this sigma would give, so **the number reads as
     * a fraction of ideal: 1.0 is a clean lone spot whatever its size or brightness.**
     *
     * All three parts of that were measured rather than chosen. The ring was
     * (srad-1)^2 <= r^2 <= (srad+1)^2 with srad = round(2*spotSigma), whose inner edge is 1 px
     * at sigma 1 but 5 px at sigma 3, so a lone spot scored exp(rmin^2/2 sigma^2) - 1.649 at
     * sigma 1 rising to 4.010 at sigma 3. A fixed threshold therefore meant something different
     * at every spot size, and the old default of 1.8 sat *above* the sigma 1 ideal: at
     * spotSigma 1.0 it rejected every real spot. Hence a thin annulus at 2 sigma and the
     * explicit normalization.
     *
     * The ring was also summarised by its brightest pixel, which is an extreme value statistic
     * over 28 to 80 pixels and so mostly noise. The 90th percentile cut the best achievable
     * error rate sum from 0.830 to 0.685 averaged over sigma 1 to 3. Doing this on the smoothed
     * foreground was tried and is worse (0.776) - smoothing blurs the neighbour into the spot,
     * which is the thing being looked for.
     *
     * The old max(1, ...) on each ring pixel was an absolute ADU floor, and it is what made the
     * filter select against bright spots on real data: a dim spot's ring is noise near zero,
     * clamped to 1, scoring a huge ratio, where a bright spot's ring is real PSF wings scoring
     * their true ratio. On hel1 that rejected 52% of spots at SNR 15-25 against 16% below 10.
     * A ring at or below zero is now simply taken as "nothing there, so no neighbour".
     */
    private double[][] spotFilterProminence(double[][] spots, ImagePlus sumImage, ImagePlus backgroundImage) {
        double[][] filteredSpots = new double[spots.length][spots[0].length+1];
        int last_col = filteredSpots[0].length-1;

        Shared fgImage = foreground(sumImage, backgroundImage);
        ImageProcessor bgProc = backgroundImage.getProcessor();
        int[][] ring = prominenceRing();
        double ideal = idealProminence(ring);
        double frames = Math.max(1, cachedFrames);
        float[] ringValues = new float[ring.length];

        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            double spotHeight = valueAt(fgImage, x, y);
            for (int j = 0; j < ring.length; j++) {
                ringValues[j] = valueAt(fgImage, x + ring[j][0], y + ring[j][1]);
            }

            // The ring pixels have had the background taken off, so their noise is zero mean and
            // the 90th percentile of it sits promNoiseBias sigmas above the true ring level. Left
            // in, that reads the ring high and the prominence low, and by more at larger spot
            // sizes where the ring is a smaller fraction of the peak - which is what stopped the
            // normalization above from actually equalizing the scores. Taking it off moves the
            // median lone spot from 0.90/0.64/0.70/0.66/0.47 at sigma 1.0 to 3.0 to
            // 1.02/0.92/1.07/1.04/0.77, which is the 1.0 it is supposed to be.
            double levelADU = Math.max(0.0, bgProc.getf(x, y) - channelCount()*cameraBlackLevel);
            double noise = Math.sqrt(levelADU / (cameraGain * frames));
            double ringLevel = percentile90(ringValues) - promNoiseBias*noise;

            double prominence = promCap;
            if (ringLevel > 0.0) {
                prominence = Math.min(promCap, (spotHeight / ringLevel) / ideal);
            }

            filteredSpots[i][last_col] = prominence;
            if (prominence < spotProminence){
                filteredSpots[i][0] = 0.0;
            }
        }

        return filteredSpots;
    }

    /**
     * Pixel offsets of the prominence ring: a thin annulus at two sigma.
     */
    private int[][] prominenceRing() {
        double radius = 2.0*spotSigma;
        int reach = (int) Math.ceil(radius + promRingHalfWidth);

        java.util.List<int[]> offsets = new java.util.ArrayList<>();
        for (int rx = -reach; rx <= reach; rx++) {
            for (int ry = -reach; ry <= reach; ry++) {
                if (Math.abs(Math.hypot(rx, ry) - radius) <= promRingHalfWidth) {
                    offsets.add(new int[] {rx, ry});
                }
            }
        }
        return offsets.toArray(new int[0][]);
    }

    /**
     * What the height over ring ratio comes to for a noiseless isolated spot of this sigma.
     *
     * Computed from the ring's own offsets with the same summary statistic the measurement uses,
     * so the two cancel exactly and a perfect lone spot scores 1.0 - rather than being a constant
     * that would have to be kept in step with the ring geometry by hand.
     */
    private double idealProminence(int[][] ring) {
        float[] profile = new float[ring.length];
        for (int i = 0; i < ring.length; i++) {
            double r2 = ring[i][0]*ring[i][0] + ring[i][1]*ring[i][1];
            profile[i] = (float) Math.exp(-r2 / (2.0*spotSigma*spotSigma));
        }
        return 1.0 / percentile90(profile);
    }

    /**
     * The 90th percentile, interpolated between the two neighbouring order statistics.
     *
     * A plain sort rather than the quickselect used for the background medians: the ring holds
     * 16 to 40 values and this runs once per spot rather than once per pixel, so the selection
     * would save nothing, and interpolating wants both neighbours rather than one.
     */
    private static double percentile90(float[] values) {
        float[] sorted = values.clone();
        java.util.Arrays.sort(sorted);

        double position = 0.9*(sorted.length - 1);
        int lower = (int) Math.floor(position);
        int upper = Math.min(lower + 1, sorted.length - 1);
        return sorted[lower] + (position - lower)*(sorted[upper] - sorted[lower]);
    }

    /**
     * Mark spots with estimated SNR less than threshold.
     *
     * The camera black level is multiplied by how many channels went into the image spots are
     * being found in - two for the sum, one for a single channel - since that is how many black
     * levels the background carries.
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
            double bg = norm*cameraGain*(valueAt(bgSmooth, x, y) - channelCount()*cameraBlackLevel);

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
    /**
     * Analyse straight away when there is no dialog to press a button in, and not otherwise.
     *
     * SciJava runs an Interactive command once as it opens the dialog, which is not what opening
     * a dialog should mean: the point of keeping it open is to set the parameters up before
     * anything happens. So an interactive session waits for findSpots.
     *
     * A macro has no button to press, and the batch macros drive this as
     * run("smFRET Spot Finder", "inputimagename=... spotthreshold=..."), so a run carrying macro
     * options has to analyse here or it would do nothing at all. ij.Macro.getOptions is non-null
     * exactly for that case. Headless is the same argument with no UI at all.
     */
    public void run() {
        if (isHeadless || (Macro.getOptions() != null)) {
            findSpots();
        }
    }

    /**
     * Rebuild the summed image, unless the inputs it depends on are unchanged since last time.
     *
     * Sets cachedTarget, cachedSource, cachedWidth and cachedHeight. The mapping is loaded here too
     * rather than on every run, which is safe because the mapping file is part of the key.
     *
     * The two halves are cached rather than the image spots are found in, so that changing
     * spotChannel does not throw away the expensive part. Comparing the same field across the
     * donor, the acceptor and the sum is exactly the thing an interactive dialog is for.
     */
    private void buildChannelImages() {
        String inputs = inputImageName + "|" + mappingFile + "|" + startSlice + "|" + endSlice;
        if (inputs.equals(cachedInputs) && (cachedTarget != null)) {
            log.info("reusing the averaged image - inputs unchanged");
            return;
        }
        cachedInputs = null;

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

        // split and transform the source half onto the target half's coordinate frame.
        log.info("split and transform");
        java.util.List<ImagePlus> images = smfcm.splitImagePlus(averageImage, true);

        cachedTarget = copyOf(ImageJFunctions.convertFloat(images.get(0)));
        cachedSource = copyOf(ImageJFunctions.convertFloat(images.get(1)));
        cachedWidth = averageImage.getWidth();
        cachedHeight = averageImage.getHeight();

        // The same clamping averageImagePlus applies, so a range that runs off the end of the
        // movie does not leave the noise estimate thinking it averaged more than it did.
        int frames = smFRETChannelMapper.frameCount(inputImage);
        cachedFrames = (frames == 1) ? 1
                : (Math.min(endSlice, frames) - Math.max(1, startSlice) + 1);
        cachedInputs = inputs;
    }

    /**
     * A Shared holding a copy of the pixels, so the cached halves are never written through.
     */
    private static Shared copyOf(RandomAccessibleInterval<FloatType> image) {
        Shared out = new Shared((int) image.dimension(0), (int) image.dimension(1));
        int[] index = {0};
        LoopBuilder.setImages(Views.flatIterable(image))
                .forEachPixel(v -> out.pixels[index[0]++] = v.get());
        return out;
    }

    /**
     * The image spots are found in, built from the cached halves according to spotChannel.
     *
     * The sum is computed in float. It was ImageCalculator's "add create", which takes the result
     * type from its first argument - the target half, which for an 8 bit movie is 8 bit, so the
     * sum of the two channels saturated at 255. On the example data that pinned 20 pixels, and
     * they were spot centres, which is exactly where the signal is. Everything downstream of here
     * - the maxima, the SNR, the prominence - was reading a clipped image.
     *
     * A fresh image every run, including for the single channel cases where it is only a copy.
     * Handing back the cached half itself would be quicker by a fraction of a millisecond and
     * would let anything downstream that ever wrote to it corrupt the cache.
     */
    private Shared channelImage() {
        Shared out = new Shared(cachedTarget.width, cachedTarget.height);
        if ("donor".equals(spotChannel)) {
            System.arraycopy(cachedTarget.pixels, 0, out.pixels, 0, out.pixels.length);
        } else if ("acceptor".equals(spotChannel)) {
            System.arraycopy(cachedSource.pixels, 0, out.pixels, 0, out.pixels.length);
        } else {
            LoopBuilder.setImages(out.img, cachedTarget.img, cachedSource.img)
                    .forEachPixel((o, t, s) -> o.set(t.get() + s.get()));
        }
        return out;
    }

    /**
     * How many camera black levels are in the image spots are found in.
     *
     * Two when the halves are added together, one otherwise. spotFilterSNR subtracts this from the
     * background before taking its square root, so getting it wrong on a single channel would
     * over-subtract a whole black level - 5 ADU on the example data - and shift every SNR.
     */
    private int channelCount() {
        return "sum".equals(spotChannel) ? 2 : 1;
    }

    /**
     * Find spots with the current settings, and write the results.
     *
     * This is both what the dialog's button calls and what run() does, so a scripted single run
     * and an interactive adjustment go down exactly the same path.
     *
     * It returns quietly when the two files have not been chosen yet, because an interactive
     * dialog exists before its inputs do and the framework is free to call run() at that point.
     */
    public void findSpots() {
        try {
            if ((inputImageName == null) || (mappingFile == null)) {
                log.info("choose an image and a mapping file, then press Find spots");
                return;
            }
	        log.info("starting spot finding on image " + inputImageName);

            // Root name to use for saving output, this is just the file name
            // without the extension.
            saveRootName = inputImageName.toString();
            int dotIndex = saveRootName.lastIndexOf('.');
            if (dotIndex > 0) {
                saveRootName = saveRootName.substring(0, dotIndex);
            }
            log.info("save root " + saveRootName);

            buildChannelImages();
            Shared sum = channelImage();
            log.info("finding spots in the " + spotChannel + " image");
            ImagePlus sumImage = new ImagePlus("spot_qc_image", sum.processor);

            // find all spots in tne sum image.
            double[][] allSpots = getMaxima(sumImage);
            log.info("initial spot number " + allSpots.length);

            // filter spots that are near the edges of either channel.
            // overlapMask because this is where the channels overlap.
            overlapMask = createOverlapMask(cachedWidth, cachedHeight);
            double[][] filteredSpots = spotFilterWithMask(allSpots, overlapMask);
            log.info("after edge proximity filter " + countGoodSpots(filteredSpots));

            // filter spots that are too close to each other.
            ImagePlus proximityMask = createSpotsNeighborhoodMask(allSpots, cachedWidth/2, cachedHeight, 2*spotSpacing);
            proximityMask = neighborhoodMaskToProximityMask(proximityMask);
            filteredSpots = spotFilterWithMask(filteredSpots, proximityMask);
            log.info("after spot proximity filter " + countGoodSpots(filteredSpots));

            // filter low SNR spots.
            backgroundMask = createSpotsNeighborhoodMask(allSpots, cachedWidth/2, cachedHeight, spotMargin());
            backgroundMask = neighborhoodMaskToBackgroundMask(backgroundMask);
            ImagePlus backgroundImage = backgroundEstimate(sumImage);
            filteredSpots = spotFilterSNR(filteredSpots, sumImage, backgroundImage);
            log.info("after SNR filter " + countGoodSpots(filteredSpots));

            // filter low prominence spots.
            filteredSpots = spotFilterProminence(filteredSpots, sumImage, backgroundImage);
            log.info("after prominence filter " + countGoodSpots(filteredSpots));

            // display as overlay on sum image.
            Overlay ov = getSpotOverlay(filteredSpots, spotMargin(), Color.GREEN);
            sumImage.setOverlay(ov);

            // One window, updated in place. Showing a new one per run would leave a window per
            // parameter adjustment, which is the opposite of what an interactive dialog is for.
            if (!isHeadless) {
                if ((qcWindow != null) && (qcWindow.getWindow() != null)) {
                    qcWindow.setProcessor(sum.processor);
                    qcWindow.setOverlay(ov);
                    qcWindow.updateAndDraw();
                } else {
                    qcWindow = sumImage;
                    ui.show(sumImage);
                }
            }

            // save analysis results.
            String masksFileName = saveRootName + "_spotf_masks.tif";
            String spotsFileName = saveRootName +  "_spotf_spots.csv";

            // JSON file w/ analysis parameters, etc.
            Map<String, Object> mapping = new HashMap<>();
            mapping.put("camera black", cameraBlackLevel);
            mapping.put("spot channel", spotChannel);
            mapping.put("camera gain", cameraGain);
            mapping.put("background kappa", clippingThreshold());
            mapping.put("edge margin", edgeMargin);
            mapping.put("end slice", endSlice);
            mapping.put("image name", inputImageName);
            mapping.put("mapping file", mappingFile);
            mapping.put("masks file", masksFileName);
            mapping.put("root name", saveRootName);
            mapping.put("spots file", spotsFileName);
            mapping.put("spot margin", spotMargin());
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
