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
import ij.process.ShortProcessor;

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
    private static final double backgroundSmoothing = 14.0;

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
        ImageProcessor original = image.getProcessor();

        // duplicate, because convertToFloatProcessor hands back the same object when the input is
        // already float, and Shared would then be aliasing the caller's pixels.
        Shared source = new Shared((FloatProcessor) original.convertToFloatProcessor().duplicate());
        float[] values = source.pixels;

        boolean[] keep = trustedPixels(source.width, source.height);
        float[] scratch = new float[values.length];
        double smoothing = backgroundSmoothing;
        double kappa = clippingThreshold();

        Shared estimate = maskedSmooth(source, keep, smoothing, scratch);
        for (int round = 0; round < backgroundClipRounds; round++) {
            float[] level = estimate.pixels;
            double spread = robustSpread(values, level, keep, scratch);
            if (spread <= 0.0) {
                break;
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
            estimate = maskedSmooth(source, keep, smoothing, scratch);
        }

        ImagePlus backgroundImage = new ImagePlus("background_image", matchType(estimate.processor, original));
        return backgroundImage;
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
        gauss(sigma, Views.extendZero(weighted.img), weightedSmooth.img);
        gauss(sigma, Views.extendZero(total.img), totalSmooth.img);

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
     * Median of the first count entries, which are sorted in place.
     */
    private static double median(float[] values, int count) {
        Arrays.sort(values, 0, count);
        int middle = count / 2;
        if ((count % 2) == 0) {
            return 0.5 * ((double) values[middle - 1] + (double) values[middle]);
        }
        return values[middle];
    }

    /**
     * Put the estimate back into the type it was measured from, so that this stays a change of
     * algorithm and not of precision. Writing it as float would drop the ~0.3 ADU of rounding
     * noise the stored background carries, but that is a separate question with its own answer.
     */
    private static ImageProcessor matchType(FloatProcessor estimate, ImageProcessor template) {
        if (template instanceof FloatProcessor) {
            return estimate;
        }
        if (template instanceof ShortProcessor) {
            return estimate.convertToShort(false);
        }
        return estimate.convertToByte(false);
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

        return new ImagePlus("overlap_mask", overlapSum.processor.convertToShort(false));
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

        return new ImagePlus("neighborhood_mask", neighborhoodMask.processor.convertToShort(false));
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

        ImagePlus mask = new ImagePlus(title, inverted.processor.convertToShort(false));
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

        // In order to calculate the SNR integrated over the spot we need the integrated magnitude
        // of the spot. We multiply by norm because blurGaussian() uses a normalized Gaussian when
        // for our purposes a unit height Gaussian would have been the correct thing to use.
        double norm = 2.0*Math.PI*spotSigma*spotSigma;
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
            if (mask.getPixel(x,y)[0] == 0){
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

            // Load the channel to channel mapping file.
            loadMappingJSON(mappingFile.toString());

            // Average image.
            log.info("average image - " + inputImage.getNSlices() + " slices");
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
