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
import ij.plugin.ImageCalculator;
import ij.plugin.filter.MaximumFinder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
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

    @Parameter (description = "radius to mask as foreground around a spot (pixels)", min = "1")
    Integer spotMargin = 4;

    @Parameter (description = "margin around the edge of the channels (pixels)", min = "1")
    Integer edgeMargin = 5;

    @Parameter (description = "background clipping threshold (robust sigmas above the estimate)", min = "0.1")
    Double backgroundKappa = 1.8;

    // Rounds of clipping in backgroundEstimate. It settles in three or four -
    // each round removes the brightest leftovers and the ones after that find
    // nothing new to remove.
    private static final int backgroundClipRounds = 4;

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
     * The default kappa of 1.8 was measured rather than chosen: on a long movie a photobleached
     * molecule must sit at zero afterwards, and 1.8 is where it does. Lower values clip too
     * hard, leave the background too low and every trace too high; higher ones keep the
     * contamination and drive traces negative, which is what the previous inpainting estimator
     * did to two thirds of them.
     *
     * It is worth re-measuring per dataset, which is why it is a parameter. How much spot light
     * reaches the trusted pixels depends on the PSF and on how crowded the field is, and the
     * right clip follows: the same movie wanted 2.5 when its spot finding had produced 488 spots
     * and 1.8 at 553. Note the circularity there - spotFilterSNR uses this estimator, so a
     * better background admits more spots, which crowds the field, which wants a tighter clip.
     * One pass either way is small; it is not worth iterating to convergence.
     */
    public ImagePlus backgroundEstimate(ImagePlus image) {
        ImageProcessor original = image.getProcessor();
        FloatProcessor pixels = original.convertToFloatProcessor();
        float[] values = (float[]) pixels.getPixels();

        boolean[] keep = trustedPixels(pixels.getWidth(), pixels.getHeight());
        float[] scratch = new float[values.length];
        double smoothing = 2.0 * (double) spotMargin;

        FloatProcessor estimate = maskedSmooth(pixels, keep, smoothing, scratch);
        for (int round = 0; round < backgroundClipRounds; round++) {
            float[] level = (float[]) estimate.getPixels();
            double spread = robustSpread(values, level, keep, scratch);
            if (spread <= 0.0) {
                break;
            }

            boolean[] fresh = new boolean[keep.length];
            boolean changed = false;
            boolean any = false;
            for (int i = 0; i < keep.length; i++) {
                fresh[i] = keep[i] && (values[i] - level[i]) < backgroundKappa * spread;
                changed |= (fresh[i] != keep[i]);
                any |= fresh[i];
            }
            if (!any || !changed) {
                break;
            }

            keep = fresh;
            estimate = maskedSmooth(pixels, keep, smoothing, scratch);
        }

        ImagePlus backgroundImage = new ImagePlus("background_image", matchType(estimate, original));
        return backgroundImage;
    }

    /**
     * Pixels the background estimate can be built from: inside the overlap and off a spot.
     */
    private boolean[] trustedPixels(int width, int height) {
        ImageProcessor overlap = overlapMask.getProcessor();
        ImageProcessor background = backgroundMask.getProcessor();

        boolean[] keep = new boolean[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                keep[y * width + x] = (overlap.get(x, y) > 0) && (background.get(x, y) > 0);
            }
        }
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
    private FloatProcessor maskedSmooth(FloatProcessor image, boolean[] keep, double sigma, float[] scratch) {
        int width = image.getWidth();
        int height = image.getHeight();
        float[] source = (float[]) image.getPixels();

        FloatProcessor weighted = new FloatProcessor(width, height);
        FloatProcessor total = new FloatProcessor(width, height);
        float[] sums = (float[]) weighted.getPixels();
        float[] counts = (float[]) total.getPixels();
        for (int i = 0; i < source.length; i++) {
            if (keep[i]) {
                sums[i] = source[i];
                counts[i] = 1.0f;
            }
        }
        weighted.blurGaussian(sigma);
        total.blurGaussian(sigma);

        // Only reached by a pixel with no kept neighbour anywhere in the kernel, which needs a
        // masked patch several times the smoothing scale across. It is a floor, not a fill.
        float floor = (float) subsetMedian(source, keep, scratch);

        FloatProcessor smoothed = new FloatProcessor(width, height);
        float[] result = (float[]) smoothed.getPixels();
        for (int i = 0; i < result.length; i++) {
            result[i] = (counts[i] > 1.0e-6f) ? (sums[i] / counts[i]) : floor;
        }
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

        // Transform and overlap allowed regions.
        java.util.List<ImagePlus> overlapImages = smfcm.splitImagePlus(overlapMask, true);
        ImagePlus overlapSumImage = ImageCalculator.run(overlapImages.get(0), overlapImages.get(1), "add create");

        // Threshold to binary.
        ImageProcessor imp = overlapSumImage.getProcessor();
        imp.threshold(190);
        imp.multiply(1.0/255.0);
        overlapSumImage.updateAndDraw();

        return overlapSumImage;
    }

    /**
     * Creates mask for (circular) neighborhood around spots.
     *
     * Spots on the pixels where this image is > 1 are filtered out.
     */
    private ImagePlus createSpotsNeighborhoodMask(double[][] spots, int imageWidth, int imageHeight, int radius){
        ImagePlus neighborhoodMask = IJ.createImage("neighborhood_mask", "16-bit black", imageWidth, imageHeight, 1);

        for (int i = 0; i < spots.length; i++) {
            int x = (int)spots[i][1];
            int y = (int)spots[i][2];
            ImagePlus temp = IJ.createImage("temp", "16-bit black", imageWidth, imageHeight, 1);
            ImageProcessor ip = temp.getProcessor();
            ip.setColor(1);
            ip.setLineWidth(0);
            ip.fillOval(x-radius, y-radius, 2*radius+1, 2*radius+1);
            temp.updateAndDraw();
            ImageCalculator.run(neighborhoodMask, temp, "add");
        }

        return neighborhoodMask;
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
        ImagePlus backgroundMask = neighborhoodMask.duplicate();
        backgroundMask.setTitle("foreground_mask");
        ImageProcessor imp = backgroundMask.getProcessor();

        for (int i = 0; i < backgroundMask.getWidth(); i++){
            for (int j = 0; j < backgroundMask.getHeight(); j++) {
                if (backgroundMask.getPixel(i,j)[0] > 0){
                    imp.putPixel(i,j,0);
                }
                else{
                    imp.putPixel(i,j,1);
                }
            }
        }
        return backgroundMask;
    }

    /**
     * converts neighborhood mask for use to detect spots that are too close to each other.
     */
    private ImagePlus neighborhoodMaskToProximityMask(ImagePlus neighborhoodMask) {
        ImagePlus proximityMask = neighborhoodMask.duplicate();
        proximityMask.setTitle("overlap_mask");
        ImageProcessor imp = proximityMask.getProcessor();

        for (int i = 0; i < proximityMask.getWidth(); i++){
            for (int j = 0; j < proximityMask.getHeight(); j++) {
                if (proximityMask.getPixel(i,j)[0] > 1){
                    imp.putPixel(i,j,0);
                }
                else{
                    imp.putPixel(i,j,1);
                }
            }
        }
        return proximityMask;
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
        ImagePlus fgImage = ImageCalculator.run(sumImage, backgroundImage, "subtract create");
        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            // Calculate spot lowest prominence over pixels in circular neighborhood.
            double spotHeight = fgImage.getPixel(x,y)[0];
            double lowestProminence = spotHeight;
            for (int rx = -srad; rx <= srad; rx += 1){
                for (int ry = -srad; ry <= srad; ry += 1){
                    if (rx*rx+ry*ry >= (srad-1)*(srad-1) && rx*rx+ry*ry <= (srad+1)*(srad+1)){
                        double pixelHeight = Math.max(1, fgImage.getPixel(x+rx,y+ry)[0]);
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

        ImagePlus fgSmooth = ImageCalculator.run(sumImage, backgroundImage, "subtract create");
        ImageProcessor impFg = fgSmooth.getProcessor();
        impFg.blurGaussian(spotSigma);
        fgSmooth.setTitle("foreground_smooth");

        ImagePlus bgSmooth = backgroundImage.duplicate();
        ImageProcessor impBg = bgSmooth.getProcessor();
        impBg.blurGaussian(spotSigma);
        bgSmooth.setTitle("background_smooth");

        // In order to calculate the SNR integrated over the spot we need the integrated magnitude
        // of the spot. We multiply by norm because blurGaussian() uses a normalized Gaussian when
        // for our purposes a unit height Gaussian would have been the correct thing to use.
        double norm = 2.0*Math.PI*spotSigma*spotSigma;
        for (int i = 0; i < spots.length; i++) {
            System.arraycopy(spots[i], 0, filteredSpots[i], 0, spots[i].length);

            int x = (int)spots[i][1];
            int y = (int)spots[i][2];

            double fg = norm*cameraGain*(fgSmooth.getPixel(x,y)[0]);
            double bg = norm*cameraGain*(bgSmooth.getPixel(x,y)[0] - 2*cameraBlackLevel);

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
            FileSaver fgSmoothImageSaver = new FileSaver(fgSmooth);
            fgSmoothImageSaver.saveAsTiff(saveRootName + "_spotf_fg_smooth.tif");

            FileSaver bgSmoothImageSaver = new FileSaver(bgSmooth);
            bgSmoothImageSaver.saveAsTiff(saveRootName + "_spotf_bg_smooth.tif");
        }

        return filteredSpots;
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

            // Load the channel to channel mapping file.
            loadMappingJSON(mappingFile.toString());

            // Average image.
            log.info("average image - " + inputImage.getNSlices() + " slices");
            ImagePlus averageImage = smfcm.averageImagePlus(inputImage, startSlice, endSlice);

            // split, transform and add the two channels together.
            log.info("split and transform");
            java.util.List<ImagePlus> images = smfcm.splitImagePlus(averageImage, true);

            ImagePlus sumImage = ImageCalculator.run(images.get(0), images.get(1), "add create");
            sumImage.setTitle("spot_qc_image");

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
            mapping.put("background kappa", backgroundKappa);
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
