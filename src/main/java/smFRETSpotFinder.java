/*
 * This class finds the single molecule spots using the mapping.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.MaximumFinder;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;


@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Spot Finder")
public class smFRETSpotFinder implements Command {
    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter
    ImagePlus inputImage;

    @Parameter (description = "first slice for averaging", min = "1")
    Integer startSlice = 1;

    @Parameter (description = "last slice for averaging", min = "1")
    Integer endSlice = 30;

    @Parameter (description = "spot size detection threshold", min = "1.0")
    Double spotThreshold = 8.0;

    @Parameter (description = "spot size (sigma)", min = "0.2")
    Double spotSigma = 1.5;

    @Parameter (description = "camera offset / black level", min = "1")
    Integer cameraBlackLevel = 1;

    @Parameter (description = "camera gain (e-/ADU)", min = "0.1")
    Double cameraGain = 1.0;

    @Parameter(description = "channel to channel mapping file", label = "mapping file", style = "open")
    File mappingFile;

    @Parameter(description = "directory to save results in", label = "save directory", style = "directory")
    File saveDirectory;

    @Parameter (description = "minimum allowed distance between spots", min = "1")
    Integer spotSpacing = 5;

    @Parameter (description = "margin around the edge of the channels", min = "1")
    Integer edgeMargin = 5;

    // Member variables.
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();

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
        java.util.List<ImagePlus> overlapImages = this.smfcm.splitImagePlus(overlapMask, true);
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
    private ImagePlus createSpotsNeighborhoodMask(Polygon spots, int imageWidth, int imageHeight, int radius){
        ImagePlus neighborhoodMask = IJ.createImage("neighborhood_mask", "16-bit black", imageWidth, imageHeight, 1);

        for (int i = 0; i < spots.npoints; i++) {
            int x = spots.xpoints[i];
            int y = spots.ypoints[i];
            ImagePlus temp = IJ.createImage("temp", "16-bit black", imageWidth, imageHeight, 1);
            ImageProcessor ip = temp.getProcessor();
            ip.setColor(1);
            ip.setLineWidth(0);
            ip.fillOval(x-radius, y-radius, 2*radius, 2*radius);
            temp.updateAndDraw();
            ImageCalculator.run(neighborhoodMask, temp, "add");
        }

        return neighborhoodMask;
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
    public static Overlay getSpotOverlay (Polygon spots, int radius, Color symbolColor) {
        Overlay ov = new Overlay();
        int diameter = 2 * radius;
        for (int i = 0; i < spots.npoints; i++) {
            int x = spots.xpoints[i];
            int y = spots.ypoints[i];
            Roi roi = new Roi(x - radius, y - radius,
                    diameter, diameter, diameter);
            roi.setStrokeColor(symbolColor);
            ov.add(roi);
        }
        return ov;
    }

    /**
     * Load an existing mapping file to initialize smFRETChannelMapper.
     */
    private smFRETChannelMapper smfcm = new smFRETChannelMapper();  // smFREChannelMapper object.
    public void loadMappingJSON(File mappingFileName){
        this.smfcm.loadMappingJSON(mappingFileName);
        this.smfcm.log = log;
    }

    /**
     * Fills in masked values by repeated Gaussian convolution.
     */
    private ImagePlus maskInpaint(ImagePlus image, ImagePlus mask, double sigma){
        ImagePlus filledImage = image.duplicate();
        ImagePlus lastFilledImage = image.duplicate();

        ImageProcessor imp = filledImage.getProcessor();
        for (int i = 0; i < 400; i++){
            imp.blurGaussian(sigma);
            if (maskInpaintDifference(filledImage, lastFilledImage, mask) < 1){
                break;
            }
            maskInpaintReset(image, filledImage, mask);
            lastFilledImage = image.duplicate();
        }
        return filledImage;
    }

    /**
     * maskInpaint helper function.
     */
    private int maskInpaintDifference(ImagePlus image, ImagePlus lastImage, ImagePlus mask) {
        int diff = 0;
        for (int i = 0; i < image.getWidth(); i++){
            for (int j = 0; j < image.getHeight(); j++) {
                if (mask.getPixel(i,j)[0] == 0){
                    int tmp = Math.abs(image.getPixel(i,j)[0] - lastImage.getPixel(i,j)[0]);
                    if (tmp > diff){
                        diff = tmp;
                    }
                }
            }
        }
        return diff;
    }

    /**
     * maskInpaint helper function.
     */
    private void maskInpaintReset(ImagePlus originalImage, ImagePlus modifiedImage, ImagePlus mask) {
        ImageProcessor imp = modifiedImage.getProcessor();

        for (int i = 0; i < originalImage.getWidth(); i++){
            for (int j = 0; j < originalImage.getHeight(); j++) {
                if (mask.getPixel(i,j)[0] == 1){
                    imp.putPixel(i,j, originalImage.getPixel(i,j)[0]);
                }
            }
        }
    }

    /**
     * Returns the minimum value in the region where overlayMask is above threshold.
     */
    /*
    private int overlayMinValue(ImagePlus image, ImagePlus mask) {
        int ovMin = 65000;
        for (int i = 0; i < image.getWidth(); i++) {
            for (int j = 0; j < image.getWidth(); j++) {
                if (mask.getPixel(i,j)[0] > smFRETSpotFinder.overlapThreshold){
                    int pixValue = image.getPixel(i,j)[0];
                    if (pixValue < ovMin) {
                        ovMin = pixValue;
                    }
                }
            }
        }
        return ovMin;
    }
     */

    /**
     *
     */
    private Polygon spotFilterWithMask(Polygon spots, ImagePlus mask, int threshold, int sign){
        Polygon filteredSpots = new Polygon();

        for (int i = 0; i < spots.npoints; i++) {
            int x = spots.xpoints[i];
            int y = spots.ypoints[i];

            if (mask.getPixel(x,y)[0]*sign > threshold*sign){
                filteredSpots.addPoint(x,y);
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
	        log.info("starting spot finding");
            loadMappingJSON(mappingFile);

            // Average image.
            log.info("average image - " + inputImage.getNSlices() + " slices");
            ImagePlus averageImage = this.smfcm.averageImagePlus(inputImage, startSlice, endSlice);

            // split, transform and add the two channels together.
            log.info("split and transform");
            java.util.List<ImagePlus> images = this.smfcm.splitImagePlus(averageImage, true);

            ImagePlus sumImage = ImageCalculator.run(images.get(0), images.get(1), "add create");
            sumImage.setTitle("channels_sum_image");

            // find all spots in tne sum image.
            ImageProcessor sumImageIp = sumImage.getProcessor();
            MaximumFinder mf = new MaximumFinder();
            Polygon allSpots = mf.getMaxima(sumImageIp, 5, true);

            // filter spots that are near the edges of either channel.
            ImagePlus overlapMask = createOverlapMask(averageImage.getWidth(), averageImage.getHeight());
            Polygon filteredSpots = spotFilterWithMask(allSpots, overlapMask, 0, 1);

            // filter spots that are too close to each other.
            ImagePlus neighborhoodMask = createSpotsNeighborhoodMask(allSpots, averageImage.getWidth()/2, averageImage.getHeight(), 2*spotSpacing);
            filteredSpots = spotFilterWithMask(filteredSpots, neighborhoodMask, 2, -1);

            // filter low SNR spots.
            ImagePlus foregroundMask = createSpotsNeighborhoodMask(allSpots, averageImage.getWidth()/2, averageImage.getHeight(), spotSpacing);

            // display as overlay on sum image.
            sumImage = maskInpaint(sumImage, overlapMask, 2.0);

            Overlay ov = getSpotOverlay(filteredSpots, 5, Color.GREEN);
            sumImage.setOverlay(ov);

            if (!this.isHeadless) {
                ui.show(overlapMask);
                ui.show(sumImage);
            }

            FileSaver sourceFile = new FileSaver(sumImage);
            sourceFile.saveAsTiff(saveDirectory.getAbsolutePath() + File.separator + "sum_image.tif");

            sourceFile = new FileSaver(foregroundMask);
            sourceFile.saveAsTiff(saveDirectory.getAbsolutePath() + File.separator + "mask_image.tif");

	        log.info("finishing spot finding");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
