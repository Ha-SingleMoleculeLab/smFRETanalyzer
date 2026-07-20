/*
 * This class handles everything related to mapping the two channels to the same space.
 * It depends on the TurboReg plugin (version 2.0.1).
 */

import com.fasterxml.jackson.databind.ObjectMapper;

import java.awt.*;
import java.io.File;
import java.lang.reflect.Method;
import java.util.*;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ImageConverter;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;


@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Channel Mapping")
public class smFRETChannelMapper implements Command {

    // Parameters.
    //@Parameter
    //OpService ops;

    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter
    ImagePlus img;

    @Parameter (description = "first slice for averaging", min = "1")
    Integer startSlice = 1;

    @Parameter (description = "last slice for averaging", min = "1")
    Integer endSlice = 30;

    @Parameter(description = "Directory to save results in", label = "Save Directory", style = "directory")
    File saveDirectory;

    // Member variables.
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private final String workingImagePathAndFileName = IJ.getDirectory("temp") + "-" + UUID.randomUUID() + "-scratch.tif";

    /**
     * Average an ImagePlus image stack.
     */
    public ImagePlus averageImagePlus(ImagePlus image, int iStart, int iEnd) {

        // FIXME: ImagePlus image could have frames instead of slices?

        // Assuming that ZProjector uses 1 based indexing.
        if (iStart < 1) { iStart = 1; }
        if (iStart > (image.getNSlices() - 1)) { iStart = image.getNSlices() - 1; }
        if (iEnd < (iStart+1)) { iEnd = iStart + 1; }
        if (iEnd > image.getNSlices()) { iEnd = image.getNSlices(); }

        log.info("averaging slices " + iStart + " to " + iEnd);
        ImagePlus averageImage = ZProjector.run(image, "ave", iStart, iEnd);
        return averageImage;
    }

    /**
     * Load an existing mapping file to initialize transformString and image size.
     */
    private String transformString = null;      // Transform string to pass to TurboReg.
    private int mapImageWidth = 0;              // Expected image width for the transform.
    private int mapImageHeight = 0;             // Expected image height for the transform.
    public void loadMappingJSON(File mappingFileName) {
        try {
            ObjectMapper mapper = new ObjectMapper();
            Map<String, Object> mapping = mapper.readValue(mappingFileName, HashMap.class);
            ArrayList<ArrayList <Double>> sourcePoints = (ArrayList) mapping.get("source points");
            ArrayList<ArrayList <Double>> targetPoints = (ArrayList) mapping.get("target points");

            this.mapImageWidth = (int) mapping.get("image width");
            this.mapImageHeight = (int) mapping.get("image height");
            this.transformString = sourcePoints.get(0).get(0) + " " + sourcePoints.get(0).get(1)
                                + " " + targetPoints.get(0).get(0) + " " + targetPoints.get(0).get(1)
                                + " " + sourcePoints.get(1).get(0) + " " + sourcePoints.get(1).get(1)
                                + " " + targetPoints.get(1).get(0) + " " + targetPoints.get(1).get(1)
                                + " " + sourcePoints.get(2).get(0) + " " + sourcePoints.get(2).get(1)
                                + " " + targetPoints.get(2).get(0) + " " + targetPoints.get(2).get(1);

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Split a (single frame) ImagePlus image vertically and return as a two element list (target, source).
     * Or maybe this would just process the first frame anyway?
     */
    public java.util.List<ImagePlus> splitImagePlus(ImagePlus image, boolean transform){

        // Record or check that the image size is correct for the current mapping.
        if (this.mapImageWidth == 0) {
            this.mapImageWidth = image.getWidth();
            this.mapImageHeight = image.getHeight();
        }
        else{
            if (this.mapImageWidth != image.getWidth()){
                throw new smFRETAnalysisException("Image width (" + image.getWidth() + ") doesn't match expected width for current mapping (" + this.mapImageWidth + ")");
            }
            if (this.mapImageHeight != image.getHeight()){
                throw new smFRETAnalysisException("Image height (" + image.getHeight() + ") doesn't match expected height for current mapping (" + this.mapImageHeight + ")");
            }
        }

        int hw = image.getWidth()/2;

        // Target image, left half, not transformed.
        image.setRoi(0,0,hw,image.getHeight());
        ImagePlus imageTarget = image.crop().duplicate();

        // Source image, right half.
        image.setRoi(hw,0,hw,image.getHeight());
        ImagePlus imageSource = image.crop().duplicate();

        // Transform source if requested.
        if (transform) {
            imageSource = transformImagePlus(imageSource);
        }

        // Image stack in order target, source.
        java.util.List<ImagePlus> images = new ArrayList<>();
        images.add(imageTarget);
        images.add(imageSource);

        return images;
    }

    /**
     * Copied from MultiStackReg_.java.
     */
    private String saveTempImageFile(ImagePlus sourceImg){
        FileSaver sourceFile = new FileSaver(sourceImg);
        String sourcePathAndFileName = IJ.getDirectory("temp") + "-" + UUID.randomUUID() + "-" + sourceImg.getTitle();
        sourceFile.saveAsTiff(sourcePathAndFileName);
        return sourcePathAndFileName;
    }

    /**
     * Transform an ImagePlus image with the current transform.
     */
    private ImagePlus transformImagePlus(ImagePlus image){
        if (this.transformString == null) {
            throw new smFRETAnalysisException("Error: Cannot transform image, transform string not set");
        }

        ImagePlus transformedImage = null;
        try {
            FileSaver sourceFile = new FileSaver(image);
            sourceFile.saveAsTiff(this.workingImagePathAndFileName);

            String options = " -transform"
                    + " -file " + this.workingImagePathAndFileName
                    + " " + image.getWidth() + " " + image.getHeight()
                    + " -affine " + this.transformString
                    + " -hideOutput";

            Object turboRegObject = IJ.runPlugIn("TurboReg_", options);
            Method method = turboRegObject.getClass().getMethod("getTransformedImage", null);
            transformedImage = (ImagePlus) method.invoke(turboRegObject, null);
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
        return transformedImage;
    }

    /**
     * Run the mapping identification process
     */
    @Override
    public void run() {
        try {
	        log.info("starting channel mapping " + img.getHeight() + " " + img.getWidth());
	    
            // Calculate average image.
            //
            // FIXME: Block or possibly handle RGB images.
            //
            log.info("calculating average image and splitting - " + img.getNSlices() + " slices");
            ImagePlus averageImage;
            if (img.getNSlices() == 1){
                averageImage = img.duplicate();
            }
            else {
                averageImage = averageImagePlus(img, startSlice, endSlice);
            }

            // Split average vertically.
            //
            // I spent a lot of time trying to figure out how to do this in memory in a way that was compatible
            // with TurboReg_ and finally gave up. It refused to use any ImagePlus objects whose titles
            // where changed, why IDK, so save images to temporary files following MultiStackReg_.
            //
            java.util.List<ImagePlus> images = splitImagePlus(averageImage, false);

            // Target image.
            ImagePlus averageImageTarget = images.get(0);
            averageImageTarget.setTitle("average_image_target");
            String averageImageTargetFilename = saveTempImageFile(averageImageTarget);
            log.info("target image " + averageImageTargetFilename);

            // Source image.
            ImagePlus averageImageSource = images.get(1);
            averageImageSource.setTitle("average_image_source");
            String averageImageSourceFilename = saveTempImageFile(averageImageSource);
            log.info("source image " + averageImageSourceFilename);

            /*
            if (!this.isHeadless) {
                ui.show(averageImageTarget);
                ui.show(averageImageSource);
            }
            */

            // Find correspondence using TurboReg.
            //
            // 'Landmark' coordinates are arranged in a triangle on the image, is this optimal?
            //
            String crop = " 0 0 " + (averageImageSource.getWidth() - 1) + " " + (averageImageSource.getHeight() - 1);
            int iw = averageImageSource.getWidth()/4;
            int ih = averageImageSource.getHeight()/4;
            String coords0 = " " + 2*iw + " " + ih + " " + 2*iw + " " + ih;
            String coords1 = " " + iw + " " + 3*ih + " " + iw + " " + 3*ih;
            String coords2 = " " + 3*iw + " " + 3*ih + " " + 3*iw + " " + 3*ih;

            log.info("identifying correspondence");
            String options = "-align"
                    + " -file " + averageImageSourceFilename
                    + crop
                    + " -file " + averageImageTargetFilename
                    + crop
                    + " -affine" + coords0 + coords1 + coords2
                    + " -hideOutput";
            Object turboRegObject = IJ.runPlugIn("TurboReg_", options);

            log.info("calculating affine transform matrix");

            // Get updated landmark points used for the best fit transform.
            Method method = turboRegObject.getClass().getMethod("getSourcePoints", null);
            double[][] sourcePoints = (double[][]) method.invoke(turboRegObject, null);
            method = turboRegObject.getClass().getMethod("getTargetPoints", null);
            double[][] targetPoints = (double[][]) method.invoke(turboRegObject, null);
            log.info(sourcePoints.length + " " + targetPoints[0].length);

            // Calculate affine transform matrix.
            /*
            double[][] sourcePointsT = MatrixUtils.createRealMatrix(sourcePoints).transpose().getData();
            double[][] targetPointsT = MatrixUtils.createRealMatrix(targetPoints).transpose().getData();
            AffineModel2D model = new AffineModel2D();

            double[] w = new double[sourcePointsT[0].length];
            Arrays.fill(w, 1.0);
            model.fit(sourcePointsT, targetPointsT, w);
            log.info(model);
             */
            log.info("saving results to " + saveDirectory);

            // Save mapping and expected image size as JSON.
            Map<String, Object> mapping = new HashMap<>();
            mapping.put("source points", sourcePoints);
            mapping.put("target points", targetPoints);
            mapping.put("image width", averageImage.getWidth());
            mapping.put("image height", averageImage.getHeight());

            ObjectMapper mapper = new ObjectMapper();
            File saveFile = new File(saveDirectory, "mapping.json");
            mapper.writeValue(saveFile, mapping);

            // We could have just loaded the transformed image from the current turboRegObject,
            // but we go the more complicated route to also test the other functionality of
            // this class.
            //
            // reload mapping to initialize transform string.
            //
            loadMappingJSON(saveFile);

            // Warp target image to source, convert Gray8 and overlay for user QC.
            log.info("calculating QC image");
            ImagePlus transformedSource = transformImagePlus(averageImageSource);

//            method = turboRegObject.getClass().getMethod("getTransformedImage", null);
//            ImagePlus transformedSource = (ImagePlus)method.invoke(turboRegObject, null);

            ImageConverter icSource = new ImageConverter(transformedSource);
            icSource.convertToGray8();

            ImageConverter icTarget = new ImageConverter(averageImageTarget);
            icTarget.convertToGray8();

            ImagePlus[] channels = new ImagePlus[3];
            channels[0] = averageImageTarget; // Red
            channels[1] = null; // Green
            channels[2] = transformedSource; // Blue

            ImagePlus rgbImageQCImage = RGBStackMerge.mergeChannels(channels, false);
            if (!this.isHeadless) {
                ui.show(rgbImageQCImage);
            }

            // Save QC image.
            FileSaver sourceFile = new FileSaver(rgbImageQCImage);
            String pathAndFileName = saveDirectory.getAbsolutePath() + File.separator + "mapping_qc_image.tif";
            log.info(pathAndFileName);
            sourceFile.saveAsTiff(pathAndFileName);

            log.info("finishing channel mapping");

            // Method debug code.
            loadMappingJSON(saveFile);
            log.info(this.transformString);

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
