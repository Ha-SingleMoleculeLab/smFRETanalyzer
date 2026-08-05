/*
 * This class handles everything related to mapping the two channels to the same space.
 * It depends on the TurboReg plugin (version 2.0.1) to *measure* a mapping. Applying one is
 * done here with mpicbg, which is a compile time dependency - see transformImagePlus.
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
import ij.process.FloatProcessor;
import ij.process.ImageConverter;

import ij.process.ImageProcessor;
import mpicbg.ij.InverseTransformMapping;
import mpicbg.models.AffineModel2D;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;


// Menu weights put the three batch plugins in the order they are meant to be run, rather than
// leaving them to sort alphabetically as they would with an unweighted menuPath.
@Plugin(type = Command.class, headless = true,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Channel Mapping", weight = 1.0)})
public class smFRETChannelMapper implements Command {

    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter(description = "image to analyze to identify the mapping", label = "Image", style = "open")
    File inputImageName;

    @Parameter (description = "first slice for averaging", label = "Start Slice", min = "1")
    Integer startSlice = 1;

    @Parameter (description = "last slice for averaging", label = "End Slice", min = "1")
    Integer endSlice = 30;

    // Member variables.
    private final boolean diagnostic_mode = true;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private int mapImageWidth = 0;              // Expected image width for the transform.
    private int mapImageHeight = 0;             // Expected image height for the transform.
    private AffineModel2D transformModel = null;    // Source to target affine, from the landmark pairs.

    /**
     * Average an ImagePlus image stack.
     */
    public ImagePlus averageImagePlus(ImagePlus image, int iStart, int iEnd) {

        // FIXME: ImagePlus image could have frames instead of slices?

        if (image.getNSlices() == 1){
            return image.duplicate();
        }

        // Assuming that ZProjector uses 1 based indexing.
        if (iStart < 1) { iStart = 1; }
        if (iStart > (image.getNSlices() - 1)) { iStart = image.getNSlices() - 1; }
        if (iEnd < (iStart+1)) { iEnd = iStart + 1; }
        if (iEnd > image.getNSlices()) { iEnd = image.getNSlices(); }

        log.info("averaging slices " + iStart + " to " + iEnd);
        return ZProjector.run(image, "ave", iStart, iEnd);
    }

    /**
     * Load an existing mapping file to initialize the transform model and image size.
     *
     * The three landmark pairs determine the affine exactly, so the fit is a solve rather than a
     * regression and the weights are all one. TurboReg derives the same affine from the same three
     * pairs - measured agreement is 2e-13 pixels - so this reads mapping files written either
     * before or after transformImagePlus stopped calling TurboReg.
     */
    public void loadMappingJSON(String mappingFileName) {
        try {
            ObjectMapper mapper = new ObjectMapper();
            File mappingFile = new File(mappingFileName);
            Map<String, Object> mapping = mapper.readValue(mappingFile, HashMap.class);
            ArrayList<ArrayList <Double>> sourcePoints = (ArrayList) mapping.get("source points");
            ArrayList<ArrayList <Double>> targetPoints = (ArrayList) mapping.get("target points");

            mapImageWidth = (int) mapping.get("image width");
            mapImageHeight = (int) mapping.get("image height");

            // mpicbg wants the coordinates grouped by axis, not by point.
            double[][] source = new double[2][3];
            double[][] target = new double[2][3];
            for (int i = 0; i < 3; i++) {
                source[0][i] = sourcePoints.get(i).get(0);
                source[1][i] = sourcePoints.get(i).get(1);
                target[0][i] = targetPoints.get(i).get(0);
                target[1][i] = targetPoints.get(i).get(1);
            }

            AffineModel2D model = new AffineModel2D();
            model.fit(source, target, new double[]{1.0, 1.0, 1.0});
            transformModel = model;

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Split a (single frame) ImagePlus image vertically and return as a two element list (target, source).
     * Or maybe this would just process the first frame anyway?
     *
     * FIXME: Can return as <ImagePlus,ImagePlus>?
     */
    public java.util.List<ImagePlus> splitImagePlus(ImagePlus image, boolean transform){

        // Record or check that the image size is correct for the current mapping.
        if (mapImageWidth == 0) {
            mapImageWidth = image.getWidth();
            mapImageHeight = image.getHeight();
        }
        else{
            if (mapImageWidth != image.getWidth()){
                throw new smFRETAnalysisException("Image width (" + image.getWidth() + ") doesn't match expected width for current mapping (" + mapImageWidth + ")");
            }
            if (mapImageHeight != image.getHeight()){
                throw new smFRETAnalysisException("Image height (" + image.getHeight() + ") doesn't match expected height for current mapping (" + mapImageHeight + ")");
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
     * https://github.com/miura/MultiStackRegistration
     */
    private String saveTempImageFile(ImagePlus sourceImg){
        FileSaver sourceFile = new FileSaver(sourceImg);
        String sourcePathAndFileName = IJ.getDirectory("temp") + "-" + UUID.randomUUID() + "-" + sourceImg.getTitle();
        sourceFile.saveAsTiff(sourcePathAndFileName);
        return sourcePathAndFileName;
    }

    /**
     * Transform an ImagePlus image with the current transform.
     *
     * This used to hand the image to TurboReg's -transform through a temporary file. mpicbg does
     * the same job in memory, which is worth doing because this runs once per frame per channel
     * from both smFRETAnalyzer stages, and because TurboReg is not a compile time dependency so
     * nothing here could be tested without a Fiji install.
     *
     * The two are not bit identical - TurboReg interpolates with a cubic spline and this is
     * bilinear - but measured against TurboReg on the example mapping and data, the geometry is
     * the same to 0.003 pixels of spot centroid, the stored 16 bit values differ for 14% of
     * pixels and almost always by exactly 1 ADU, and per-spot intensities after the smoothing
     * that precedes trace measurement differ by a median of 0.45%. The rounding in
     * convertToShort below has an rms of 0.26 ADU against 0.24 ADU between the interpolators, so
     * the choice of interpolator perturbs the stored image less than storing it does.
     *
     * Out of range pixels are left at zero, which is what createOverlapMask relies on to find
     * the region the two channels share.
     */
    private ImagePlus transformImagePlus(ImagePlus image){
        if (transformModel == null) {
            throw new smFRETAnalysisException("Error: Cannot transform image, transform model not set");
        }

        // Interpolate in float regardless of the input type, so the result is not quantized twice.
        FloatProcessor source = image.getProcessor().convertToFloatProcessor();
        source.setInterpolationMethod(ImageProcessor.BILINEAR);

        FloatProcessor target = new FloatProcessor(image.getWidth(), image.getHeight());
        new InverseTransformMapping<AffineModel2D>(transformModel).mapInterpolated(source, target);

        // convert to 16 bit.
        return new ImagePlus("transformed image", target.convertToShort(false));
    }

    /**
     * Run the mapping identification process
     */
    @Override
    public void run() {
        try {
	        log.info("starting channel mapping from image " + inputImageName);

            // Root name to use for saving output, this is just the file name
            // without the extension.
            String saveRootName = inputImageName.toString();
            int dotIndex = saveRootName.lastIndexOf('.');
            if (dotIndex > 0) {
                saveRootName = saveRootName.substring(0, dotIndex);
            }
            log.info("save root " + saveRootName);

            // Load the image to process.
            //Opener sourceOpener = new Opener();
            //ImagePlus inputImage = sourceOpener.openImage(inputImageName.toString());
            ImagePlus inputImage = new ImagePlus(inputImageName.toString());

            log.info("starting channel mapping " + inputImage.getHeight() + " " + inputImage.getWidth());

            // Calculate average image.
            //
            // FIXME: Block or possibly handle RGB images.
            //
            log.info("calculating average image and splitting - " + inputImage.getNSlices() + " slices");
            ImagePlus averageImage = averageImagePlus(inputImage, startSlice, endSlice);

            FileSaver averageImageFileSaver = new FileSaver(averageImage);
            if (diagnostic_mode){
                averageImageFileSaver.saveAsTiff(saveRootName + "_mapping_average_image.tif");
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

            FileSaver targetImageFileSaver = new FileSaver(averageImageTarget);
            if (diagnostic_mode){
                targetImageFileSaver.saveAsTiff(saveRootName + "_mapping_average_target.tif");
            }

            // Source image.
            ImagePlus averageImageSource = images.get(1);
            averageImageSource.setTitle("average_image_source");
            String averageImageSourceFilename = saveTempImageFile(averageImageSource);
            log.info("source image " + averageImageSourceFilename);

            FileSaver sourceImageFileSaver = new FileSaver(averageImageSource);
            if (diagnostic_mode){
                sourceImageFileSaver.saveAsTiff(saveRootName + "_mapping_average_source.tif");
            }

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
            log.info("saving results");

            // Save mapping and expected image size as JSON.
            Map<String, Object> mapping = new HashMap<>();
            mapping.put("source points", sourcePoints);
            mapping.put("target points", targetPoints);
            mapping.put("image width", averageImage.getWidth());
            mapping.put("image height", averageImage.getHeight());

            ObjectMapper mapper = new ObjectMapper();
            File saveFile = new File(saveRootName + "_mapping.json");
            mapper.writeValue(saveFile, mapping);

            // We could have just loaded the transformed image from the current turboRegObject,
            // but we go the more complicated route to also test the other functionality of
            // this class.
            //
            // reload mapping to initialize transform string.
            //
            loadMappingJSON(saveFile.toString());

            // Warp target image to source, convert Gray8 and overlay for user QC.
            log.info("calculating QC image");
            ImagePlus transformedSource = transformImagePlus(averageImageSource);

            ImageConverter icSource = new ImageConverter(transformedSource);
            icSource.convertToGray8();

            ImageConverter icTarget = new ImageConverter(averageImageTarget);
            icTarget.convertToGray8();

            ImagePlus[] channels = new ImagePlus[3];
            channels[0] = averageImageTarget; // Red
            channels[1] = null; // Green
            channels[2] = transformedSource; // Blue

            ImagePlus rgbImageQCImage = RGBStackMerge.mergeChannels(channels, false);
            if (!isHeadless) {
                ui.show(averageImage);
                ui.show(rgbImageQCImage);
            }
            rgbImageQCImage.setTitle("mapping QC image");

            // Save QC image.
            FileSaver qcImageFileSaver = new FileSaver(rgbImageQCImage);
            String pathAndFileName = saveRootName + "_mapping_qc_image.tif";
            log.info(pathAndFileName);
            qcImageFileSaver.saveAsTiff(pathAndFileName);

            log.info("finishing channel mapping");

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
