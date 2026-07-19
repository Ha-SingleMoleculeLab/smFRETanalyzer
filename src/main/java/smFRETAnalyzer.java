
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ImageConverter;

import mpicbg.models.AffineModel2D;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.IllDefinedDataPointsException;

import net.imagej.ops.OpService;

import org.apache.commons.math3.linear.MatrixUtils;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.menu.MenuConstants;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;


/*
@Plugin(type = Command.class, label = "smFRET Channel Mapping", menu = {
        @Menu(label = MenuConstants.PLUGINS_LABEL, weight = 1, mnemonic = MenuConstants.PLUGINS_MNEMONIC),
        @Menu(label = "smFRET Analysis", weight = 1, mnemonic = 'r'),
        @Menu(label = "smFRET Channel Mapping", weight = 1) })
 */
@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Channel Mapping")
public class smFRETAnalyzer implements Command {

    //@Parameter
    //OpService ops;

    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter
    ImagePlus img;

    @Parameter
    Integer startFrame;

    @Parameter
    Integer endFrame;

    @Parameter(label = "Save Directory", style = "directory")
    File saveDirectory;

    /**
     * Load an existing mapping file to initialize transformString.
     */
    private String transformString = null;
    public void loadMappingJSON(File mappingFileName) {
        try {
            ObjectMapper mapper = new ObjectMapper();

            Map<String, Object> mapping = mapper.readValue(mappingFileName, HashMap.class);
            ArrayList<ArrayList <Double>> sourcePoints = (ArrayList) mapping.get("source points");
            ArrayList<ArrayList <Double>> targetPoints = (ArrayList) mapping.get("target points");

            this.transformString = sourcePoints.get(0).get(0) + " " + sourcePoints.get(0).get(1)
                                + " " + targetPoints.get(0).get(0) + " " + targetPoints.get(0).get(1)
                                + " " + sourcePoints.get(1).get(0) + " " + sourcePoints.get(1).get(1)
                                + " " + targetPoints.get(1).get(0) + " " + targetPoints.get(1).get(1)
                                + " " + sourcePoints.get(2).get(0) + " " + sourcePoints.get(2).get(1)
                                + " " + targetPoints.get(2).get(0) + " " + targetPoints.get(2).get(1);
        } catch (Exception e) {
            log.info(e);
        }
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
     * Run the mapping identification process
     */
    @Override
    public void run() {
        //RandomAccessibleInterval<FloatType> imgFloat = ops.convert().float32(img);

        log.info("starting channel mapping " + img.getHeight() + " " + img.getWidth());

        try {
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
                // Assuming that ZProjector uses 1 based indexing.
                if (startFrame < 1){ startFrame = 1;}
                if (startFrame > (img.getNFrames()-1)){ startFrame = img.getNFrames() - 1;}
                if (endFrame < startFrame){ endFrame = startFrame + 1;}
                if (endFrame > img.getNFrames()){ endFrame = img.getNFrames();}

                log.info("averaging slices " + startFrame + " to " + endFrame);
                averageImage = ZProjector.run(img, "ave", startFrame, endFrame);
            }

            // Split average vertically.
            //
            // I spent a lot of time trying to figure out how to do this in memory in a way that was compatible
            // with TurboReg_ and finally gave up. It refused to use any ImagePlus objects whose titles
            // where changed, why IDK, so save images to temporary files following MultiStackReg_.
            //
            int hw = averageImage.getWidth()/2;

            // Target image, left half, not transformed.
            averageImage.setRoi(0,0,hw,averageImage.getHeight());
            ImagePlus averageImageTarget = averageImage.crop().duplicate();
            averageImageTarget.setTitle("average_image_target");
            String averageImageTargetFilename = saveTempImageFile(averageImageTarget);
            log.info("target image " + averageImageTargetFilename);

            // Source image, right half, transformed.
            averageImage.setRoi(hw,0,hw,averageImage.getHeight());
            ImagePlus averageImageSource = averageImage.crop().duplicate();
            averageImageSource.setTitle("average_image_source");
            String averageImageSourceFilename = saveTempImageFile(averageImageSource);
            log.info("source image " + averageImageSourceFilename);

            ui.show(averageImageTarget);
            ui.show(averageImageSource);

            // Find correspondence using TurboReg.
            //
            // 'Landmark' coordinates are arranged in a triangle on the image, is this optimal?
            //
            String crop = " 0 0 " + (averageImageSource.getWidth() - 1) + " " + (averageImageSource.getHeight() - 1);
            int iw = averageImageSource.getWidth()/4;
            int ih = averageImageSource.getWidth()/4;
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

            // Calculate affine transform matrix.
            Method method = turboRegObject.getClass().getMethod("getSourcePoints", null);
            double[][] sourcePoints = (double[][]) method.invoke(turboRegObject, null);
            method = turboRegObject.getClass().getMethod("getTargetPoints", null);
            double[][] targetPoints = (double[][]) method.invoke(turboRegObject, null);
            log.info(sourcePoints.length + " " + targetPoints[0].length);

            /*
            double[][] sourcePointsT = MatrixUtils.createRealMatrix(sourcePoints).transpose().getData();
            double[][] targetPointsT = MatrixUtils.createRealMatrix(targetPoints).transpose().getData();
            AffineModel2D model = new AffineModel2D();

            double[] w = new double[sourcePointsT[0].length];
            Arrays.fill(w, 1.0);
            model.fit(sourcePointsT, targetPointsT, w);
            log.info(model);
             */

            // Warp target image to source, convert Gray8 and overlay for user QC.
            log.info("calculating QC image");
            method = turboRegObject.getClass().getMethod("getTransformedImage", null);
            ImagePlus transformedSource = (ImagePlus)method.invoke(turboRegObject, null);

            ImageConverter icSource = new ImageConverter(transformedSource);
            icSource.convertToGray8();

            ImageConverter icTarget = new ImageConverter(averageImageTarget);
            icTarget.convertToGray8();

            ImagePlus[] channels = new ImagePlus[3];
            channels[0] = averageImageTarget; // Red
            channels[1] = null; // Green
            channels[2] = transformedSource; // Blue

            ImagePlus rgbImageQCImage = RGBStackMerge.mergeChannels(channels, false);
            ui.show(rgbImageQCImage);

            log.info("saving results to " + saveDirectory);

            // Save QC image.
            FileSaver sourceFile = new FileSaver(rgbImageQCImage);
            String pathAndFileName = saveDirectory.getAbsolutePath() + File.separator + "mapping_qc_image.tif";
            log.info(pathAndFileName);
            sourceFile.saveAsTiff(pathAndFileName);

            // Save mapping as JSON.
            Map<String, Object> mapping = new HashMap<>();
            mapping.put("source points", sourcePoints);
            mapping.put("target points", targetPoints);

            ObjectMapper mapper = new ObjectMapper();
            File saveFile = new File(saveDirectory, "mapping.json");
            mapper.writeValue(saveFile, mapping);

            log.info("ending channel mapping");

            // Method debug code.
            loadMappingJSON(saveFile);
            log.info(this.transformString);

            // UI for user verification of found transform.
        } catch (Exception e) {
            log.info(e);
        }
/*
        } catch (NoSuchMethodException |
                 IllegalAccessException |
                 InvocationTargetException |
                 // NotEnoughDataPointsException
                 // IllDefinedDataPointsException
                 IOException e){
            log.info(e);
        }
*/
    }
}
