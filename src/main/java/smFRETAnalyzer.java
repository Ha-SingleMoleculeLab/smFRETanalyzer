
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.ZProjector;

import mpicbg.models.AffineModel2D;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.IllDefinedDataPointsException;

import net.imagej.ops.OpService;

import org.apache.commons.math3.linear.MatrixUtils;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.UUID;

@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Channel Mapping")
public class smFRETAnalyzer implements Command {

    @Parameter
    OpService ops;

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

    /**
     * Copied from MultiStackReg_.java.
     */
    String saveTempImageFile(ImagePlus sourceImg){
        FileSaver sourceFile = new FileSaver(sourceImg);
        String sourcePathAndFileName = IJ.getDirectory("temp") + UUID.randomUUID() + sourceImg.getTitle();
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
            ImagePlus averageImage = ZProjector.run(img, "ave", startFrame, endFrame);

            // Split average vertically.
            log.info("calculating average image");
            int hw = averageImage.getWidth()/2;

            averageImage.setRoi(0,0,hw,averageImage.getHeight());
            ImagePlus averageImageSource = averageImage.crop().duplicate();
            averageImageSource.setTitle("average_image_source");
            String averageImageSourceFilename = saveTempImageFile(averageImageSource);
            log.info("source image " + averageImageSourceFilename);

            averageImage.setRoi(hw,0,hw,averageImage.getHeight());
            ImagePlus averageImageTarget = averageImage.crop().duplicate();
            averageImageTarget.setTitle("average_image_target");
            String averageImageTargetFilename = saveTempImageFile(averageImageTarget );
            log.info("target image " + averageImageTargetFilename);

            if (true) {
                ui.show(averageImageSource);
                ui.show(averageImageTarget);
            }

            // Find correspondence using TurboReg.
            // I much preferred version 2.0.0 where you didn't have to specify landmarks.
            //
            int iw = averageImageSource.getWidth();
            int ih = averageImageSource.getWidth();

            log.info("identifying correspondence");
            String options = "-align"
                    + " -file " + averageImageSourceFilename
                    + " 0 0 " + (averageImageSource.getWidth() - 1) + " " + (averageImageSource.getHeight() - 1)
                    + " -file " + averageImageTargetFilename
                    + " 0 0 " + (averageImageTarget.getWidth() - 1) + " " + (averageImageTarget.getHeight() - 1)
                    + " -affine 10 10 10 10 " + (ih - 10) + " 10 " + (ih - 10) + " 10 10 " + (iw - 10) + " 10 " + (iw - 10)
                    + " -hideOutput";
            Object turboRegObject = IJ.runPlugIn("TurboReg_", options);

            log.info("calculating affine transform matrix");

            // Calculate affine transform matrix.
            Method method = turboRegObject.getClass().getMethod("getSourcePoints", null);
            double[][] sourcePoints = (double[][]) method.invoke(turboRegObject, null);
            method = turboRegObject.getClass().getMethod("getTargetPoints", null);
            double[][] targetPoints = (double[][]) method.invoke(turboRegObject, null);
            log.info(sourcePoints.length + " " + targetPoints[0].length);

            double[][] sourcePointsT = MatrixUtils.createRealMatrix(sourcePoints).transpose().getData();
            double[][] targetPointsT = MatrixUtils.createRealMatrix(targetPoints).transpose().getData();

            AffineModel2D model = new AffineModel2D();
            double[] w = new double[sourcePointsT[0].length];
            Arrays.fill(w, 1.0);
            model.fit(sourcePointsT, targetPointsT, w);
            log.info(model);
            //for (e in b) {
            //    println(model.applyInverse((double[]) e));
            //}

            log.info("ending channel mapping");

            // UI for user verification of found transform.
        } catch (NoSuchMethodException |
                 IllegalAccessException |
                 InvocationTargetException |
                 NotEnoughDataPointsException |
                 IllDefinedDataPointsException e){
            log.info(e);
        }
    }
}
