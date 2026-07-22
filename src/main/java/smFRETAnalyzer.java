/*
 * This class handles extracts the time traces using the mapping and spot locations.
 */

import com.fasterxml.jackson.databind.ObjectMapper;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Filters3D;
import ij.process.ImageProcessor;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Time Traces")
public class smFRETAnalyzer implements Command {

    // Parameters.
    @Parameter
    OpService ops;

    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter(description = "spot finding JSON output file", label = "Spot Finder JSON file", style = "open")
    File spotJSONFile;

    @Parameter (description = "number of frames to use for background averaging estimation", min = "1")
    Integer backgroundAverageNFrames = 30;

    // Member variables.
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private smFRETSpotFinder smfsf = new smFRETSpotFinder();

    /**
     * Stack background estimation.
     * This is just a temporal boxcar filter.
     */
    public java.util.List<ImagePlus> backGroundEstimation(ImagePlus image) {

        // Mean filter on z axis.
        Filters3D flt3D = new Filters3D();
        ImageStack imageS = image.getStack();
        ImageStack imageZFlt = flt3D.filter(imageS, Filters3D.MEAN, 1, 1, backgroundAverageNFrames);
        ImagePlus imageZFltImp = new ImagePlus("time smoothed", imageZFlt);

        // Split into two separate stacks, one for each channel.
        ImageStack targetBg = new ImageStack();
        ImageStack sourceBg = new ImageStack();

        IJ.showStatus("Background estimation..");
        for (int i = 1; i <= imageZFlt.size(); i++){
            //IJ.showProgress(i, imageZFlt.size());
            IJ.showProgress((double)i/((double)imageZFlt.size()));

            // Get slice from stack.
            ImageProcessor ipZ = imageZFltImp.getStack().getProcessor(i);
            ImagePlus tmp = new ImagePlus("tmp", ipZ);

            // Split and transform source to target.
            java.util.List<ImagePlus> splitImages = smfsf.splitImagePlus(tmp);
            ImagePlus targetImg = splitImages.get(0);
            ImagePlus sourceImg = splitImages.get(1);

            // Estimate background in source and target.
            ImagePlus targetImgEst = smfsf.backgroundEstimate(targetImg);
            ImagePlus sourceImgEst = smfsf.backgroundEstimate(sourceImg);

            // Add to stack.
            targetBg.addSlice(targetImgEst.getProcessor());
            sourceBg.addSlice(sourceImgEst.getProcessor());
        }

        ImagePlus targetBgImp = new ImagePlus("target background estimate", targetBg);
        ImagePlus sourceBgImp = new ImagePlus("source background estimate", sourceBg);

        java.util.List<ImagePlus> images = new ArrayList<>();
        images.add(targetBgImp);
        images.add(sourceBgImp);

        IJ.showStatus("Processing complete.");
        return images;
    }

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
            log.info("starting time trace measurement");

            // Load spot JSON file w/ the analysis parameters.
            ObjectMapper mapper = new ObjectMapper();
            Map<String, Object> mapping = mapper.readValue(spotJSONFile, HashMap.class);
            String saveRootName = (String) mapping.get("root name");
            String inputImageName = (String) mapping.get("image name");
            String mappingFileName = (String) mapping.get("mapping file");
            String masksFileName = (String) mapping.get("masks file");
            String spotsFileName = (String) mapping.get("spots file");
            double spotSigma = ((Number) mapping.get("spot sigma")).doubleValue();
            double cameraBlackLevel = ((Number) mapping.get("camera black")).doubleValue();
            double cameraGain = ((Number) mapping.get("camera gain")).doubleValue();

            log.info("initializing spot finder and loading images");

            // Initialize spot finder object.
            smfsf.log = log;
            smfsf.loadMappingJSON(mappingFileName);
            smfsf.loadMasks(masksFileName);
            Polygon spots = smfsf.loadSpotLocations(spotsFileName);
            log.info("loaded " + spots.npoints + " spots");

            // Load image to process.
            ImagePlus image = new ImagePlus(inputImageName);

            // Estimate background in the two image channels.
            log.info("estimating background in channels");
            java.util.List<ImagePlus> bgEstimates = backGroundEstimation(image);

            ui.show(bgEstimates.get(0));
            ui.show(bgEstimates.get(1));

            // Measure spot time traces.
            log.info("measuring time traces");

            // Save time traces.

            /*
            if (!this.isHeadless) {
                ui.show(smfsf.overlapMask);
                ui.show(smfsf.backgroundMask);
            }
            */

            log.info("finishing time trace measurement");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}


            /*
            float[] data = new float[10];
            for(int i = 0; i < 10; i++){
                data[i] = (float)0.1;
            }
            long[] dimensions = new long[]{1, 1, 10};
            RandomAccessibleInterval<FloatType> kernel = ArrayImgs.floats(data, dimensions);

            Img<FloatType> imageW = ImageJFunctions.wrap(image);
            //RandomAccessible<FloatType> imageWExt = Views.extendMirrorSingle(imageW);
            Img<FloatType> imageZBoxcar = ops.create().img(imageW);
            //Img<FloatType>imageZBoxcar = (Img<FloatType>) ops.filter().convolve(imageW, kernel);
            ops.filter().convolve(imageZBoxcar, imageW, kernel);
            */
