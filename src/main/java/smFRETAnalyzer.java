/*
 * This class handles extracts the time traces using the mapping and spot locations.
 */

import com.fasterxml.jackson.databind.ObjectMapper;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.Filters3D;
import ij.process.ImageProcessor;
import net.imagej.ops.OpService;
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
    private final boolean diagnostic_mode = true;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private String saveRootName;
    private final smFRETSpotFinder smfsf = new smFRETSpotFinder();

    /**
     * Stack background estimation.
     * This is just a temporal boxcar filter.
     */
    public java.util.List<ImagePlus> backGroundEstimation(ImagePlus image) {

        // Mean filter on z axis.
        // FIXME: Use 0.5 for x/y radius?
        // FIXME: Not sure how edge pixels are handled. Could be an issue if there is a
        //        a large delta w/ time at the beginning/end of the movie.
        ImageStack imageZFlt = Filters3D.filter(image.getStack(), Filters3D.MEAN, 1, 1, backgroundAverageNFrames);

        // Split into two separate stacks, one for each channel.
        ImageStack targetBg = new ImageStack();
        ImageStack sourceBg = new ImageStack();

        IJ.showStatus("Background estimation..");
        ImagePlus targetImgEst = null;
        ImagePlus sourceImgEst = null;
        for (int i = 1; i <= imageZFlt.size(); i++){
            IJ.showProgress((double)i/((double)imageZFlt.size()));

            // Get slice from stack.
            ImageProcessor ipZ = imageZFlt.getProcessor(i);
            ImagePlus tmp = new ImagePlus("tmp", ipZ);

            // Split and transform source to target.
            java.util.List<ImagePlus> splitImages = smfsf.splitImagePlus(tmp);
            ImagePlus targetImg = splitImages.get(0);
            ImagePlus sourceImg = splitImages.get(1);

            // Estimate background in source and target.
            targetImgEst = smfsf.backgroundEstimate(targetImg, targetImgEst);
            sourceImgEst = smfsf.backgroundEstimate(sourceImg, sourceImgEst);

            // Add to stack.
            targetBg.addSlice(targetImgEst.getProcessor());
            sourceBg.addSlice(sourceImgEst.getProcessor());

            if(diagnostic_mode){
                log.info("");
            }
        }

        ImagePlus targetBgImp = new ImagePlus("target background estimate", targetBg);
        ImagePlus sourceBgImp = new ImagePlus("source background estimate", sourceBg);

        if (diagnostic_mode) {
            FileSaver fgSmoothImageSaver = new FileSaver(targetBgImp);
            fgSmoothImageSaver.saveAsTiff(saveRootName + "_fret_target_bg.tif");

            FileSaver bgSmoothImageSaver = new FileSaver(sourceBgImp);
            bgSmoothImageSaver.saveAsTiff(saveRootName + "_fret_source_bg.tif");
        }

        java.util.List<ImagePlus> images = new ArrayList<>();
        images.add(targetBgImp);
        images.add(sourceBgImp);

        IJ.showStatus("Processing complete.");
        return images;
    }

    /**
     * Trace extraction.
     *
     * Result data format is [time , spot], w/ per spot data in order target, source.
     */
    public double[][] measureTimeTraces(ImagePlus image, java.util.List<ImagePlus> bgEstimates, Polygon spots, double spotSigma, double cameraGain){
        ImageStack imageS = image.getStack();
        // Duplicating so we don't modify the input image.
        ImageStack targetBg = bgEstimates.get(0).getStack().duplicate();
        ImageStack sourceBg = bgEstimates.get(1).getStack().duplicate();
        double[][] timeTraces = new double[imageS.getSize()][2*spots.npoints];

        log.info("measuring time traces");
        double norm = 2.0*Math.PI*spotSigma*spotSigma;
        for (int i = 1; i <= imageS.size(); i++) {
            IJ.showProgress((double) i / ((double) imageS.size()));

            // Get slice from stack.
            ImageProcessor ipZ = imageS.getProcessor(i);
            ImagePlus tmp = new ImagePlus("tmp", ipZ);

            // Split, transform source to target and Gaussian smoothing.
            java.util.List<ImagePlus> splitImages = smfsf.splitImagePlus(tmp);
            ImagePlus targetImgI = splitImages.get(0);
            ImageProcessor targetImgIImp = targetImgI.getProcessor();
            targetImgIImp.blurGaussian(spotSigma);

            ImagePlus sourceImgI = splitImages.get(1);
            ImageProcessor sourceImgIImp = sourceImgI.getProcessor();
            sourceImgIImp.blurGaussian(spotSigma);

            // Gaussian smoothing of background.
            ImagePlus targetImgIBg = new ImagePlus("tmp_fg", targetBg.getProcessor(i));
            ImageProcessor targetImgIBgImp = targetImgIBg.getProcessor();
            targetImgIBgImp.blurGaussian(spotSigma);

            ImagePlus sourceImgIBg = new ImagePlus("tmp_bg", sourceBg.getProcessor(i));
            ImageProcessor sourceImgIBgImp = sourceImgIBg.getProcessor();
            sourceImgIBgImp.blurGaussian(spotSigma);

            // Record spot intensities in both channels.
            for (int j = 0; j < spots.npoints; j++) {
                int x = spots.xpoints[j];
                int y = spots.ypoints[j];

                timeTraces[i-1][2 * j] = norm * cameraGain * ((double) targetImgI.getPixel(x, y)[0] - (double) targetImgIBg.getPixel(x, y)[0]);
                timeTraces[i-1][2 * j + 1] = norm * cameraGain * ((double) sourceImgI.getPixel(x, y)[0] - (double) sourceImgIBg.getPixel(x, y)[0]);
            }
        }

        return timeTraces;
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
            saveRootName = (String) mapping.get("root name");
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
            smfsf.spotMargin = (Integer) mapping.get("spot margin");
            smfsf.edgeMargin = (Integer) mapping.get("edge margin");
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
            double[][] timeTraces = measureTimeTraces(image, bgEstimates, spots, spotSigma, cameraGain);
            log.info(timeTraces.length + " " + timeTraces[0].length);

            // Save time traces.

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
