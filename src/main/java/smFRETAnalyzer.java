/*
 * This class handles extracts the time traces using the mapping and spot locations.
 */

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.io.LittleEndianDataOutputStream;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.Filters3D;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import net.imagej.ops.OpService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;
import java.io.FileOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, headless = true,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Time Traces", weight = 3.0)})
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
    private final boolean saveAsTraces = true;
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
        int stepSize = Math.max(imageZFlt.size()/100, 1);
        for (int i = 1; i <= imageZFlt.size(); i++){
            if (i%stepSize == 0){
                log.info("processing " + i + " of " + imageZFlt.size() + " images");
            }
            IJ.showProgress((double)i/((double)imageZFlt.size()));

            // Get slice from stack.
            ImageProcessor ipZ = imageZFlt.getProcessor(i);
            ImagePlus tmp = new ImagePlus("tmp", ipZ);

            // Split and transform source to target.
            java.util.List<ImagePlus> splitImages = smfsf.splitImagePlus(tmp);
            ImagePlus targetImg = splitImages.get(0);
            ImagePlus sourceImg = splitImages.get(1);

            // Estimate background in source and target. Each frame is estimated from scratch;
            // the estimator settles in a fixed few rounds and has nothing to carry over from
            // the previous frame the way the old inpainting fill did.
            ImagePlus targetImgEst = smfsf.backgroundEstimate(targetImg);
            ImagePlus sourceImgEst = smfsf.backgroundEstimate(sourceImg);

            // Add to stack.
            targetBg.addSlice(targetImgEst.getProcessor());
            sourceBg.addSlice(sourceImgEst.getProcessor());

            //log.info("");
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
     * Result data format is [spot, time], w/ per spot time data in order target, source.
     */
    public double[][] measureTimeTraces(ImagePlus image, java.util.List<ImagePlus> bgEstimates, double[][] spots, double spotSigma, double cameraGain){
        ImageStack imageS = image.getStack();

        // Not duplicated any more: smoothed() copies before it blurs, so nothing here writes to
        // the stacks it was handed.
        ImageStack targetBg = bgEstimates.get(0).getStack();
        ImageStack sourceBg = bgEstimates.get(1).getStack();
        double[][] timeTraces = new double[2*spots.length][imageS.getSize()];

        IJ.showStatus("Measuring time traces..");
        double norm = 2.0*Math.PI*spotSigma*spotSigma;
        int stepSize = Math.max(imageS.size()/100, 1);
        for (int i = 1; i <= imageS.size(); i++) {
            if (i%stepSize == 0){
                log.info("processing " + i + " of " + imageS.size() + " images");
            }
            IJ.showProgress((double) i / ((double) imageS.size()));

            // Get slice from stack.
            ImageProcessor ipZ = imageS.getProcessor(i);
            ImagePlus tmp = new ImagePlus("tmp", ipZ);

            // Split, transform source to target and Gaussian smoothing.
            java.util.List<ImagePlus> splitImages = smfsf.splitImagePlus(tmp);
            smFRETSpotFinder.Shared targetImgI = smoothed(splitImages.get(0), spotSigma);
            smFRETSpotFinder.Shared sourceImgI = smoothed(splitImages.get(1), spotSigma);

            // Gaussian smoothing of background.
            smFRETSpotFinder.Shared targetImgIBg = smoothed(targetBg.getProcessor(i), spotSigma);
            smFRETSpotFinder.Shared sourceImgIBg = smoothed(sourceBg.getProcessor(i), spotSigma);

            // Record spot intensities in both channels. Reading the backing array directly is the
            // same value ImageProcessor.getValue returned, without going through ImagePlus - but
            // getValue handed back a double, so the subtraction has to be widened by hand. Left
            // in float it lands 1e-4 away, which is nothing next to the signal and still enough
            // to make this migration visible in the traces.
            for (int j = 0; j < spots.length; j++) {
                int x = (int)spots[j][0];
                int y = (int)spots[j][1];
                int index = y * targetImgI.width + x;

                timeTraces[2 * j][i - 1] = norm * cameraGain
                        * ((double) targetImgI.pixels[index] - (double) targetImgIBg.pixels[index]);
                timeTraces[2 * j + 1][i - 1] = norm * cameraGain
                        * ((double) sourceImgI.pixels[index] - (double) sourceImgIBg.pixels[index]);
            }
        }

        return timeTraces;
    }

    /**
     * A float copy of an image, smoothed at the spot scale.
     *
     * The blur is ImageJ1's, deliberately: Gauss3 differs from it by enough to move these traces,
     * measured at 0.04 ADU at this sigma. Shared lets the blur run over the same array the rest
     * of the pipeline reads as imglib2, so there is no conversion around it. The copy is what
     * makes it safe to hand this a processor belonging to a stack.
     */
    private static smFRETSpotFinder.Shared smoothed(ImagePlus image, double sigma) {
        return smoothed(image.getProcessor(), sigma);
    }

    private static smFRETSpotFinder.Shared smoothed(ImageProcessor image, double sigma) {
        smFRETSpotFinder.Shared copy = new smFRETSpotFinder.Shared(
                (FloatProcessor) image.convertToFloatProcessor().duplicate());
        copy.processor.blurGaussian(sigma);
        return copy;
    }

    /**
     * Save in HDF5 format.
     */
    private void saveToHDF5File(String hdf5FileName, double [][] timeTraces, double[][] spots) {
        try {
            IHDF5Writer writer = HDF5Factory.configure(hdf5FileName).writer();

            // Some analysis metadata.
            writer.writeString("spot-json-file", spotJSONFile.toString());
            writer.writeInt("background-average-n-frames", backgroundAverageNFrames.shortValue());

            String spotJSONFileContents = new String(Files.readAllBytes(Paths.get(spotJSONFile.toString())));
            writer.writeString("spot-json-file-contents", spotJSONFileContents);

            // Save spot locations (in target channel).
            writer.writeString("spots-fields", smfsf.columnHeaders.toString());
            writer.writeDoubleMatrix("spots", spots);

            // Split time trace data into target, source so that indexing matches spots.
            float[][] targetTraces = new float[timeTraces.length / 2][timeTraces[0].length];
            float[][] sourceTraces = new float[timeTraces.length / 2][timeTraces[0].length];
            for (int i = 0; i < targetTraces.length; i++){
                for (int j = 0; j < targetTraces[0].length; j++) {
                    targetTraces[i][j] = (float)timeTraces[2*i][j];
                    sourceTraces[i][j] = (float)timeTraces[2*i+1][j];
                }
            }
            writer.writeFloatMatrix("target-traces", targetTraces);
            writer.writeFloatMatrix("source-traces", sourceTraces);

            writer.close();
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Save in the .traces format, which is:
     *
     * int32 (Int) - length of traces.
     * int16 (Short) - number of traces.
     * int16 (Short) - 2 * (length of traces) * (number of traces) in order [time][trace number].
     */
    private void saveToTracesFile(String tracesFileName, double [][] timeTraces){
        try (LittleEndianDataOutputStream tracesDos = new LittleEndianDataOutputStream(new FileOutputStream(tracesFileName))) {
            tracesDos.writeInt(timeTraces[0].length);   // Length of traces.
            tracesDos.writeShort(timeTraces.length / 2); // Number of traces.

            for (int j = 0; j < timeTraces[0].length; j++) {
                for (int i = 0; i < timeTraces.length; i++) {
                    tracesDos.writeShort((short) Math.round(timeTraces[i][j]));
                }
            }
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
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
            // Neither the spot margin nor the edge margin is set here any more.
            // Both only ever mattered to background estimation, and that now
            // takes its mask from the masks file and its smoothing from a
            // constant, so nothing downstream of loadMasks reads either one.
            // Absent from spot finder JSON written before background estimation changed, in
            // which case the estimator's own default stands.
            Object kappa = mapping.get("background kappa");
            if (kappa != null) {
                smfsf.backgroundKappa = ((Number) kappa).doubleValue();
            }
            smfsf.loadMappingJSON(mappingFileName);
            smfsf.loadMasks(masksFileName);
            double[][] spots = smfsf.loadSpotLocations(spotsFileName);
            log.info("loaded " + spots.length + " spots");

            // Load image to process.
            ImagePlus image = new ImagePlus(inputImageName);

            // Estimate background in the two image channels.
            log.info("estimating background in channels");
            java.util.List<ImagePlus> bgEstimates = backGroundEstimation(image);

            if (!isHeadless) {
                ui.show(bgEstimates.get(0));
                ui.show(bgEstimates.get(1));
            }

            // Measure spot time traces.
            log.info("measuring time traces");
            double[][] timeTraces = measureTimeTraces(image, bgEstimates, spots, spotSigma, cameraGain);
            log.info(timeTraces.length + " " + timeTraces[0].length);

            // Save time traces in '.traces' format:
            if (saveAsTraces){
                saveToTracesFile(saveRootName + ".traces", timeTraces);
            }

            // Save time traces, spot locations and metadata to .h5 file.
            saveToHDF5File(saveRootName + ".h5", timeTraces, spots);

            log.info("finishing time trace measurement");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
