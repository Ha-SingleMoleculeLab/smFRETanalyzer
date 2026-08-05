/*
 * This class handles everything related to mapping the two channels to the same space.
 *
 * The pixel work - projecting, splitting, warping - is imglib2. Two things are not, and cannot
 * be: TurboReg, which *measures* a mapping and only speaks ImagePlus through files, and the
 * ImagePlus in and out of the public methods, which is what smFRETSpotFinder and smFRETAnalyzer
 * still pass. File writing stays ImageJ1 as well, deliberately - swapping FileSaver for SCIFIO
 * would change the TIFFs that existing analyses read back.
 */

import com.fasterxml.jackson.databind.ObjectMapper;

import java.awt.*;
import java.io.File;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.util.*;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.RGBStackMerge;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.converter.Converters;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
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
    private AffineTransform2D transformModel = null;    // Source to target affine, from the landmark pairs.

    /**
     * Average an ImagePlus image stack.
     *
     * Replaces ZProjector, and reproduces it rather than improving on it: the mean is accumulated
     * in float and the result is rounded back into the input's type, so an 8 bit movie still
     * averages to an 8 bit image. Keeping the mean as float would be defensible on its own terms
     * but would change every spot the next stage finds.
     */
    public ImagePlus averageImagePlus(ImagePlus image, int iStart, int iEnd) {

        // FIXME: ImagePlus image could have frames instead of slices?

        if (image.getNSlices() == 1){
            return image.duplicate();
        }

        // Assuming 1 based indexing, as ZProjector used.
        if (iStart < 1) { iStart = 1; }
        if (iStart > (image.getNSlices() - 1)) { iStart = image.getNSlices() - 1; }
        if (iEnd < (iStart+1)) { iEnd = iStart + 1; }
        if (iEnd > image.getNSlices()) { iEnd = image.getNSlices(); }

        log.info("averaging slices " + iStart + " to " + iEnd);

        Img<FloatType> stack = ImageJFunctions.convertFloat(image);
        Img<FloatType> mean = ArrayImgs.floats(image.getWidth(), image.getHeight());
        for (int slice = iStart - 1; slice < iEnd; slice++) {
            RandomAccessibleInterval<FloatType> plane = Views.hyperSlice(stack, 2, slice);
            LoopBuilder.setImages(mean, plane).forEachPixel((m, p) -> m.set(m.get() + p.get()));
        }

        final float n = iEnd - iStart + 1;
        LoopBuilder.setImages(mean).forEachPixel(m -> m.set(m.get() / n));

        ImagePlus averaged = new ImagePlus("average", toProcessor(mean, image.getBitDepth()));
        averaged.setCalibration(image.getCalibration());
        return averaged;
    }

    /**
     * Round a float result back into an ImageJ1 processor of the given depth.
     *
     * The arithmetic is ImageJ1's own float conversion - add a half and truncate, clamping at the
     * type's limits - because these values are compared against, and written into, files produced
     * before any of this was imglib2.
     */
    private static ImageProcessor toProcessor(RandomAccessibleInterval<FloatType> image, int bitDepth) {
        int width = (int) image.dimension(0);
        int height = (int) image.dimension(1);

        float[] raw = flatten(image, width * height);

        if (bitDepth == 32) {
            return new FloatProcessor(width, height, raw, null);
        }

        if (bitDepth == 8) {
            byte[] pixels = new byte[raw.length];
            for (int i = 0; i < raw.length; i++) {
                double value = raw[i] + 0.5;
                if (value < 0.0) { value = 0.0; }
                if (value > 255.0) { value = 255.0; }
                pixels[i] = (byte) ((int) value & 255);
            }
            return new ByteProcessor(width, height, pixels, null);
        }

        short[] pixels = new short[raw.length];
        for (int i = 0; i < raw.length; i++) {
            double value = raw[i] + 0.5;
            if (value < 0.0) { value = 0.0; }
            if (value > 65535.0) { value = 65535.0; }
            pixels[i] = (short) (int) value;
        }
        return new ShortProcessor(width, height, pixels, null);
    }

    /**
     * Flatten an image into a row major array, which is the layout ImageJ1 processors want.
     */
    private static float[] flatten(RandomAccessibleInterval<FloatType> image, int size) {
        float[] out = new float[size];
        int i = 0;
        for (FloatType value : Views.flatIterable(image)) {
            out[i++] = value.get();
        }
        return out;
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

            double[][] source = new double[3][2];
            double[][] target = new double[3][2];
            for (int i = 0; i < 3; i++) {
                source[i][0] = sourcePoints.get(i).get(0);
                source[i][1] = sourcePoints.get(i).get(1);
                target[i][0] = targetPoints.get(i).get(0);
                target[i][1] = targetPoints.get(i).get(1);
            }

            double[] affine = solveAffine(source, target);
            AffineTransform2D model = new AffineTransform2D();
            model.set(affine[0], affine[1], affine[2], affine[3], affine[4], affine[5]);
            transformModel = model;

        } catch (Exception e) {
            log.info(e);
            log.error(e);
        }
    }

    /**
     * The affine carrying the three source landmarks onto the three target landmarks.
     *
     * Three pairs determine an affine exactly, so this is a solve and not a fit - two 3x3 systems
     * sharing a matrix, by Cramer's rule. mpicbg's AffineModel2D.fit was doing the same thing for
     * the same three points; the results agreed to 2e-13 pixels, and this drops the dependency.
     */
    private static double[] solveAffine(double[][] source, double[][] target) {
        double[][] m = {
            {source[0][0], source[0][1], 1.0},
            {source[1][0], source[1][1], 1.0},
            {source[2][0], source[2][1], 1.0}};

        double det = determinant(m);
        if (Math.abs(det) < 1.0e-12) {
            throw new smFRETAnalysisException("Error: mapping landmarks are collinear, no affine exists");
        }

        double[] affine = new double[6];
        for (int axis = 0; axis < 2; axis++) {
            double[] rhs = {target[0][axis], target[1][axis], target[2][axis]};
            for (int column = 0; column < 3; column++) {
                double[][] substituted = new double[3][];
                for (int row = 0; row < 3; row++) {
                    substituted[row] = m[row].clone();
                    substituted[row][column] = rhs[row];
                }
                affine[axis * 3 + column] = determinant(substituted) / det;
            }
        }
        return affine;
    }

    /**
     * solveAffine helper.
     */
    private static double determinant(double[][] m) {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
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
        int height = image.getHeight();

        // Views on the two halves. Note the halves do not end up the same type: the target keeps
        // whatever the input was, while a transformed source comes back as 16 bit. That asymmetry
        // predates imglib2 and is preserved here rather than quietly fixed.
        Img<FloatType> whole = ImageJFunctions.convertFloat(image);
        RandomAccessibleInterval<FloatType> targetHalf =
                Views.zeroMin(Views.interval(whole, new long[]{0, 0}, new long[]{hw - 1, height - 1}));
        RandomAccessibleInterval<FloatType> sourceHalf =
                Views.zeroMin(Views.interval(whole, new long[]{hw, 0}, new long[]{image.getWidth() - 1, height - 1}));

        ImagePlus imageTarget = new ImagePlus("target", toProcessor(targetHalf, image.getBitDepth()));

        // Transform source if requested.
        ImagePlus imageSource;
        if (transform) {
            imageSource = transformImagePlus(sourceHalf, hw, height);
        }
        else {
            imageSource = new ImagePlus("source", toProcessor(sourceHalf, image.getBitDepth()));
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
    private String saveTempImageFile(ImagePlus sourceImg) throws java.io.IOException {
        java.nio.file.Path directory = Files.createTempDirectory("smfret-");
        directory.toFile().deleteOnExit();
        String sourcePathAndFileName = directory.resolve(sourceImg.getTitle() + ".tif").toString();
        new FileSaver(sourceImg).saveAsTiff(sourcePathAndFileName);
        new File(sourcePathAndFileName).deleteOnExit();
        return sourcePathAndFileName;
    }

    /**
     * Transform an ImagePlus image with the current transform.
     *
     * This was TurboReg's -transform through a temporary file, then mpicbg in memory, and is now
     * imglib2. Each move was measured against the one before it.
     *
     * Against TurboReg the difference is real but small, and it is interpolation order rather
     * than geometry: TurboReg uses a cubic spline where this is n-linear, spot centroids agree to
     * 0.003 px, the stored 16 bit values differ for 14% of pixels and almost always by exactly
     * 1 ADU, and per-spot intensities after the smoothing that precedes trace measurement differ
     * by a median of 0.45%. Rounding into 16 bit below has an rms of 0.26 ADU against 0.24 ADU
     * between the interpolators, so storing the image perturbs it more than the interpolator
     * choice does.
     *
     * Against mpicbg there is no difference at all: n-linear here and bilinear there agreed on
     * every pixel of the example data, which is why this move needed no re-validation of anything
     * downstream.
     *
     * Out of range pixels are zero, which is what createOverlapMask relies on to find the region
     * the two channels share. The test is on the inverse mapped coordinate, not on the sampled
     * value: a target pixel is written only when it lands strictly inside the source, and is left
     * at zero otherwise. That is what mpicbg's mapInterpolated did, and reproducing it is the
     * difference between this migration changing nothing and changing a 3 pixel band around the
     * frame - which was measured, and reached up to 19 ADU before this was matched. Extending
     * with zero and letting the interpolator blend across the boundary would be defensible on its
     * own terms, but it is a change of behaviour and not of implementation.
     */
    private ImagePlus transformImagePlus(RandomAccessibleInterval<FloatType> image, int width, int height){
        if (transformModel == null) {
            throw new smFRETAnalysisException("Error: Cannot transform image, transform model not set");
        }

        // Interpolate in double and store the result as float. ImageJ1's bilinear interpolation,
        // which mpicbg called and which produced every result this pipeline has been validated
        // against, accumulates in double and lands in a FloatProcessor. Interpolating in float
        // instead disagrees by one ADU on about five pixels per million - sparse, but enough to
        // make this migration show up in the traces.
        RandomAccessibleInterval<DoubleType> precise =
                Converters.convert(image, (in, out) -> out.set(in.getRealDouble()), new DoubleType());

        // extendBorder is never actually reached, because of the bounds test below. It is here so
        // that the interpolator has something defined to read at the last row and column.
        RealRandomAccess<DoubleType> source = Views.interpolate(
                Views.extendBorder(precise), new NLinearInterpolatorFactory<DoubleType>()).realRandomAccess();

        Img<FloatType> transformed = ArrayImgs.floats(width, height);
        Cursor<FloatType> target = transformed.localizingCursor();
        double[] targetPosition = new double[2];
        double[] sourcePosition = new double[2];

        while (target.hasNext()) {
            target.fwd();
            targetPosition[0] = target.getDoublePosition(0);
            targetPosition[1] = target.getDoublePosition(1);
            transformModel.applyInverse(sourcePosition, targetPosition);

            if ((sourcePosition[0] < 0.0) || (sourcePosition[0] > (image.dimension(0) - 1))
                    || (sourcePosition[1] < 0.0) || (sourcePosition[1] > (image.dimension(1) - 1))) {
                target.get().set(0.0f);
                continue;
            }

            source.setPosition(sourcePosition);
            target.get().set((float) source.get().get());
        }

        // convert to 16 bit.
        return new ImagePlus("transformed image", toProcessor(transformed, 16));
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
            ImagePlus transformedSource = transformImagePlus(
                    ImageJFunctions.convertFloat(averageImageSource),
                    averageImageSource.getWidth(), averageImageSource.getHeight());

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
            log.error(e);
        }
    }
}
