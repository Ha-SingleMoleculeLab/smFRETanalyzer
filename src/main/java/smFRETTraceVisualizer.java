/*
 * This class lets a single spot's time trace be inspected alongside the images it was measured
 * from.
 *
 * Like smFRETTraceHistogram it is interactive rather than batch and writes nothing. It differs
 * in using two windows: the field the spots were found in, which is where a spot is picked, and
 * a second window holding everything that follows from that pick - the traces above, and the
 * donor and acceptor channels zoomed in on the spot side by side below. The frame the two
 * images show is chosen with a slider, so they can be walked through while the traces stay in
 * view and a cursor tracks the position along them.
 *
 * Its input is the spot finder JSON, which names everything else - the image, the mapping, the
 * spot table - and whose own path gives the root name the trace H5 was written under.
 */

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

import com.fasterxml.jackson.databind.ObjectMapper;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashMap;
import java.util.Map;


// The single-type import of org.scijava.plugin.Menu shadows java.awt.Menu from the wildcard
// import above, so Menu here is the SciJava annotation.
@Plugin(type = Command.class,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Trace Visualizer", weight = 5.0)})
public class smFRETTraceVisualizer implements Command {

    @Parameter
    LogService log;

    @Parameter(description = "JSON file written by smFRET Spot Finder", label = "Spot Finder JSON file", style = "open")
    File spotJSONFile;

    // The suffix smFRETSpotFinder appends to its root name when it writes the JSON. Stripping it
    // recovers the root, which is what every other output is named from.
    private static final String JSON_SUFFIX = "_spotf_finding.json";

    // FRET is drawn over a fixed range, matching smFRETTraceHistogram so that the two views of
    // the same data agree.
    private static final double FRET_MIN = -0.2;
    private static final double FRET_MAX = 1.2;

    // Display range for the zoom panels is measured from this many frames, evenly spaced. Every
    // sample costs a split and a warp of one frame, so measuring all of them would stall a click
    // for a second or so on a long movie for a range that barely moves.
    private static final int ZOOM_RANGE_SAMPLES = 48;

    private static final Color ACCEPTOR_COLOR = new Color(200, 60, 40);
    private static final Color DONOR_COLOR = new Color(30, 140, 60);
    private static final Color FRET_COLOR = new Color(70, 115, 175);

    // Member variables.
    private ZoomPanel acceptorPanel;
    private java.util.List<ImagePlus> cachedSplit;
    private int cachedSplitFrame = -1;
    private int currentFrame = 1;
    private ZoomPanel donorPanel;
    private ImagePlus fieldImage;
    private FieldPanel fieldPanel;
    private JSlider frameSlider;
    private final java.util.List<JFrame> frames = new java.util.ArrayList<>();
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private ImagePlus movie;
    private int nFrames = 0;
    private int nSpots = 0;
    private int selectedSpot = -1;
    private final smFRETSpotFinder smfsf = new smFRETSpotFinder();
    private float[][] sourceTraces;      // [spot][frame], acceptor.
    private double[][] spots;            // [spot][x, y, snr, prominence] - the reloaded layout.
    private double spotSigma = 2.0;
    private JLabel statusLabel;
    private float[][] targetTraces;      // [spot][frame], donor.
    private TracePanel tracePanel;
    private int zoomHalfWidth = 8;
    private double zoomHigh = 1.0;
    private double zoomLow = 0.0;

    /**
     * Everything the plugin needs, found from the spot finder JSON.
     *
     * The root name recorded inside the JSON is relative to whatever directory the spot finder
     * ran in, so it does not survive the files being opened from anywhere else. The JSON's own
     * path does, and it was written as root + JSON_SUFFIX, so the root is recovered from that
     * and the recorded value is only a fallback.
     */
    void load(File jsonFile) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        Map<String, Object> mapping = mapper.readValue(jsonFile, HashMap.class);

        String path = jsonFile.getAbsolutePath();
        String root = path.endsWith(JSON_SUFFIX)
                ? path.substring(0, path.length() - JSON_SUFFIX.length())
                : (String) mapping.get("root name");
        File jsonDir = jsonFile.getAbsoluteFile().getParentFile();

        spotSigma = ((Number) mapping.get("spot sigma")).doubleValue();
        zoomHalfWidth = Math.max(6, (int) Math.round(4.0 * spotSigma));

        File spotsFile = locate((String) mapping.get("spots file"), jsonDir, root + "_spotf_spots.csv");
        File imageFile = locate((String) mapping.get("image name"), jsonDir, null);
        File mappingFile = locate((String) mapping.get("mapping file"), jsonDir, null);
        File fieldFile = new File(root + "_spotf_qc_image.tif");
        File traceFile = new File(root + ".h5");

        if (!traceFile.exists()) {
            throw new smFRETAnalysisException("No traces at " + traceFile
                    + " - run smFRET Time Traces on this spot finder JSON first.");
        }
        for (File needed : new File[] {spotsFile, imageFile, mappingFile, fieldFile}) {
            if ((needed == null) || !needed.exists()) {
                throw new smFRETAnalysisException("Could not find " + needed
                        + ", named by " + jsonFile);
            }
        }

        smfsf.log = log;
        smfsf.loadMappingJSON(mappingFile.toString());
        spots = smfsf.loadSpotLocations(spotsFile.toString());
        nSpots = spots.length;

        try (IHDF5Reader reader = HDF5Factory.openForReading(traceFile)) {
            targetTraces = reader.readFloatMatrix("target-traces");
            sourceTraces = reader.readFloatMatrix("source-traces");
        }
        if (targetTraces.length != nSpots) {

            // Row j of the trace matrices is row j of the spot table - measureTimeTraces walks
            // the spots in the order it loaded them - so a mismatch means the two files are from
            // different runs and every trace would be attributed to the wrong spot.
            throw new smFRETAnalysisException("The spot table has " + nSpots + " spots but "
                    + traceFile + " has " + targetTraces.length
                    + " traces - they are from different runs.");
        }
        nFrames = Math.min(targetTraces[0].length, sourceTraces[0].length);

        fieldImage = new ImagePlus(fieldFile.toString());
        movie = new ImagePlus(imageFile.toString());
        smFRETChannelMapper.requireGrayscale(movie, "the image " + imageFile);
        smFRETChannelMapper.requireTimeStack(movie, "the image " + imageFile);

        log.info("loaded " + nSpots + " spots, " + nFrames + " frames, from " + root);
    }

    /**
     * The recorded path if it is still there, otherwise the same file beside the JSON.
     */
    private static File locate(String recorded, File jsonDir, String derived) {
        if (recorded != null) {
            File asRecorded = new File(recorded);
            if (asRecorded.exists()) {
                return asRecorded;
            }
            File beside = new File(jsonDir, new File(recorded).getName());
            if (beside.exists()) {
                return beside;
            }
        }
        return (derived == null) ? new File(String.valueOf(recorded)) : new File(derived);
    }

    /**
     * The donor and acceptor halves of one frame, the acceptor warped onto the donor's frame.
     *
     * Only the frame being displayed is converted, rather than the movie being converted up
     * front the way smFRETAnalyzer does it - a viewer that turned a 1295 frame movie into float
     * on open would want gigabytes to show one spot.
     */
    private java.util.List<ImagePlus> splitFrame(int frame) {
        if ((frame == cachedSplitFrame) && (cachedSplit != null)) {
            return cachedSplit;
        }
        ImageProcessor slice = movie.getStack().getProcessor(frame);
        cachedSplit = smfsf.splitImagePlus(new ImagePlus("frame", slice));
        cachedSplitFrame = frame;
        return cachedSplit;
    }

    /**
     * The square of pixels around a spot that the zoom panels show.
     */
    private float[] crop(ImagePlus half, int centreX, int centreY) {
        int size = 2 * zoomHalfWidth + 1;
        float[] out = new float[size * size];

        // getf, not getPixel: on a float image getPixel hands back the raw bits of the float.
        ImageProcessor processor = half.getProcessor();
        for (int dy = 0; dy < size; dy++) {
            for (int dx = 0; dx < size; dx++) {
                int x = centreX - zoomHalfWidth + dx;
                int y = centreY - zoomHalfWidth + dy;
                boolean inside = (x >= 0) && (y >= 0)
                        && (x < processor.getWidth()) && (y < processor.getHeight());
                out[dy * size + dx] = inside ? processor.getf(x, y) : 0.0f;
            }
        }
        return out;
    }

    /**
     * One display range for both zoom panels over the whole movie, so that brightness means the
     * same thing from frame to frame and from one channel to the other - which is what makes
     * bleaching and FRET anticorrelation visible rather than being normalized away.
     */
    private void measureZoomRange() {
        int x = (int) spots[selectedSpot][0];
        int y = (int) spots[selectedSpot][1];
        int samples = Math.min(nFrames, ZOOM_RANGE_SAMPLES);

        double low = Double.MAX_VALUE;
        double high = -Double.MAX_VALUE;
        for (int i = 0; i < samples; i++) {
            int frame = 1 + (int) Math.round((double) i * (nFrames - 1) / Math.max(1, samples - 1));
            for (ImagePlus half : splitFrame(frame)) {
                for (float value : crop(half, x, y)) {
                    if (value < low) { low = value; }
                    if (value > high) { high = value; }
                }
            }
        }
        if (high <= low) {
            high = low + 1.0;
        }
        zoomLow = low;
        zoomHigh = high;
    }

    /**
     * Grey levels from float pixels, everything outside the range flattened to black or white.
     */
    private static BufferedImage render(float[] pixels, int width, int height, double low, double high) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        double scale = 255.0 / (high - low);
        for (int i = 0; i < pixels.length; i++) {
            int level = (int) Math.round((pixels[i] - low) * scale);
            level = Math.max(0, Math.min(255, level));
            image.setRGB(i % width, i / width, (level << 16) | (level << 8) | level);
        }
        return image;
    }

    /**
     * A display range that leaves the spots visible.
     *
     * Plain minimum to maximum is no use on a field of spots: one bright spot sets the top and
     * everything else collapses into the bottom few levels. These are percentiles instead, which
     * is what ImageJ's own auto contrast does.
     */
    private static double[] displayRange(float[] pixels, double lowFraction, double highFraction) {
        float[] sorted = pixels.clone();
        java.util.Arrays.sort(sorted);
        double low = sorted[(int) Math.min(sorted.length - 1, lowFraction * sorted.length)];
        double high = sorted[(int) Math.min(sorted.length - 1, highFraction * sorted.length)];
        if (high <= low) {
            high = low + 1.0;
        }
        return new double[] {low, high};
    }

    /**
     * The spots window: the image the spots were found in, with the spots on it.
     *
     * This is the spot finder's QC image rather than a frame of the movie, so it is the averaged
     * image at whatever spotChannel was chosen, and the spots sit where they were actually
     * found. A single frame would be far noisier and most spots would not be visible at all.
     */
    private class FieldPanel extends JPanel {

        private final BufferedImage image;
        private double scale = 1.0;
        private int offsetX = 0;
        private int offsetY = 0;

        FieldPanel() {
            ImageProcessor processor = fieldImage.getProcessor();
            int width = processor.getWidth();
            int height = processor.getHeight();
            float[] pixels = new float[width * height];
            for (int i = 0; i < pixels.length; i++) {
                pixels[i] = processor.getf(i % width, i / width);
            }
            double[] range = displayRange(pixels, 0.005, 0.9995);
            image = render(pixels, width, height, range[0], range[1]);

            setPreferredSize(new Dimension(Math.min(560, width * 2), Math.min(760, height * 2)));
            setBackground(Color.BLACK);
            addMouseListener(new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    selectNearest((e.getX() - offsetX) / scale, (e.getY() - offsetY) / scale);
                }
            });
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g.create();

            scale = Math.min((double) getWidth() / image.getWidth(),
                    (double) getHeight() / image.getHeight());
            int drawWidth = (int) Math.round(image.getWidth() * scale);
            int drawHeight = (int) Math.round(image.getHeight() * scale);
            offsetX = (getWidth() - drawWidth) / 2;
            offsetY = (getHeight() - drawHeight) / 2;
            g2.drawImage(image, offsetX, offsetY, drawWidth, drawHeight, null);

            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            int radius = Math.max(3, (int) Math.round(2.0 * spotSigma * scale));
            for (int i = 0; i < nSpots; i++) {
                int x = offsetX + (int) Math.round(spots[i][0] * scale);
                int y = offsetY + (int) Math.round(spots[i][1] * scale);
                if (i == selectedSpot) {
                    g2.setColor(Color.YELLOW);
                    g2.setStroke(new BasicStroke(2.0f));
                    g2.drawOval(x - radius - 2, y - radius - 2, 2 * radius + 4, 2 * radius + 4);
                } else {
                    g2.setColor(new Color(70, 200, 120, 170));
                    g2.setStroke(new BasicStroke(1.0f));
                    g2.drawOval(x - radius, y - radius, 2 * radius, 2 * radius);
                }
            }
            g2.dispose();
        }
    }

    /**
     * Upper panel of the traces window: the selected spot's traces.
     *
     * Two panels sharing the frame axis rather than one panel with two vertical scales - the
     * intensities and the ratio have nothing to do with each other numerically, and FRET is
     * worth reading off to a couple of decimals, which it cannot be if it is squeezed against
     * two intensity traces.
     */
    private class TracePanel extends JPanel {

        private static final int GAP = 26;
        private static final int MARGIN_BOTTOM = 40;
        private static final int MARGIN_LEFT = 68;
        private static final int MARGIN_RIGHT = 14;
        private static final int MARGIN_TOP = 14;

        TracePanel() {
            setPreferredSize(new Dimension(640, 420));
            setBackground(Color.WHITE);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (selectedSpot < 0) {
                g.setColor(Color.GRAY);
                g.drawString("Click a spot in the field window", 20, 30);
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int plotWidth = getWidth() - MARGIN_LEFT - MARGIN_RIGHT;
            int total = getHeight() - MARGIN_TOP - MARGIN_BOTTOM - GAP;
            if ((plotWidth < 20) || (total < 40)) {
                g2.dispose();
                return;
            }

            // The intensities get the larger share: they carry two traces and their range is not
            // known in advance, where FRET is always inside a fixed band.
            int upperHeight = (int) Math.round(total * 0.62);
            int lowerHeight = total - upperHeight;
            int upperTop = MARGIN_TOP;
            int lowerTop = MARGIN_TOP + upperHeight + GAP;

            float[] donor = targetTraces[selectedSpot];
            float[] acceptor = sourceTraces[selectedSpot];

            double low = Double.MAX_VALUE;
            double high = -Double.MAX_VALUE;
            for (int t = 0; t < nFrames; t++) {
                low = Math.min(low, Math.min(donor[t], acceptor[t]));
                high = Math.max(high, Math.max(donor[t], acceptor[t]));
            }
            double pad = 0.05 * Math.max(1.0, high - low);
            low -= pad;
            high += pad;

            drawAxes(g2, upperTop, upperHeight, plotWidth, low, high, "intensity");
            drawTrace(g2, donor, upperTop, upperHeight, plotWidth, low, high, DONOR_COLOR);
            drawTrace(g2, acceptor, upperTop, upperHeight, plotWidth, low, high, ACCEPTOR_COLOR);

            float[] fret = new float[nFrames];
            for (int t = 0; t < nFrames; t++) {
                double sum = (double) donor[t] + (double) acceptor[t];

                // A near zero total makes the ratio meaningless rather than merely noisy. Parking
                // it at the bottom of the axis keeps the line continuous and obviously invalid.
                fret[t] = (Math.abs(sum) < 1.0e-9) ? (float) FRET_MIN : (float) (acceptor[t] / sum);
            }
            drawAxes(g2, lowerTop, lowerHeight, plotWidth, FRET_MIN, FRET_MAX, "FRET");
            drawTrace(g2, fret, lowerTop, lowerHeight, plotWidth, FRET_MIN, FRET_MAX, FRET_COLOR);

            // The frame the two zoom panels are showing, marked on both plots so the slider
            // position can be read against the trace rather than off the slider alone.
            int cursorX = MARGIN_LEFT + (int) Math.round((currentFrame - 1.0) * plotWidth
                    / Math.max(1, nFrames - 1));
            g2.setColor(new Color(230, 120, 20));
            g2.setStroke(new BasicStroke(1.5f));
            g2.drawLine(cursorX, upperTop, cursorX, upperTop + upperHeight);
            g2.drawLine(cursorX, lowerTop, cursorX, lowerTop + lowerHeight);

            // Legend.
            FontMetrics fm = g2.getFontMetrics();
            int legendX = MARGIN_LEFT + 6;
            g2.setColor(DONOR_COLOR);
            g2.drawString("donor", legendX, upperTop + fm.getAscent() + 2);
            g2.setColor(ACCEPTOR_COLOR);
            g2.drawString("acceptor", legendX + fm.stringWidth("donor") + 12,
                    upperTop + fm.getAscent() + 2);

            g2.setColor(Color.DARK_GRAY);
            String xTitle = "frame";
            g2.drawString(xTitle, MARGIN_LEFT + (plotWidth - fm.stringWidth(xTitle)) / 2,
                    getHeight() - 8);
            g2.dispose();
        }

        private void drawAxes(Graphics2D g2, int top, int height, int plotWidth,
                              double low, double high, String title) {
            g2.setColor(Color.DARK_GRAY);
            g2.setStroke(new BasicStroke(1.0f));
            g2.drawLine(MARGIN_LEFT, top + height, MARGIN_LEFT + plotWidth, top + height);
            g2.drawLine(MARGIN_LEFT, top, MARGIN_LEFT, top + height);

            FontMetrics fm = g2.getFontMetrics();
            for (int i = 0; i <= 2; i++) {
                double value = low + (high - low) * i / 2.0;
                int y = top + height - (int) Math.round((double) i * height / 2.0);
                g2.drawLine(MARGIN_LEFT - 4, y, MARGIN_LEFT, y);
                String label = (Math.abs(high - low) < 10.0)
                        ? String.format("%.2f", value) : String.format("%.0f", value);
                g2.drawString(label, MARGIN_LEFT - 8 - fm.stringWidth(label), y + fm.getAscent() / 2 - 1);
            }

            // Frame ticks on the shared axis.
            for (int i = 0; i <= 4; i++) {
                int frame = 1 + (int) Math.round((double) i * (nFrames - 1) / 4.0);
                int x = MARGIN_LEFT + (int) Math.round((double) i * plotWidth / 4.0);
                g2.drawLine(x, top + height, x, top + height + 4);
                String label = Integer.toString(frame);
                g2.drawString(label, x - fm.stringWidth(label) / 2, top + height + 16);
            }

            Graphics2D g2r = (Graphics2D) g2.create();
            g2r.rotate(-Math.PI / 2.0, 14, top + height / 2.0);
            g2r.drawString(title, 14 - fm.stringWidth(title) / 2, top + height / 2.0f);
            g2r.dispose();
        }

        private void drawTrace(Graphics2D g2, float[] values, int top, int height, int plotWidth,
                               double low, double high, Color color) {
            g2.setColor(color);
            g2.setStroke(new BasicStroke(1.2f));
            int previousX = 0;
            int previousY = 0;
            for (int t = 0; t < nFrames; t++) {
                int x = MARGIN_LEFT + (int) Math.round((double) t * plotWidth / Math.max(1, nFrames - 1));
                double fraction = (values[t] - low) / (high - low);
                fraction = Math.max(0.0, Math.min(1.0, fraction));
                int y = top + height - (int) Math.round(fraction * height);
                if (t > 0) {
                    g2.drawLine(previousX, previousY, x, y);
                }
                previousX = x;
                previousY = y;
            }
        }
    }

    /**
     * Lower panels of the traces window: one channel of the selected spot at one frame.
     */
    private class ZoomPanel extends JPanel {

        private final String title;
        private float[] pixels;

        ZoomPanel(String title) {
            this.title = title;
            setPreferredSize(new Dimension(260, 280));
            setBackground(Color.BLACK);
        }

        void setPixels(float[] pixels) {
            this.pixels = pixels;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (pixels == null) {
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            int size = 2 * zoomHalfWidth + 1;
            BufferedImage image = render(pixels, size, size, zoomLow, zoomHigh);

            int side = Math.min(getWidth(), getHeight() - 22);
            int offsetX = (getWidth() - side) / 2;
            int offsetY = 22 + (getHeight() - 22 - side) / 2;

            // Nearest neighbour: at this magnification the pixels are the data, and smoothing
            // them would invent detail the camera never recorded.
            g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
            g2.drawImage(image, offsetX, offsetY, side, side, null);

            // Where the spot was measured, which is the centre pixel by construction.
            double pixelSide = (double) side / size;
            int centre = (int) Math.round(offsetX + (zoomHalfWidth + 0.5) * pixelSide);
            int centreY = (int) Math.round(offsetY + (zoomHalfWidth + 0.5) * pixelSide);
            int marker = (int) Math.round(2.0 * spotSigma * pixelSide);
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setColor(new Color(255, 230, 40, 190));
            g2.drawOval(centre - marker, centreY - marker, 2 * marker, 2 * marker);

            g2.setColor(Color.WHITE);
            g2.drawString(title + " - frame " + currentFrame, 8, 15);
            g2.dispose();
        }
    }

    /**
     * Select the spot nearest a point in field image coordinates, if one is close enough.
     */
    private void selectNearest(double x, double y) {
        double limit = Math.max(4.0, 2.0 * spotSigma);
        double bestDistance = Double.MAX_VALUE;
        int best = -1;
        for (int i = 0; i < nSpots; i++) {
            double distance = Math.hypot(spots[i][0] - x, spots[i][1] - y);
            if ((distance < bestDistance) && (distance <= limit)) {
                bestDistance = distance;
                best = i;
            }
        }
        if (best >= 0) {
            selectSpot(best);
        }
    }

    /**
     * Show a spot: its traces, and its images at the current frame.
     */
    private void selectSpot(int index) {
        if ((index < 0) || (index >= nSpots)) {
            return;
        }
        selectedSpot = index;
        measureZoomRange();
        updateZoom();
        fieldPanel.repaint();
        tracePanel.repaint();
        updateStatus();
    }

    /**
     * Redraw the two zoom panels for the current spot and frame.
     */
    private void updateZoom() {
        if (selectedSpot < 0) {
            return;
        }
        int x = (int) spots[selectedSpot][0];
        int y = (int) spots[selectedSpot][1];
        java.util.List<ImagePlus> split = splitFrame(currentFrame);
        donorPanel.setPixels(crop(split.get(0), x, y));
        acceptorPanel.setPixels(crop(split.get(1), x, y));
    }

    private void updateStatus() {
        if (selectedSpot < 0) {
            statusLabel.setText("no spot selected");
            return;
        }
        double donor = targetTraces[selectedSpot][currentFrame - 1];
        double acceptor = sourceTraces[selectedSpot][currentFrame - 1];
        double sum = donor + acceptor;

        // A total near zero - the illumination off at the end of a movie, say - makes the ratio
        // blow up rather than merely get noisy. The plot clamps such a value to the edge of its
        // fixed range, so the number is flagged here to say why the two do not look alike.
        String fret;
        if (Math.abs(sum) < 1.0e-9) {
            fret = "n/a";
        } else {
            double value = acceptor / sum;
            fret = String.format("%.3f", value)
                    + (((value < FRET_MIN) || (value > FRET_MAX)) ? " (off scale)" : "");
        }
        statusLabel.setText(String.format(
                "spot %,d of %,d at (%d, %d) · SNR %.1f · frame %d: D %,.0f  A %,.0f  FRET %s",
                selectedSpot + 1, nSpots, (int) spots[selectedSpot][0], (int) spots[selectedSpot][1],
                spots[selectedSpot][2], currentFrame, donor, acceptor, fret));
    }

    /**
     * Move to a different frame, which moves the zoom panels and the trace cursor together.
     */
    private void setFrame(int frame) {
        currentFrame = Math.max(1, Math.min(nFrames, frame));
        updateZoom();
        tracePanel.repaint();
        updateStatus();
    }

    /**
     * Build the two windows.
     *
     * The field keeps a window of its own because it is 256 x 512 and wants the height, and
     * because clicking an individual spot in a crowded field needs it drawn large. Everything
     * that follows from a selection - the traces and the two channel images - is one window,
     * since those are always read together.
     */
    private void showWindows() {
        String name = spotJSONFile.getName();

        fieldPanel = new FieldPanel();
        JFrame fieldFrame = frame("smFRET spots - " + name, fieldPanel);

        tracePanel = new TracePanel();

        donorPanel = new ZoomPanel("donor");
        acceptorPanel = new ZoomPanel("acceptor");
        JPanel zoomRow = new JPanel(new GridLayout(1, 2, 6, 0));
        zoomRow.add(donorPanel);
        zoomRow.add(acceptorPanel);

        // A split rather than a fixed division, so that either half can be given the space when
        // one of them is what is being looked at. The weight sends resizing to the traces: the
        // images are square and stop gaining from extra height, where a longer trace does not.
        JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT, tracePanel, zoomRow);
        split.setResizeWeight(1.0);
        split.setBorder(null);

        frameSlider = new JSlider(1, Math.max(1, nFrames), 1);
        frameSlider.addChangeListener(e -> setFrame(frameSlider.getValue()));

        JButton previousSpot = new JButton("< spot");
        previousSpot.addActionListener(e -> selectSpot(selectedSpot - 1));
        JButton nextSpot = new JButton("spot >");
        nextSpot.addActionListener(e -> selectSpot(selectedSpot + 1));

        JPanel controls = new JPanel(new BorderLayout(8, 0));
        controls.setBorder(new EmptyBorder(2, 10, 6, 10));
        JPanel sliderRow = new JPanel(new BorderLayout(6, 0));
        sliderRow.add(new JLabel("Frame"), BorderLayout.WEST);
        sliderRow.add(frameSlider, BorderLayout.CENTER);
        JPanel spotButtons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 4, 0));
        spotButtons.add(previousSpot);
        spotButtons.add(nextSpot);
        sliderRow.add(spotButtons, BorderLayout.EAST);

        statusLabel = new JLabel(" ");
        controls.add(sliderRow, BorderLayout.NORTH);
        controls.add(statusLabel, BorderLayout.SOUTH);

        JPanel traceContent = new JPanel(new BorderLayout());
        traceContent.add(split, BorderLayout.CENTER);
        traceContent.add(controls, BorderLayout.SOUTH);

        JFrame traceFrame = frame("smFRET traces - " + name, traceContent);

        // Either window closes both. Neither is any use on its own - the field has nothing to
        // report a selection to, and the traces have no way to change which spot they show.
        WindowAdapter closeBoth = new WindowAdapter() {
            @Override
            public void windowClosed(WindowEvent e) {
                for (JFrame other : frames) {
                    if (other != e.getWindow()) {
                        other.dispose();
                    }
                }
            }
        };
        fieldFrame.addWindowListener(closeBoth);
        traceFrame.addWindowListener(closeBoth);

        Rectangle screen = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        int left = screen.x + 20;
        int top = screen.y + 20;
        fieldFrame.setLocation(left, top);
        traceFrame.setLocation(left + fieldFrame.getWidth() + 12, top);

        for (JFrame each : frames) {
            each.setVisible(true);
        }

        updateStatus();
        if (nSpots > 0) {
            selectSpot(0);
        }
    }

    /**
     * showWindows() helper, one packed frame registered for group disposal.
     */
    private JFrame frame(String title, JComponent content) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.getContentPane().add(content);
        frame.pack();
        frames.add(frame);
        return frame;
    }

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
            if (isHeadless) {
                log.info("smFRET Trace Visualizer is interactive and cannot run headless");
                return;
            }

            log.info("loading " + spotJSONFile);
            load(spotJSONFile);

            SwingUtilities.invokeLater(this::showWindows);

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
