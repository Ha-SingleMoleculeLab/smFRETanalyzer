/*
 * This class plots histograms of the time traces measured by smFRETAnalyzer.
 *
 * Unlike the other plugins in this package it is interactive rather than batch - it opens a
 * window whose histogram is recomputed as the controls are adjusted, so that thresholds can
 * be chosen by eye. It reads the '.h5' file written by smFRETAnalyzer.
 */

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

import ij.IJ;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;


// The single-type import of org.scijava.plugin.Menu shadows java.awt.Menu from the wildcard
// import above, so Menu here is the SciJava annotation.
@Plugin(type = Command.class,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Trace Histograms", weight = 4.0)})
public class smFRETTraceHistogram implements Command {

    @Parameter
    LogService log;

    @Parameter(description = "H5 file written by smFRET Time Traces", label = "Trace H5 file", style = "open")
    File h5File;

    // Histogram types, indices match the order of the type radio buttons.
    private static final int TYPE_FRET = 0;
    private static final int TYPE_DONOR = 1;
    private static final int TYPE_ACCEPTOR = 2;
    private static final int TYPE_TOTAL = 3;
    private static final String[] TYPE_NAMES = {"FRET efficiency", "Donor (target)", "Acceptor (source)", "Total (D+A)"};

    // Quantity the intensity threshold is applied to. These are mutually exclusive, so they are a
    // combo box beside a single slider rather than one slider each.
    private static final int FILTER_TOTAL = 0;
    private static final int FILTER_DONOR = 1;
    private static final int FILTER_ACCEPTOR = 2;
    private static final String[] FILTER_NAMES = {"Total (D+A)", "Donor (target)", "Acceptor (source)"};

    // FRET efficiency is plotted over a fixed range, slightly wider than [0,1] so that the
    // noise skirts either side of the physical range stay visible.
    private static final double FRET_MIN = -0.2;
    private static final double FRET_MAX = 1.2;

    // Member variables.
    private JSlider binsSlider;
    private JComboBox<String> filterCombo;
    private final double[] filterMax = new double[FILTER_NAMES.length];
    private final double[] filterMin = new double[FILTER_NAMES.length];
    private JSlider firstFrameSlider;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private JSlider lastFrameSlider;
    private JSlider minValueSlider;
    private int nFrames = 0;
    private int nSpots = 0;
    private HistogramPanel plotPanel;
    private Histogram result;
    private float[][] sourceTraces;      // [spot][frame], acceptor.
    private JLabel statusLabel;
    private boolean suspendUpdates = false;
    private float[][] targetTraces;      // [spot][frame], donor.
    private JRadioButton[] typeButtons;

    /**
     * The result of binning the traces, everything the plot panel needs to draw itself.
     */
    static class Histogram {
        double binWidth;
        int[] counts;
        double lo;
        int maxCount;
        int nOutside;       // Traces dropped for falling outside [lo,hi].
        int nPoints;        // Traces actually binned.
        int nSpotsUsed;     // Traces passing the intensity threshold.
        String valueLabel;
    }

    /**
     * Bin the loaded traces. Takes its settings as arguments rather than reading the controls
     * directly so that the binning can be exercised without a GUI.
     */
    Histogram computeHistogram(int type, int firstFrame, int lastFrame,
                               int filterType, double minValue, int nBins) {

        // One point per trace, the average over the selected interval. For FRET the donor and
        // acceptor are averaged first and the ratio taken from those averages - averaging the
        // per frame ratios instead would not give the same answer.
        double[] values = new double[nSpots];
        int nValues = 0;
        int nSpotsUsed = 0;
        int nIntervalFrames = lastFrame - firstFrame + 1;

        for (int i = 0; i < nSpots; i++) {
            double donorSum = 0.0;
            double acceptorSum = 0.0;
            double lowestFrameValue = Double.MAX_VALUE;
            for (int t = firstFrame - 1; t < lastFrame; t++) {
                double frameDonor = targetTraces[i][t];
                double frameAcceptor = sourceTraces[i][t];
                donorSum += frameDonor;
                acceptorSum += frameAcceptor;

                double frameValue = filterValue(filterType, frameDonor, frameAcceptor,
                        frameDonor + frameAcceptor);
                if (frameValue < lowestFrameValue) {
                    lowestFrameValue = frameValue;
                }
            }

            // The whole trace goes if any single frame in the interval is below the threshold, so
            // a molecule that bleaches part way through contributes nothing rather than a diluted
            // average.
            if (lowestFrameValue < minValue) {
                continue;
            }

            double donor = donorSum / nIntervalFrames;
            double acceptor = acceptorSum / nIntervalFrames;
            double total = donor + acceptor;

            double value;
            if (type == TYPE_FRET) {
                // A near zero total makes the ratio meaningless, not just noisy.
                if (Math.abs(total) < 1.0e-9) {
                    continue;
                }
                value = acceptor / total;
            } else if (type == TYPE_DONOR) {
                value = donor;
            } else if (type == TYPE_ACCEPTOR) {
                value = acceptor;
            } else {
                value = total;
            }

            values[nValues++] = value;
            nSpotsUsed += 1;
        }

        Histogram hist = new Histogram();
        hist.counts = new int[nBins];
        hist.nSpotsUsed = nSpotsUsed;
        hist.valueLabel = TYPE_NAMES[type];

        // Fixed range for FRET efficiency, auto range for the intensity histograms.
        double lo;
        double hi;
        if (type == TYPE_FRET) {
            lo = FRET_MIN;
            hi = FRET_MAX;
        } else {
            lo = Double.MAX_VALUE;
            hi = -Double.MAX_VALUE;
            for (int i = 0; i < nValues; i++) {
                if (values[i] < lo) { lo = values[i]; }
                if (values[i] > hi) { hi = values[i]; }
            }
            if (nValues == 0) {
                lo = 0.0;
                hi = 1.0;
            }
        }
        if (hi <= lo) {
            hi = lo + 1.0;
        }

        hist.lo = lo;
        hist.binWidth = (hi - lo) / nBins;

        double hiEdge = lo + hist.binWidth * nBins;
        for (int i = 0; i < nValues; i++) {

            // Range test the value rather than the bin index. A cast truncates toward zero, so a
            // value just below lo gives bin 0 and would be silently folded into the first bin
            // instead of being counted as out of range.
            if ((values[i] < lo) || (values[i] > hiEdge)) {
                hist.nOutside += 1;
                continue;
            }

            int bin = (int) ((values[i] - lo) / hist.binWidth);

            // The largest value lands one past the last bin, keep it rather than dropping it.
            if (bin >= nBins) {
                bin = nBins - 1;
            }
            hist.counts[bin] += 1;
            hist.nPoints += 1;
        }

        for (int count : hist.counts) {
            if (count > hist.maxCount) {
                hist.maxCount = count;
            }
        }

        return hist;
    }

    /**
     * The intensity the threshold slider is currently applied to.
     */
    private static double filterValue(int filterType, double donor, double acceptor, double total) {
        if (filterType == FILTER_DONOR) {
            return donor;
        }
        if (filterType == FILTER_ACCEPTOR) {
            return acceptor;
        }
        return total;
    }

    /**
     * Panel that draws the current histogram.
     */
    private class HistogramPanel extends JPanel {

        private static final int MARGIN_BOTTOM = 46;
        private static final int MARGIN_LEFT = 66;
        private static final int MARGIN_RIGHT = 18;
        private static final int MARGIN_TOP = 18;

        HistogramPanel() {
            setPreferredSize(new Dimension(660, 340));
            setBackground(Color.WHITE);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (result == null) {
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int plotWidth = getWidth() - MARGIN_LEFT - MARGIN_RIGHT;
            int plotHeight = getHeight() - MARGIN_TOP - MARGIN_BOTTOM;
            if ((plotWidth < 10) || (plotHeight < 10)) {
                g2.dispose();
                return;
            }

            int nBins = result.counts.length;
            double yScale = (result.maxCount > 0) ? ((double) plotHeight / (double) result.maxCount) : 0.0;

            // Bars.
            g2.setColor(new Color(70, 115, 175));
            for (int i = 0; i < nBins; i++) {
                if (result.counts[i] == 0) {
                    continue;
                }
                int x0 = MARGIN_LEFT + (int) Math.round((double) i * plotWidth / nBins);
                int x1 = MARGIN_LEFT + (int) Math.round((double) (i + 1) * plotWidth / nBins);
                int h = (int) Math.round(result.counts[i] * yScale);
                int barWidth = Math.max(1, x1 - x0 - 1);
                g2.fillRect(x0, MARGIN_TOP + plotHeight - h, barWidth, h);
            }

            // Axes.
            g2.setColor(Color.DARK_GRAY);
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP + plotHeight, MARGIN_LEFT + plotWidth, MARGIN_TOP + plotHeight);
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + plotHeight);

            FontMetrics fm = g2.getFontMetrics();

            // X ticks, placed on round multiples of a step rather than at fixed fractions of the
            // axis. For the fixed FRET range this works out as -0.2, 0.0, 0.2 ... 1.2.
            double hi = result.lo + result.binWidth * nBins;
            double range = hi - result.lo;
            double step = niceTickStep(range);
            int decimals = tickDecimals(step);
            double eps = step * 1.0e-6;

            for (int k = (int) Math.ceil(result.lo / step - eps); ; k++) {
                double value = k * step;
                if (value > (hi + eps)) {
                    break;
                }

                // Multiplying out can leave a tiny residue where the tick should be exactly zero.
                if (Math.abs(value) < eps) {
                    value = 0.0;
                }

                int x = MARGIN_LEFT + (int) Math.round((value - result.lo) / range * plotWidth);
                g2.drawLine(x, MARGIN_TOP + plotHeight, x, MARGIN_TOP + plotHeight + 4);
                String label = formatTick(value, decimals);
                g2.drawString(label, x - fm.stringWidth(label) / 2, MARGIN_TOP + plotHeight + 18);
            }

            // Y ticks.
            for (int i = 0; i <= 4; i++) {
                double frac = i / 4.0;
                int y = MARGIN_TOP + plotHeight - (int) Math.round(frac * plotHeight);
                g2.drawLine(MARGIN_LEFT - 4, y, MARGIN_LEFT, y);
                String label = Integer.toString((int) Math.round(frac * result.maxCount));
                g2.drawString(label, MARGIN_LEFT - 8 - fm.stringWidth(label), y + fm.getAscent() / 2 - 1);
            }

            // Axis titles.
            String xTitle = result.valueLabel;
            g2.drawString(xTitle, MARGIN_LEFT + (plotWidth - fm.stringWidth(xTitle)) / 2, getHeight() - 8);

            Graphics2D g2r = (Graphics2D) g2.create();
            g2r.rotate(-Math.PI / 2.0, 16, MARGIN_TOP + plotHeight / 2.0);
            g2r.drawString("counts", 16 - fm.stringWidth("counts") / 2, MARGIN_TOP + plotHeight / 2.0f);
            g2r.dispose();

            g2.dispose();
        }
    }

    /**
     * A round tick spacing (1, 2, 2.5 or 5 times a power of ten) for an axis of this range. The
     * FRET range of 1.4 gives 0.2.
     */
    private static double niceTickStep(double range) {
        if (!(range > 0.0)) {
            return 1.0;
        }

        double raw = range / 7.0;
        double magnitude = Math.pow(10.0, Math.floor(Math.log10(raw)));
        double normalized = raw / magnitude;

        double step;
        if (normalized <= 1.0) {
            step = 1.0;
        } else if (normalized <= 2.0) {
            step = 2.0;
        } else if (normalized <= 2.5) {
            step = 2.5;
        } else if (normalized <= 5.0) {
            step = 5.0;
        } else {
            step = 10.0;
        }
        return step * magnitude;
    }

    /**
     * How many decimals a tick label needs to distinguish one step from the next.
     */
    private static int tickDecimals(double step) {
        if (step >= 1.0) {
            return 0;
        }
        if (step >= 0.1) {
            return 1;
        }
        if (step >= 0.01) {
            return 2;
        }
        return 3;
    }

    /**
     * Short tick labels, the intensity histograms can run to large values.
     */
    private static String formatTick(double value, int decimals) {
        if (Math.abs(value) >= 100000.0) {
            return String.format("%.1e", value);
        }
        return String.format("%." + decimals + "f", value);
    }

    /**
     * Load the trace matrices from an smFRETAnalyzer H5 file.
     */
    void loadTraces(File file) {
        try (IHDF5Reader reader = HDF5Factory.openForReading(file)) {
            targetTraces = reader.readFloatMatrix("target-traces");
            sourceTraces = reader.readFloatMatrix("source-traces");
        }

        if ((targetTraces.length == 0) || (sourceTraces.length == 0)) {
            throw new smFRETAnalysisException("No traces in " + file);
        }
        if (targetTraces.length != sourceTraces.length) {
            throw new smFRETAnalysisException("Target and source trace counts differ ("
                    + targetTraces.length + " vs " + sourceTraces.length + ") in " + file);
        }

        nSpots = targetTraces.length;
        nFrames = Math.min(targetTraces[0].length, sourceTraces[0].length);

        // Range of each quantity the threshold can be applied to, this sets the slider limits.
        for (int f = 0; f < FILTER_NAMES.length; f++) {
            filterMin[f] = Double.MAX_VALUE;
            filterMax[f] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < nSpots; i++) {
            for (int t = 0; t < nFrames; t++) {
                double donor = targetTraces[i][t];
                double acceptor = sourceTraces[i][t];
                for (int f = 0; f < FILTER_NAMES.length; f++) {
                    double value = filterValue(f, donor, acceptor, donor + acceptor);
                    if (value < filterMin[f]) { filterMin[f] = value; }
                    if (value > filterMax[f]) { filterMax[f] = value; }
                }
            }
        }
        for (int f = 0; f < FILTER_NAMES.length; f++) {
            if (filterMax[f] <= filterMin[f]) {
                filterMax[f] = filterMin[f] + 1.0;
            }
        }

        log.info("loaded " + nSpots + " traces of " + nFrames + " frames from " + file);
    }

    /**
     * Prompt for a different H5 file and reload.
     */
    private void onBrowse(JFrame frame) {
        JFileChooser chooser = new JFileChooser(h5File.getParentFile());
        chooser.setDialogTitle("Select an smFRET trace H5 file");
        if (chooser.showOpenDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try {
            File selected = chooser.getSelectedFile();
            loadTraces(selected);
            h5File = selected;
            frame.setTitle("smFRET Trace Histograms - " + h5File.getName());
            resetSliderRanges();
            update();
        } catch (Exception e) {
            log.info(e);
            JOptionPane.showMessageDialog(frame, "Could not read traces:\n" + e.getMessage(),
                    "smFRET Trace Histograms", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * Write the current histogram out as a CSV table.
     */
    private void onSaveCsv(JFrame frame) {
        JFileChooser chooser = new JFileChooser(h5File.getParentFile());
        chooser.setSelectedFile(new File(stripExtension(h5File) + "_histogram.csv"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try (PrintWriter writer = new PrintWriter(chooser.getSelectedFile())) {
            writer.println("# " + result.valueLabel + " from " + h5File);
            writer.println("# frames " + firstFrameSlider.getValue() + "-" + lastFrameSlider.getValue()
                    + ", min " + FILTER_NAMES[filterCombo.getSelectedIndex()] + " " + minValueSlider.getValue()
                    + ", " + result.nPoints + " of " + result.nSpotsUsed + " traces in range");
            writer.println("bin_center,count");
            for (int i = 0; i < result.counts.length; i++) {
                writer.println((result.lo + (i + 0.5) * result.binWidth) + "," + result.counts[i]);
            }
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Write the current plot out as a PNG.
     */
    private void onSavePng(JFrame frame) {
        JFileChooser chooser = new JFileChooser(h5File.getParentFile());
        chooser.setSelectedFile(new File(stripExtension(h5File) + "_histogram.png"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try {
            ImageIO.write(renderPlotImage(), "png", chooser.getSelectedFile());
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Render the current plot for saving, titled with the H5 file name.
     *
     * The window itself shows the file name in its header row, but a saved plot travels on its
     * own and would otherwise lose track of which data it came from, so the title is added here
     * rather than being drawn into the panel on screen.
     */
    BufferedImage renderPlotImage() {
        String title = h5File.getName();
        int titleHeight = 30;

        BufferedImage image = new BufferedImage(plotPanel.getWidth(),
                plotPanel.getHeight() + titleHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = image.createGraphics();
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        g2.setColor(plotPanel.getBackground());
        g2.fillRect(0, 0, image.getWidth(), titleHeight);

        g2.setColor(Color.DARK_GRAY);
        g2.setFont(g2.getFont().deriveFont(Font.BOLD, 14.0f));
        FontMetrics fm = g2.getFontMetrics();
        g2.drawString(title, (image.getWidth() - fm.stringWidth(title)) / 2,
                (titleHeight + fm.getAscent()) / 2 - 2);

        g2.translate(0, titleHeight);
        plotPanel.paint(g2);
        g2.dispose();

        return image;
    }

    /**
     * Point the frame and threshold sliders at the currently loaded traces.
     */
    private void resetSliderRanges() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            firstFrameSlider.setMaximum(nFrames);
            firstFrameSlider.setValue(1);
            lastFrameSlider.setMaximum(nFrames);
            lastFrameSlider.setValue(nFrames);
            resetFilterSliderRange();
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /**
     * Point the threshold slider at the range of the currently selected filter quantity. The
     * threshold starts at the low end so that nothing is hidden until the user asks for it.
     */
    private void resetFilterSliderRange() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            int filterType = filterCombo.getSelectedIndex();
            int lo = (int) Math.floor(filterMin[filterType]);
            int hi = (int) Math.ceil(filterMax[filterType]);

            // Widen the maximum first, otherwise setMinimum() can drag the old maximum along with it.
            minValueSlider.setMaximum(Math.max(hi, minValueSlider.getMaximum()));
            minValueSlider.setMinimum(lo);
            minValueSlider.setMaximum(hi);
            minValueSlider.setValue(lo);
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /**
     * Which histogram type is currently selected.
     */
    private int selectedType() {
        for (int i = 0; i < typeButtons.length; i++) {
            if (typeButtons[i].isSelected()) {
                return i;
            }
        }
        return TYPE_FRET;
    }

    /**
     * File name without its extension, used to suggest save names.
     */
    private static String stripExtension(File file) {
        String name = file.toString();
        int dotIndex = name.lastIndexOf('.');
        if (dotIndex > 0) {
            name = name.substring(0, dotIndex);
        }
        return name;
    }

    /**
     * Recompute the histogram and redraw. Called whenever a control changes.
     */
    private void update() {
        if (suspendUpdates) {
            return;
        }

        // The frame sliders are independent, so keep them from crossing over.
        if (firstFrameSlider.getValue() > lastFrameSlider.getValue()) {
            if (firstFrameSlider.getValueIsAdjusting()) {
                lastFrameSlider.setValue(firstFrameSlider.getValue());
            } else {
                firstFrameSlider.setValue(lastFrameSlider.getValue());
            }
        }

        result = computeHistogram(selectedType(),
                firstFrameSlider.getValue(),
                lastFrameSlider.getValue(),
                filterCombo.getSelectedIndex(),
                minValueSlider.getValue(),
                binsSlider.getValue());

        String status = String.format("%,d of %,d traces · frames %d-%d",
                result.nSpotsUsed, nSpots,
                firstFrameSlider.getValue(), lastFrameSlider.getValue());
        if (result.nOutside > 0) {
            status += String.format(" · %,d outside range", result.nOutside);
        }
        statusLabel.setText(status);

        plotPanel.repaint();
    }

    /**
     * Build the window.
     */
    private void showWindow() {
        JFrame frame = new JFrame("smFRET Trace Histograms - " + h5File.getName());
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        // Histogram type.
        typeButtons = new JRadioButton[TYPE_NAMES.length];
        ButtonGroup typeGroup = new ButtonGroup();
        JPanel typePanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 12, 2));
        typePanel.add(new JLabel("Histogram:"));
        for (int i = 0; i < TYPE_NAMES.length; i++) {
            typeButtons[i] = new JRadioButton(TYPE_NAMES[i], i == TYPE_FRET);
            typeButtons[i].addActionListener(e -> update());
            typeGroup.add(typeButtons[i]);
            typePanel.add(typeButtons[i]);
        }

        // File row.
        JPanel filePanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));
        JButton browseButton = new JButton("Browse...");
        browseButton.addActionListener(e -> onBrowse(frame));
        filePanel.add(new JLabel("H5 file:"));
        JLabel fileLabel = new JLabel(h5File.getName());
        filePanel.add(fileLabel);
        filePanel.add(browseButton);

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
        topPanel.add(filePanel);
        topPanel.add(typePanel);

        // Plot.
        plotPanel = new HistogramPanel();
        plotPanel.setBorder(new EmptyBorder(6, 6, 6, 6));

        // Controls. The sliders are created before resetSliderRanges() fills in their limits.
        // One point per trace means far fewer points than the old per frame histogram, so the
        // default bin count is correspondingly lower.
        binsSlider = new JSlider(10, 200, 30);
        firstFrameSlider = new JSlider(1, Math.max(1, nFrames), 1);
        lastFrameSlider = new JSlider(1, Math.max(1, nFrames), Math.max(1, nFrames));
        minValueSlider = new JSlider(0, 1, 0);

        // The threshold applies to one intensity at a time, chosen here. Switching rescales the
        // slider to the new quantity and clears the threshold, since the ranges are unrelated.
        filterCombo = new JComboBox<>(FILTER_NAMES);
        filterCombo.setSelectedIndex(FILTER_TOTAL);
        filterCombo.addActionListener(e -> {
            resetFilterSliderRange();
            update();
        });

        JPanel filterLabelPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        filterLabelPanel.add(new JLabel("Min"));
        filterLabelPanel.add(filterCombo);

        JPanel controlPanel = new JPanel(new GridBagLayout());
        controlPanel.setBorder(new EmptyBorder(4, 8, 4, 8));
        addSliderRow(controlPanel, 0, new JLabel("Bins"), binsSlider, () -> Integer.toString(binsSlider.getValue()));
        addSliderRow(controlPanel, 1, new JLabel("First frame"), firstFrameSlider, () -> Integer.toString(firstFrameSlider.getValue()));
        addSliderRow(controlPanel, 2, new JLabel("Last frame"), lastFrameSlider, () -> Integer.toString(lastFrameSlider.getValue()));
        addSliderRow(controlPanel, 3, filterLabelPanel, minValueSlider, () -> Integer.toString(minValueSlider.getValue()));

        // Status and save buttons.
        statusLabel = new JLabel(" ");
        JButton saveCsvButton = new JButton("Save CSV...");
        saveCsvButton.addActionListener(e -> onSaveCsv(frame));
        JButton savePngButton = new JButton("Save PNG...");
        savePngButton.addActionListener(e -> onSavePng(frame));

        JPanel bottomPanel = new JPanel(new BorderLayout(8, 0));
        bottomPanel.setBorder(new EmptyBorder(2, 10, 8, 10));
        bottomPanel.add(statusLabel, BorderLayout.CENTER);
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 0));
        buttonPanel.add(saveCsvButton);
        buttonPanel.add(savePngButton);
        bottomPanel.add(buttonPanel, BorderLayout.EAST);

        JPanel southPanel = new JPanel(new BorderLayout());
        southPanel.add(controlPanel, BorderLayout.CENTER);
        southPanel.add(bottomPanel, BorderLayout.SOUTH);

        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(topPanel, BorderLayout.NORTH);
        frame.getContentPane().add(plotPanel, BorderLayout.CENTER);
        frame.getContentPane().add(southPanel, BorderLayout.SOUTH);

        resetSliderRanges();
        update();

        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    /**
     * showWindow() helper, adds one labelled slider row with a live value readout.
     */
    private void addSliderRow(JPanel parent, int row, JComponent label, JSlider slider,
                              java.util.function.Supplier<String> valueText) {
        GridBagConstraints c = new GridBagConstraints();
        c.gridy = row;
        c.insets = new Insets(1, 2, 1, 6);

        c.gridx = 0;
        c.anchor = GridBagConstraints.WEST;
        parent.add(label, c);

        c.gridx = 1;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        parent.add(slider, c);

        c.gridx = 2;
        c.fill = GridBagConstraints.NONE;
        c.weightx = 0.0;
        JLabel valueLabel = new JLabel(valueText.get());
        valueLabel.setPreferredSize(new Dimension(58, valueLabel.getPreferredSize().height));
        parent.add(valueLabel, c);

        slider.addChangeListener(e -> {
            valueLabel.setText(valueText.get());
            update();
        });
    }

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
            if (isHeadless) {
                log.info("smFRET Trace Histograms is interactive and cannot run headless");
                return;
            }

            log.info("loading traces from " + h5File);
            loadTraces(h5File);

            SwingUtilities.invokeLater(this::showWindow);

        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
