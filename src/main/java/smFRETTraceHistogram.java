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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;


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

    // Quantity the intensity range is applied to. These are mutually exclusive, so they are a
    // combo box beside a single slider rather than one slider each.
    private static final int FILTER_TOTAL = 0;
    private static final int FILTER_DONOR = 1;
    private static final int FILTER_ACCEPTOR = 2;
    private static final String[] FILTER_NAMES = {"Total (D+A)", "Donor (target)", "Acceptor (source)"};

    // FRET efficiency is plotted over a fixed range, slightly wider than [0,1] so that the
    // noise skirts either side of the physical range stay visible.
    private static final double FRET_MIN = -0.2;
    private static final double FRET_MAX = 1.2;

    // The traces as smFRETAnalyzer wrote them, no corrections applied.
    private static final Corrections NO_CORRECTIONS = new Corrections(0.0, 0.0, 0.0);

    // Member variables.
    private JSpinner acceptorBaselineSpinner;
    private JSlider binsSlider;
    private JSpinner donorBaselineSpinner;
    private JComboBox<String> filterCombo;
    private final double[] filterMax = new double[FILTER_NAMES.length];
    private final double[] filterMin = new double[FILTER_NAMES.length];
    private RangeSlider frameRangeSlider;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private JSpinner leakageSpinner;
    private RangeSlider valueRangeSlider;
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
        int nSpotsUsed;     // Traces inside the intensity range.
        String valueLabel;
    }

    /**
     * The per frame corrections applied to the traces before anything is measured from them.
     *
     * Both baselines are *subtracted*, so a negative value adds one. Leakage is the fraction of
     * the donor that appears in the acceptor channel, and it is taken off the acceptor using the
     * *baseline corrected* donor: leakage is a fraction of real donor emission, and an offset
     * that survived into the trace is not donor emission.
     *
     * Grouped rather than passed as three more doubles because computeHistogram() already takes
     * seven positional arguments, and three unlabelled doubles in a row are easy to transpose.
     */
    static class Corrections {

        final double acceptorBaseline;
        final double donorBaseline;
        final double leakage;

        Corrections(double donorBaseline, double acceptorBaseline, double leakage) {
            this.acceptorBaseline = acceptorBaseline;
            this.donorBaseline = donorBaseline;
            this.leakage = leakage;
        }

        /**
         * Takes the *corrected* donor rather than the raw one, both because that is the quantity
         * the leakage is a fraction of and so that callers do not correct the donor twice.
         */
        double correctAcceptor(double rawAcceptor, double correctedDonor) {
            return rawAcceptor - acceptorBaseline - leakage * correctedDonor;
        }

        double correctDonor(double rawDonor) {
            return rawDonor - donorBaseline;
        }

        String describe() {
            return "baseline D " + compact(donorBaseline) + " / A " + compact(acceptorBaseline)
                    + ", leakage " + compact(leakage);
        }

        boolean isIdentity() {
            return (acceptorBaseline == 0.0) && (donorBaseline == 0.0) && (leakage == 0.0);
        }
    }

    /**
     * Bin the loaded traces. Takes its settings as arguments rather than reading the controls
     * directly so that the binning can be exercised without a GUI.
     */
    Histogram computeHistogram(int type, int firstFrame, int lastFrame,
                               int filterType, double minValue, double maxValue, int nBins,
                               Corrections corrections) {

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
            double highestFrameValue = -Double.MAX_VALUE;
            for (int t = firstFrame - 1; t < lastFrame; t++) {
                double frameDonor = corrections.correctDonor(targetTraces[i][t]);
                double frameAcceptor = corrections.correctAcceptor(sourceTraces[i][t], frameDonor);
                donorSum += frameDonor;
                acceptorSum += frameAcceptor;

                double frameValue = filterValue(filterType, frameDonor, frameAcceptor,
                        frameDonor + frameAcceptor);
                if (frameValue < lowestFrameValue) {
                    lowestFrameValue = frameValue;
                }
                if (frameValue > highestFrameValue) {
                    highestFrameValue = frameValue;
                }
            }

            // The whole trace goes if any single frame in the interval falls outside the range,
            // so a molecule that bleaches part way through contributes nothing rather than a
            // diluted average. The maximum works the same way by design, which makes it strict:
            // one bright frame is enough to drop a trace. That is what catches an aggregate, and
            // it also means a single spike will do it.
            if ((lowestFrameValue < minValue) || (highestFrameValue > maxValue)) {
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
     * The intensity the range slider is currently applied to.
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
     * Two handle slider for choosing the frame interval to average over.
     *
     * Swing has no range slider. One ordinary slider per end can express the same interval, but
     * makes the common operation - moving a window of fixed width through the movie - awkward,
     * because both ends have to be dragged separately and kept in step. Here dragging either
     * handle resizes the interval and dragging the bar between them slides it at constant width.
     *
     * Only low <= high states are reachable, so callers do not have to order the two values.
     */
    private static class RangeSlider extends JComponent {

        private static final int BAR_HEIGHT = 6;
        private static final int THUMB_SIZE = 13;

        private static final int DRAG_NONE = 0;
        private static final int DRAG_LOW = 1;
        private static final int DRAG_HIGH = 2;
        private static final int DRAG_BAR = 3;

        private int dragMode = DRAG_NONE;
        private int grabHigh;               // Interval at the start of a bar drag, so that the
        private int grabLow;                // width is preserved exactly however far it is
        private int grabValue;              // dragged, rather than drifting by a rounding error.
        private int high;
        private final ArrayList<ChangeListener> listeners = new ArrayList<>();
        private int low;
        private int maximum;
        private int minimum;
        private boolean valueIsAdjusting;

        RangeSlider(int minimum, int maximum) {
            this.minimum = minimum;
            this.maximum = Math.max(maximum, minimum);
            low = this.minimum;
            high = this.maximum;

            setPreferredSize(new Dimension(200, THUMB_SIZE + 8));
            setFocusable(true);

            // The keys are worth more than the mouse on a long movie - a track a few hundred
            // pixels wide cannot address every frame of it - and there is no way to discover
            // them, nor that the slider has to be clicked first to take the focus, other than
            // being told here. HTML because a tooltip does not otherwise wrap.
            //
            // Worded in steps rather than frames: both the frame interval and the intensity
            // range are this same component.
            setToolTipText("<html>Drag either end to resize the interval, "
                    + "drag the middle to slide it.<br>"
                    + "Click first, then:<br>"
                    + "&nbsp;&nbsp;<b>\u2190</b> <b>\u2192</b> move the interval by one step<br>"
                    + "&nbsp;&nbsp;<b>PgUp</b> <b>PgDn</b> move it a whole interval<br>"
                    + "&nbsp;&nbsp;<b>Home</b> <b>End</b> jump to the start or the end<br>"
                    + "&nbsp;&nbsp;<b>\u2191</b> <b>\u2193</b> resize it</html>");

            MouseAdapter mouse = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    requestFocusInWindow();
                    onPress(e.getX());
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    onDrag(e.getX());
                }

                @Override
                public void mouseReleased(MouseEvent e) {
                    dragMode = DRAG_NONE;

                    // Fire once more on release so listeners that only act on settled values
                    // see the final interval.
                    if (valueIsAdjusting) {
                        valueIsAdjusting = false;
                        fireChange();
                    }
                }
            };
            addMouseListener(mouse);
            addMouseMotionListener(mouse);

            addKeyListener(new KeyAdapter() {
                @Override
                public void keyPressed(KeyEvent e) {
                    onKey(e.getKeyCode());
                }
            });

            addFocusListener(new FocusAdapter() {
                @Override
                public void focusGained(FocusEvent e) {
                    repaint();
                }

                @Override
                public void focusLost(FocusEvent e) {
                    repaint();
                }
            });
        }

        void addChangeListener(ChangeListener listener) {
            listeners.add(listener);
        }

        private int clamp(int value) {
            if (value < minimum) {
                return minimum;
            }
            if (value > maximum) {
                return maximum;
            }
            return value;
        }

        private void fireChange() {
            ChangeEvent event = new ChangeEvent(this);
            for (ChangeListener listener : listeners) {
                listener.stateChanged(event);
            }
        }

        int getHigh() {
            return high;
        }

        int getLow() {
            return low;
        }

        private void onDrag(int x) {
            if (dragMode == DRAG_NONE) {
                return;
            }
            valueIsAdjusting = true;
            int value = xToValue(x);

            if (dragMode == DRAG_LOW) {
                setValues(Math.min(value, high), high);
            } else if (dragMode == DRAG_HIGH) {
                setValues(low, Math.max(value, low));
            } else {
                slide(grabLow + (value - grabValue), grabHigh - grabLow);
            }
        }

        private void onKey(int keyCode) {
            int width = high - low;
            if (keyCode == KeyEvent.VK_LEFT) {
                slide(low - 1, width);
            } else if (keyCode == KeyEvent.VK_RIGHT) {
                slide(low + 1, width);
            } else if (keyCode == KeyEvent.VK_PAGE_UP) {

                // A whole interval per press rather than a fixed number of frames, so each press
                // lands on the next window that does not overlap this one. A fixed step would be
                // either useless on a 1295 frame movie or far too coarse on a 30 frame one; this
                // scales with whatever is being looked at.
                //
                // width + 1, not width: the interval is inclusive, so 10..14 is five frames and
                // stepping by four would leave frame 14 in both windows.
                slide(low - (width + 1), width);
            } else if (keyCode == KeyEvent.VK_PAGE_DOWN) {
                slide(low + (width + 1), width);
            } else if (keyCode == KeyEvent.VK_HOME) {
                slide(minimum, width);
            } else if (keyCode == KeyEvent.VK_END) {
                slide(maximum - width, width);
            } else if (keyCode == KeyEvent.VK_UP) {
                setValues(low, high + 1);
            } else if (keyCode == KeyEvent.VK_DOWN) {
                setValues(low, Math.max(high - 1, low));
            }
        }

        private void onPress(int x) {
            int lowX = valueToX(low);
            int highX = valueToX(high);
            int lowDistance = Math.abs(x - lowX);
            int highDistance = Math.abs(x - highX);

            if ((lowDistance <= THUMB_SIZE) || (highDistance <= THUMB_SIZE)) {

                // Which handle was grabbed. When the interval is empty the two coincide, and the
                // side of the handle that was clicked decides - otherwise a collapsed interval
                // could only ever be reopened in one direction.
                if ((lowDistance < highDistance) || ((lowDistance == highDistance) && (x <= lowX))) {
                    dragMode = DRAG_LOW;
                } else {
                    dragMode = DRAG_HIGH;
                }
            } else if ((x > lowX) && (x < highX)) {
                dragMode = DRAG_BAR;
                grabValue = xToValue(x);
                grabLow = low;
                grabHigh = high;
            } else {

                // On the track outside the interval, bring the nearer end out to meet the click.
                dragMode = (xToValue(x) < low) ? DRAG_LOW : DRAG_HIGH;
            }
            onDrag(x);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int barY = (getHeight() - BAR_HEIGHT) / 2;
            int lowX = valueToX(low);
            int highX = valueToX(high);

            g2.setColor(new Color(205, 205, 205));
            g2.fillRoundRect(trackX(), barY, trackWidth(), BAR_HEIGHT, BAR_HEIGHT, BAR_HEIGHT);

            // The selected interval, in the same blue as the histogram bars.
            g2.setColor(new Color(70, 115, 175));
            g2.fillRoundRect(lowX, barY, Math.max(1, highX - lowX), BAR_HEIGHT, BAR_HEIGHT, BAR_HEIGHT);

            paintThumb(g2, lowX);
            paintThumb(g2, highX);

            g2.dispose();
        }

        private void paintThumb(Graphics2D g2, int x) {
            int y = (getHeight() - THUMB_SIZE) / 2;
            g2.setColor(Color.WHITE);
            g2.fillOval(x - THUMB_SIZE / 2, y, THUMB_SIZE, THUMB_SIZE);
            g2.setColor(isFocusOwner() ? new Color(40, 80, 140) : new Color(90, 90, 90));
            g2.drawOval(x - THUMB_SIZE / 2, y, THUMB_SIZE - 1, THUMB_SIZE - 1);
        }

        /**
         * Set the limits, keeping the interval inside them.
         */
        void setRange(int newMinimum, int newMaximum) {
            minimum = newMinimum;
            maximum = Math.max(newMaximum, newMinimum);
            setValues(low, high);
            repaint();
        }

        /**
         * Set the interval, clamped to the limits and to low <= high.
         */
        void setValues(int newLow, int newHigh) {
            int clampedLow = clamp(newLow);
            int clampedHigh = Math.max(clamp(newHigh), clampedLow);
            if ((clampedLow == low) && (clampedHigh == high)) {
                return;
            }

            low = clampedLow;
            high = clampedHigh;
            repaint();
            fireChange();
        }

        /**
         * Move the interval to start at newLow without changing its width, stopping at the ends
         * rather than letting the width shrink against them.
         */
        private void slide(int newLow, int width) {
            int clampedLow = Math.min(Math.max(newLow, minimum), maximum - width);
            setValues(clampedLow, clampedLow + width);
        }

        // Half a handle of padding at each end, so that the handles stay inside the component
        // when the interval is at either limit.
        private int trackWidth() {
            return Math.max(1, getWidth() - THUMB_SIZE);
        }

        private int trackX() {
            return THUMB_SIZE / 2;
        }

        private int valueToX(int value) {
            if (maximum == minimum) {
                return trackX();
            }
            return trackX() + (int) Math.round((double) (value - minimum) * trackWidth() / (maximum - minimum));
        }

        private int xToValue(int x) {
            if (maximum == minimum) {
                return minimum;
            }
            double fraction = (double) (x - trackX()) / trackWidth();
            return clamp(minimum + (int) Math.round(fraction * (maximum - minimum)));
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
     * A correction value for the status line and the CSV header, without the trailing zeroes.
     *
     * Not DecimalFormat, which would put a comma where a locale wants a decimal separator and
     * make the CSV header lie. %f always emits a decimal point, so stripping trailing zeroes
     * cannot eat a digit of an integer the way it would on a bare "100".
     */
    private static String compact(double value) {
        String text = String.format(java.util.Locale.US, "%.6f", value);
        while (text.endsWith("0")) {
            text = text.substring(0, text.length() - 1);
        }
        return text.endsWith(".") ? text.substring(0, text.length() - 1) : text;
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
        computeFilterBounds(corrections());

        log.info("loaded " + nSpots + " traces of " + nFrames + " frames from " + file);
    }

    /**
     * Range of each quantity the intensity range slider can be applied to, which sets its limits.
     *
     * Measured on the corrected traces, so it has to be redone whenever a correction changes -
     * subtracting a baseline moves the data out from under bounds taken on the raw traces, and an
     * end of the slider parked at a stale extreme would then silently exclude traces.
     */
    private void computeFilterBounds(Corrections corrections) {
        for (int f = 0; f < FILTER_NAMES.length; f++) {
            filterMin[f] = Double.MAX_VALUE;
            filterMax[f] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < nSpots; i++) {
            for (int t = 0; t < nFrames; t++) {
                double donor = corrections.correctDonor(targetTraces[i][t]);
                double acceptor = corrections.correctAcceptor(sourceTraces[i][t], donor);
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
    }

    /**
     * The corrections the spinboxes currently ask for.
     *
     * loadTraces() runs before the window is built, and again from Browse... after it is, so this
     * has to answer in both states - uncorrected before there are controls to read.
     */
    private Corrections corrections() {
        if (donorBaselineSpinner == null) {
            return NO_CORRECTIONS;
        }
        return new Corrections(((Number) donorBaselineSpinner.getValue()).doubleValue(),
                ((Number) acceptorBaselineSpinner.getValue()).doubleValue(),
                ((Number) leakageSpinner.getValue()).doubleValue());
    }

    /**
     * A correction changed: the traces have moved, so the intensity range slider is rescaled to
     * them and reopened, exactly as switching the quantity it applies to already does.
     */
    private void onCorrectionChanged() {
        if (suspendUpdates) {
            return;
        }
        computeFilterBounds(corrections());
        resetFilterSliderRange();
        update();
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
            writer.println("# frames " + frameRangeSlider.getLow() + "-" + frameRangeSlider.getHigh()
                    + ", " + FILTER_NAMES[filterCombo.getSelectedIndex()] + " "
                    + valueRangeSlider.getLow() + "-" + valueRangeSlider.getHigh()
                    + ", " + result.nPoints + " of " + result.nSpotsUsed + " traces in range");

            // Written even when they are all zero, so that a saved histogram says what was done
            // to the traces rather than leaving it to be inferred from the absence of a line.
            writer.println("# corrections: " + corrections().describe());
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
     * Point the frame and intensity sliders at the currently loaded traces.
     */
    private void resetSliderRanges() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            frameRangeSlider.setRange(1, Math.max(1, nFrames));
            frameRangeSlider.setValues(1, Math.max(1, nFrames));
            resetFilterSliderRange();
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /**
     * Point the intensity range slider at the range of the currently selected filter quantity.
     * Both handles start at the extremes so that nothing is hidden until the user asks for it.
     *
     * The bounds are the per-frame minimum and maximum rather than the range of the trace
     * averages, which is what guarantees neither end can silently exclude a trace when it is
     * parked at its extreme.
     */
    private void resetFilterSliderRange() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            int filterType = filterCombo.getSelectedIndex();
            int lo = (int) Math.floor(filterMin[filterType]);
            int hi = (int) Math.ceil(filterMax[filterType]);

            valueRangeSlider.setRange(lo, hi);
            valueRangeSlider.setValues(lo, hi);
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

        Corrections corrections = corrections();
        result = computeHistogram(selectedType(),
                frameRangeSlider.getLow(),
                frameRangeSlider.getHigh(),
                filterCombo.getSelectedIndex(),
                valueRangeSlider.getLow(),
                valueRangeSlider.getHigh(),
                binsSlider.getValue(),
                corrections);

        String status = String.format("%,d of %,d traces · frames %d-%d (%,d wide)",
                result.nSpotsUsed, nSpots,
                frameRangeSlider.getLow(), frameRangeSlider.getHigh(),
                frameRangeSlider.getHigh() - frameRangeSlider.getLow() + 1);
        if (result.nOutside > 0) {
            status += String.format(" · %,d outside range", result.nOutside);
        }

        // Only when they are doing something, so that the common uncorrected case reads as it
        // did before rather than carrying three zeroes around.
        if (!corrections.isIdentity()) {
            status += " · " + corrections.describe();
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
        frameRangeSlider = new RangeSlider(1, Math.max(1, nFrames));
        valueRangeSlider = new RangeSlider(0, 1);

        // The range applies to one intensity at a time, chosen here. Switching rescales the
        // slider to the new quantity and reopens it to the full range, since the ranges of the
        // three quantities are unrelated.
        filterCombo = new JComboBox<>(FILTER_NAMES);
        filterCombo.setSelectedIndex(FILTER_TOTAL);
        filterCombo.addActionListener(e -> {
            resetFilterSliderRange();
            update();
        });

        JPanel filterLabelPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        filterLabelPanel.add(new JLabel("Range"));
        filterLabelPanel.add(filterCombo);

        // Corrections. Spin boxes rather than sliders because these are set to a measured number -
        // a baseline read off a blank region, a leakage measured on a donor only sample - rather
        // than dialled in by eye, and typing the number is the point.
        donorBaselineSpinner = correctionSpinner(new SpinnerNumberModel(0.0, -1.0e9, 1.0e9, 1.0), "0.###",
                "Subtracted from every donor value. Negative adds one.");
        acceptorBaselineSpinner = correctionSpinner(new SpinnerNumberModel(0.0, -1.0e9, 1.0e9, 1.0), "0.###",
                "Subtracted from every acceptor value. Negative adds one.");
        // Four decimals rather than the baselines' three: a leakage coefficient is a small
        // fraction, and the editor commits what it displays, so the format is the precision.
        leakageSpinner = correctionSpinner(new SpinnerNumberModel(0.0, 0.0, 1.0, 0.01), "0.####",
                "Fraction of the baseline corrected donor removed from the acceptor.");

        JPanel correctionPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 1));
        correctionPanel.add(new JLabel("Baseline   donor"));
        correctionPanel.add(donorBaselineSpinner);
        correctionPanel.add(new JLabel("acceptor"));
        correctionPanel.add(acceptorBaselineSpinner);
        correctionPanel.add(Box.createHorizontalStrut(14));
        correctionPanel.add(new JLabel("Donor leakage"));
        correctionPanel.add(leakageSpinner);

        JPanel controlPanel = new JPanel(new GridBagLayout());
        controlPanel.setBorder(new EmptyBorder(4, 8, 4, 8));
        addSliderRow(controlPanel, 0, new JLabel("Bins"), binsSlider, () -> Integer.toString(binsSlider.getValue()));
        addSliderRow(controlPanel, 1, new JLabel("Frames"), frameRangeSlider,
                () -> frameRangeSlider.getLow() + "-" + frameRangeSlider.getHigh());
        addSliderRow(controlPanel, 2, filterLabelPanel, valueRangeSlider,
                () -> valueRangeSlider.getLow() + "-" + valueRangeSlider.getHigh());

        GridBagConstraints correctionConstraints = new GridBagConstraints();
        correctionConstraints.anchor = GridBagConstraints.WEST;
        correctionConstraints.gridwidth = 3;
        correctionConstraints.gridx = 0;
        correctionConstraints.gridy = 3;
        correctionConstraints.insets = new Insets(1, 0, 1, 6);
        controlPanel.add(correctionPanel, correctionConstraints);

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
     * showWindow() helper, one correction spin box.
     */
    private JSpinner correctionSpinner(SpinnerNumberModel model, String format, String tip) {
        JSpinner spinner = new JSpinner(model);
        JSpinner.NumberEditor editor = new JSpinner.NumberEditor(spinner, format);
        editor.getTextField().setColumns(6);
        spinner.setEditor(editor);
        spinner.setToolTipText(tip);
        spinner.addChangeListener(e -> onCorrectionChanged());
        return spinner;
    }

    /**
     * showWindow() helper, adds one labelled slider row with a live value readout.
     *
     * JSlider and RangeSlider have the same addChangeListener() method but no common supertype
     * that declares it, hence the pair of overloads over a shared layout helper.
     */
    private void addSliderRow(JPanel parent, int row, JComponent label, JSlider slider,
                              java.util.function.Supplier<String> valueText) {
        JLabel valueLabel = addControlRow(parent, row, label, slider, valueText);
        slider.addChangeListener(e -> {
            valueLabel.setText(valueText.get());
            update();
        });
    }

    /**
     * showWindow() helper, the range slider flavour of addSliderRow().
     */
    private void addSliderRow(JPanel parent, int row, JComponent label, RangeSlider slider,
                              java.util.function.Supplier<String> valueText) {
        JLabel valueLabel = addControlRow(parent, row, label, slider, valueText);
        slider.addChangeListener(e -> {
            valueLabel.setText(valueText.get());
            update();
        });
    }

    /**
     * addSliderRow() helper, lays out one row and returns its value readout.
     */
    private JLabel addControlRow(JPanel parent, int row, JComponent label, JComponent control,
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
        parent.add(control, c);

        c.gridx = 2;
        c.fill = GridBagConstraints.NONE;
        c.weightx = 0.0;
        JLabel valueLabel = new JLabel(valueText.get());

        // Wide enough for the frame row's "1234-5678", the longest readout of the three.
        valueLabel.setPreferredSize(new Dimension(78, valueLabel.getPreferredSize().height));
        parent.add(valueLabel, c);

        return valueLabel;
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
