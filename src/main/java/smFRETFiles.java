/*
 * Type checks for the files the user chooses in the parameter dialogs.
 *
 * The checks live here rather than in the widget because the widget cannot do them. SciJava's
 * file style does support "extensions:tif/tiff", but SwingFileWidget applies that filter only in
 * its drag and drop handler - the Browse button calls UIService.chooseFile with the style
 * collapsed back to a bare "open", and there is no single file chooseFile overload that takes a
 * filter. A path typed into the text field or passed by a macro bypasses the widget entirely in
 * any case, and the batch macros are how this pipeline is normally driven.
 *
 * The point of all of it is the message. Two of the four file kinds here are JSON written by
 * this pipeline, and from the outside they are indistinguishable - both sit beside the movie with
 * a '.json' extension. Handing the spot finder's own output back to it as a mapping used to parse
 * cleanly and then fail several steps later as "Cannot transform image, transform model not set",
 * which names neither the file nor the mistake. So every check here reports what the file
 * actually is as well as what was wanted.
 */

import com.fasterxml.jackson.databind.ObjectMapper;

import ij.IJ;
import ij.ImagePlus;

import org.scijava.log.LogService;

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

final class smFRETFiles {

    private smFRETFiles() {
    }

    // A key each of the two JSON files has and the other does not, which is the only thing that
    // tells them apart. The two key sets are disjoint - checked against the example data - so
    // either key identifies its file on its own.
    private static final String MAPPING_KEY = "source points";
    private static final String FINDING_KEY = "spots file";

    // What each kind is called, for the "is X" half of a message.
    private static final String IMAGE_IS = "an image";
    private static final String TIFF_IS = "a TIFF image";
    private static final String HDF5_IS = "an HDF5 trace file";
    private static final String MAPPING_IS = "a channel mapping JSON file";
    private static final String FINDING_IS = "a spot finder JSON file";
    private static final String TABLE_IS = "a spot table CSV";

    // Which stage writes the file that was wanted, for the sentence that says how to fix it.
    static final String MOVIE_REMEDY = " Choose the two channel movie, a TIFF stack.";
    static final String MAPPING_REMEDY =
            " Choose the '_mapping.json' written by smFRET Channel Mapper.";
    static final String FINDING_REMEDY =
            " Choose the '_spotf_finding.json' written by smFRET Spot Finder.";
    static final String HDF5_REMEDY = " Choose the '.h5' written by smFRET Time Traces.";

    // The two files named by the spot finder JSON rather than chosen, so the fix is not a choice.
    static final String TABLE_REMEDY =
            " The spot table is the '_spotf_spots.csv' written by smFRET Spot Finder.";
    static final String MASKS_REMEDY =
            " The masks are the '_spotf_masks.tif' written by smFRET Spot Finder.";

    // Enough to carry every signature tested for, and a header line for the CSV.
    private static final int HEAD_BYTES = 512;

    /**
     * Open a movie parameter, or say what the file is instead of an image.
     *
     * The sniff comes first so that the common mistakes - a JSON or a CSV where the movie was
     * wanted - never reach ImageJ's opener, whose own complaint is a dialog reading "Unsupported
     * format or file not found" and says nothing about which of this plugin's inputs was wrong.
     */
    static ImagePlus openImage(String path, String role) {
        return openImage((path == null) ? null : new File(path), role);
    }

    static ImagePlus openImage(File file, String role) {
        requireReadable(file, IMAGE_IS);

        String actual = describe(file);
        if ((actual != null) && !TIFF_IS.equals(actual)) {
            throw new smFRETAnalysisException("Error: " + role + " " + file + " is " + actual
                    + ", not " + IMAGE_IS + "." + MOVIE_REMEDY);
        }

        ImagePlus image = read(file);
        if (image.getProcessor() == null) {
            throw new smFRETAnalysisException("Error: could not read " + role + " " + file
                    + " - ImageJ does not recognize the format." + MOVIE_REMEDY);
        }

        smFRETChannelMapper.requireGrayscale(image, role + " " + file);
        smFRETChannelMapper.requireTimeStack(image, role + " " + file);
        return image;
    }

    /**
     * Read a file through ImageJ with its error dialog sent to the Log window instead.
     *
     * A format ImageJ does not know raises IJ.error inside the opener, which would put a dialog up
     * before this class gets to say anything more useful. redirectErrorMessages is a one shot
     * flag cleared by the error it redirects, so it has to be cleared again for the case where no
     * error fired and the flag would otherwise be left set for something unrelated.
     */
    static ImagePlus read(File file) {
        IJ.redirectErrorMessages();
        try {
            return new ImagePlus(file.toString());
        }
        finally {
            IJ.redirectErrorMessages(false);
        }
    }

    /**
     * The mapping written by smFRET Channel Mapper, as a map, or a message saying what this is.
     */
    static Map<String, Object> readMappingJSON(String path) {
        return readJSON((path == null) ? null : new File(path),
                MAPPING_KEY, MAPPING_IS, MAPPING_REMEDY);
    }

    /**
     * The parameters written by smFRET Spot Finder, as a map, or a message saying what this is.
     */
    static Map<String, Object> readSpotFinderJSON(File file) {
        return readJSON(file, FINDING_KEY, FINDING_IS, FINDING_REMEDY);
    }

    private static Map<String, Object> readJSON(File file, String key, String wanted, String remedy) {
        requireReadable(file, wanted);

        Map<String, Object> json = looksJSON(head(file)) ? parseJSON(file) : null;
        if ((json != null) && json.containsKey(key)) {
            return json;
        }

        String problem;
        if (json == null) {
            problem = isNot(file, wanted);
        }
        else {
            String kind = jsonKind(json);
            problem = (kind == null)
                    ? "is a JSON file but not one this pipeline wrote - it has no \"" + key + "\" entry"
                    : "is " + kind + ", not " + wanted;
        }
        throw new smFRETAnalysisException("Error: " + file + " " + problem + "." + remedy);
    }

    /**
     * Refuse anything that is not HDF5 before the HDF5 library is asked to open it.
     *
     * The signature is only looked for at the start of the file. HDF5 allows it at 512, 1024 and
     * so on, for a file with something else prepended, but nothing writes the traces that way.
     */
    static void requireHDF5(File file) {
        requireReadable(file, HDF5_IS);
        if (looksHDF5(head(file))) {
            return;
        }
        throw new smFRETAnalysisException("Error: " + file + " " + isNot(file, HDF5_IS)
                + "." + HDF5_REMEDY);
    }

    /**
     * Refuse anything that is not the spot table before ResultsTable is asked to parse it.
     *
     * ResultsTable.open2 does not fail on a binary file - it returns a table of whatever it made
     * of the bytes - so without this the caller's check for a missing column fires and blames a
     * file written by an older spot finder, which is a different problem with a different fix.
     */
    static void requireSpotTable(File file) {
        requireReadable(file, TABLE_IS);
        if (looksSpotTable(head(file))) {
            return;
        }
        throw new smFRETAnalysisException("Error: " + file + " " + isNot(file, TABLE_IS)
                + "." + TABLE_REMEDY);
    }

    /**
     * The "is a TIFF image, not an HDF5 trace file" clause, for a message that needs to say both.
     */
    static String isNot(File file, String wanted) {
        String actual = describe(file);
        return (actual == null) ? "is not " + wanted : "is " + actual + ", not " + wanted;
    }

    /**
     * The checks that come before knowing anything about the contents.
     *
     * Separate messages for each, because they are separate mistakes: a file that is not there is
     * usually a moved output, one that cannot be read is a permissions problem, and an empty one
     * is a run that died partway through writing it.
     */
    static void requireReadable(File file, String wanted) {
        if (file == null) {
            throw new smFRETAnalysisException("Error: no file was given where " + wanted
                    + " is needed.");
        }
        if (!file.exists()) {
            throw new smFRETAnalysisException("Error: " + file + " does not exist.");
        }
        if (file.isDirectory()) {
            throw new smFRETAnalysisException("Error: " + file + " is a directory, not " + wanted + ".");
        }
        if (!file.canRead()) {
            throw new smFRETAnalysisException("Error: " + file
                    + " cannot be read - check its permissions.");
        }
        if (file.length() == 0) {
            throw new smFRETAnalysisException("Error: " + file + " is empty.");
        }
    }

    /**
     * Show one of this plugin's own validation messages, without a stack trace.
     *
     * smFRETAnalysisException is raised only by this plugin's own checks, so it always carries a
     * message written for the person who chose the file. Everything else reaching a run() catch is
     * unexpected and its stack trace is the useful part, which is why the two are separated at the
     * catch rather than here.
     *
     * A dialog only when there is somebody to dismiss it. IJ.error under a macro aborts the macro,
     * which is right, but it puts up a modal dialog first, and the batch macros run unattended.
     */
    static void report(LogService log, smFRETAnalysisException e) {
        String message = (e.getMessage() == null) ? e.toString() : e.getMessage();
        log.info(e);
        if (GraphicsEnvironment.isHeadless() || IJ.isMacro()) {
            log.error(message);
        }
        else {
            IJ.error("smFRET", message);
        }
    }

    /**
     * What the file actually is, as a phrase, or null if it is nothing this pipeline recognizes.
     */
    private static String describe(File file) {
        byte[] head = head(file);
        if (head.length == 0) {
            return null;
        }
        if (looksTIFF(head)) {
            return TIFF_IS;
        }
        if (looksHDF5(head)) {
            return HDF5_IS;
        }
        if (looksJSON(head)) {
            Map<String, Object> json = parseJSON(file);
            if (json != null) {
                String kind = jsonKind(json);
                return (kind == null) ? "a JSON file, but not one this pipeline wrote" : kind;
            }
        }
        if (looksSpotTable(head)) {
            return TABLE_IS;
        }
        return null;
    }

    private static String jsonKind(Map<String, Object> json) {
        if (json.containsKey(MAPPING_KEY)) {
            return MAPPING_IS;
        }
        if (json.containsKey(FINDING_KEY)) {
            return FINDING_IS;
        }
        return null;
    }

    /**
     * Byte order mark plus magic number, for classic TIFF and for BigTIFF.
     */
    private static boolean looksTIFF(byte[] head) {
        if (head.length < 4) {
            return false;
        }
        boolean little = (head[0] == 'I') && (head[1] == 'I');
        boolean big = (head[0] == 'M') && (head[1] == 'M');
        int magic = little ? ((head[3] & 0xff) << 8 | (head[2] & 0xff))
                : ((head[2] & 0xff) << 8 | (head[3] & 0xff));
        return (little || big) && ((magic == 42) || (magic == 43));
    }

    private static final byte[] HDF5_SIGNATURE =
            {(byte) 0x89, 'H', 'D', 'F', '\r', '\n', 0x1a, '\n'};

    private static boolean looksHDF5(byte[] head) {
        if (head.length < HDF5_SIGNATURE.length) {
            return false;
        }
        for (int i = 0; i < HDF5_SIGNATURE.length; i++) {
            if (head[i] != HDF5_SIGNATURE[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Whether it is worth handing the file to Jackson at all.
     *
     * Both JSON files here are objects, so the first thing that is not whitespace is '{'. Testing
     * that first keeps a movie from being streamed through a JSON parser to learn what its first
     * two bytes already said.
     */
    private static boolean looksJSON(byte[] head) {
        for (byte b : head) {
            if ((b == ' ') || (b == '\t') || (b == '\r') || (b == '\n')) {
                continue;
            }
            return b == '{';
        }
        return false;
    }

    /**
     * The spot table's own header line - x and y first, whether or not ResultsTable wrote a row
     * number column ahead of them.
     */
    private static boolean looksSpotTable(byte[] head) {
        int end = 0;
        while ((end < head.length) && (head[end] != '\n') && (head[end] != '\r')) {
            end++;
        }
        String[] columns = new String(head, 0, end, java.nio.charset.StandardCharsets.ISO_8859_1)
                .split(",", -1);
        int first = ((columns.length > 0) && columns[0].trim().isEmpty()) ? 1 : 0;
        return (columns.length >= first + 2)
                && "x".equals(columns[first].trim())
                && "y".equals(columns[first + 1].trim());
    }

    private static Map<String, Object> parseJSON(File file) {
        try {
            @SuppressWarnings("unchecked")
            Map<String, Object> json = new ObjectMapper().readValue(file, HashMap.class);
            return json;
        }
        catch (Exception e) {
            return null;
        }
    }

    private static byte[] head(File file) {
        byte[] buffer = new byte[HEAD_BYTES];
        try (InputStream in = new FileInputStream(file)) {
            int read = 0;
            while (read < HEAD_BYTES) {
                int n = in.read(buffer, read, HEAD_BYTES - read);
                if (n < 0) {
                    break;
                }
                read += n;
            }
            return Arrays.copyOf(buffer, read);
        }
        catch (IOException e) {
            return new byte[0];
        }
    }
}
