/*
 * This class finds the single molecule spots using the mapping.
 */

import ij.IJ;
import ij.ImagePlus;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;


@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Spot Finder")
public class smFRETSpotFinder implements Command {
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
     * Run ...
     */
    @Override
    public void run() {
        try {
	    log.info("starting spot finding");
	    log.info("finishing spot finding");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
