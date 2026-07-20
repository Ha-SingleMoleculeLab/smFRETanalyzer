/*
 * This class finds the single molecule spots using the mapping.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.ImageCalculator;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;


@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>smFRET>smFRET Spot Finder")
public class smFRETSpotFinder implements Command {
    @Parameter
    LogService log;

    @Parameter
    UIService ui;

    @Parameter
    ImagePlus img;

    @Parameter (description = "first slice for averaging", min = "1")
    Integer startSlice = 1;

    @Parameter (description = "last slice for averaging", min = "1")
    Integer endSlice = 30;

    @Parameter(description = "Channel to channel mapping file", label = "Mapping file", style = "open")
    File mappingFile;

    @Parameter(description = "Directory to save results in", label = "Save Directory", style = "directory")
    File saveDirectory;

    // Member variables.
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
	        log.info("starting spot finding");
            smFRETChannelMapper smfcm = new smFRETChannelMapper();
            smfcm.log = log;

            // Average image.
            log.info("average image - " + img.getNSlices() + " slices");
            ImagePlus averageImage;
            if (img.getNSlices() == 1){
                averageImage = img.duplicate();
            }
            else {
                averageImage = smfcm.averageImagePlus(img, startSlice, endSlice);
            }

            // split, transform and add the two channels together.
            log.info("split and transform");
            smfcm.loadMappingJSON(mappingFile);
            java.util.List<ImagePlus> images = smfcm.splitImagePlus(averageImage, true);

            ImagePlus sumImage = ImageCalculator.run(images.get(0), images.get(1), "add create");
            if (!this.isHeadless) {
                ui.show(sumImage);
            }

            // find peaks in sum image.

	        log.info("finishing spot finding");
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
