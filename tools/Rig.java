import org.scijava.Context;
import org.scijava.log.LogService;
import java.io.File;

/** Rebuild stage 2 and stage 3 outputs for one movie. args: <dir> <image.tif> <mapping.json> */
public class Rig {
    public static void main(String[] args) throws Exception {
        Context context = new Context(LogService.class);
        LogService log = context.getService(LogService.class);

        smFRETSpotFinder s2 = new smFRETSpotFinder();
        s2.log = log;
        s2.inputImageName = new File(args[0], args[1]);
        s2.mappingFile = new File(args[2]);
        s2.startSlice = 1; s2.endSlice = 30;
        s2.spotThreshold = 6.0; s2.spotTolerance = 5.0; s2.spotProminence = 0.4;
        s2.spotSigma = 2.0; s2.cameraBlackLevel = 5; s2.cameraGain = 1.0;
        s2.spotSpacing = 2; s2.edgeMargin = 5; s2.backgroundKappa = 0.0;
        s2.run();

        String root = new File(args[0], args[1]).toString();
        root = root.substring(0, root.lastIndexOf('.'));

        smFRETAnalyzer s3 = new smFRETAnalyzer();
        s3.log = log;
        s3.spotJSONFile = new File(root + "_spotf_finding.json");
        s3.backgroundAverageNFrames = 30;
        s3.run();

        System.out.println("### rig built for " + root);
        System.exit(0);
    }
}
