import org.scijava.Context;
import org.scijava.log.LogService;
import java.io.File;

/**
 * Stage 2 on a simulated movie with the tolerance and prominence forced, for the sweeps.
 *
 * Pass a prominence below zero to disable that filter, which is what the prominence sweeps do:
 * it is the last stage and affects nothing upstream, so one run with it off gives the score of
 * every surviving spot and any threshold can then be applied offline.
 *
 * args: &lt;dir&gt; &lt;sigma&gt; &lt;tolerance&gt; &lt;prominence&gt; &lt;frames&gt; [threshold]
 */
public class RunTune {
    public static void main(String[] args) throws Exception {
        String dir = args[0];

        Context context = new Context(LogService.class);
        LogService log = context.getService(LogService.class);

        smFRETSpotFinder s = new smFRETSpotFinder();
        s.log = log;
        s.inputImageName = new File(dir, "sim.tif");
        s.mappingFile = new File(dir, "sim_mapping.json");
        s.startSlice = 1;
        s.endSlice = Integer.parseInt(args[4]);
        s.spotSigma = Double.parseDouble(args[1]);
        s.spotTolerance = Double.parseDouble(args[2]);
        s.spotProminence = Double.parseDouble(args[3]);
        s.spotThreshold = (args.length > 5) ? Double.parseDouble(args[5]) : 6.0;
        s.cameraBlackLevel = 5; s.cameraGain = 1.0;
        s.spotSpacing = 3; s.edgeMargin = 5; s.backgroundKappa = 0.0;
        s.run();
        System.exit(0);
    }
}
