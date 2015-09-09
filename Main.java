package nhs.genetics.cardiff;

import java.io.File;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 3){
            System.err.println("Usage: <GATKCoverageFile> <BED> <Min>");
            System.exit(1);
        }

        //calculate coverage from GATK depth of coverage
        Coverage coverage = new Coverage(new File(args[0]), new File(args[1]), Integer.parseInt(args[2]));
        coverage.populateTargetBases();
        coverage.extractBasesPassingMinDepth();
        coverage.extractBasesFailingMinDepth();
        coverage.windowMissingBases();
        coverage.printCoverageMetrics();

    }
}
