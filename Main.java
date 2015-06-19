package nhs.genetics.cardiff;

import java.io.File;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 3){
            System.err.println("Usage: <File> <List> <Min>");
            System.exit(1);
        }

        //get file extension
        String[] fields = args[0].split("\\.");

        //read bed list
        ListReader bedFilePathList = new ListReader(new File(args[1]));
        bedFilePathList.parseListReader();

        if (fields[1].equals("vcf")){ //haplotypecaller

            //calculate coverage from GATK HaplotypeCaller
            Coverage coverage = new Coverage(new File(args[0]), bedFilePathList, Integer.parseInt(args[2]));
            coverage.populateTargetBases();
            coverage.extractBasesPassingMinThresholdFromVCF();
            coverage.extractBasesFailingMinThreshold();
            coverage.windowMissingBases();
            coverage.printCoverageMetrics();

        } else { //depth of coverage

            //calculate coverage from GATK depth of coverage
            Coverage coverage = new Coverage(new File(args[0]), bedFilePathList, Integer.parseInt(args[2]));
            coverage.populateTargetBases();
            coverage.extractBasesPassingMinThresholdFromDepthOfCoverage();
            coverage.extractBasesFailingMinThreshold();
            coverage.windowMissingBases();
            coverage.printCoverageMetrics();

        }

    }
}
