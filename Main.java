package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.TabixReader;
import nhs.genetics.cardiff.framework.GFF3Record;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Program for reporting coverage metrics with gene context
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-03-15
 */
public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final String program = "CoverageCalculator";
    private static final String version = "3.0.0";

    public static void main(String[] args) {

        //parse command line
        CommandLineParser commandLineParser = new DefaultParser();
        CommandLine commandLine = null;
        HelpFormatter formatter = new HelpFormatter();
        Options options = new Options();

        options.addOption("T", "Transcript", true, "Path to GTF/GFF file");
        options.addOption("C", "Coverage", true, "Path to per base coverage");
        options.addOption("R", "Regions", true, "Path to target BED file");
        options.addOption("P", "Padding", true, "Base-pair padding around GTF entries");
        options.addOption("D", "Depth", true, "Minimum depth to report against");

        try {
            commandLine = commandLineParser.parse(options, args);

            if (!commandLine.hasOption("T") || !commandLine.hasOption("C") || !commandLine.hasOption("R")){
                throw new NullPointerException("Incorrect arguments");
            }

        } catch (ParseException | NullPointerException e){
            formatter.printHelp(program + " " + version, options);
            log.log(Level.SEVERE, e.getMessage());
            System.exit(-1);
        }

        //set args
        File transcriptFile = new File (commandLine.getOptionValue("T"));
        File coverageFile = new File (commandLine.getOptionValue("C"));
        File regionsFile = new File (commandLine.getOptionValue("R"));
        int padding = commandLine.hasOption("P") ? Integer.parseInt(commandLine.getOptionValue("P")) : 5;
        int minDepth = commandLine.hasOption("D") ? Integer.parseInt(commandLine.getOptionValue("D")) : 20;

        log.log(Level.INFO, "Padding: " + padding);
        log.log(Level.INFO, "Minimum Depth: " + minDepth);
        log.log(Level.INFO, "Transcript file: " + transcriptFile.getAbsolutePath());

        //read input BED file
        try (AbstractFeatureReader bedReader = AbstractFeatureReader.getFeatureReader(regionsFile.toString(), new BEDCodec(BEDCodec.StartOffset.ONE), false)) {
            Iterable<BEDFeature> bedIter = bedReader.iterator();

            // extract overlapping GTF/GFF entries
            for (BEDFeature bedFeature : bedIter){

                TabixReader tabixReader = new TabixReader(transcriptFile.toString());
                TabixReader.Iterator results = tabixReader.query(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd());
                String line;

                //read over lapping records
                while ((line = results.next()) != null) {
                    GFF3Record gff3Record = new GFF3Record(line);
                    gff3Record.parseGFFRecord();

                }
            }

        } catch (IOException e){
            log.log(Level.SEVERE, "Coult not read BED file: " + e.getMessage());
            System.exit(-1);
        }

    }

}
