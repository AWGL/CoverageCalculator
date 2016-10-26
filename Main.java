package nhs.genetics.cardiff;

import nhs.genetics.cardiff.framework.*;
import nhs.genetics.cardiff.framework.GFF3Record;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * Program for reporting coverage metrics with gene context
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-03-15
 */
public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final String version = "2.0.1";

    public static void main(String[] args) {

        int minDepth = 30;
        int padding = 20;

        if (args.length < 3 || args.length > 5){
            System.err.println("CoverageCalculator v" + version);
            System.err.println("Usage: <PerBaseCoverage> <Gene/ENSTList> <GTF/GFF>");
            System.err.println("-d Minimum depth [" + minDepth + "]");
            System.err.println("-p Interval padding [" + padding + "]");
            System.exit(1);
        }

        //get cmd line args
        for (String arg : args){
            if (Pattern.matches("^-d.*", arg)){
                String[] fields = arg.split("d");
                minDepth = Integer.parseInt(fields[1]);
            }
            if (Pattern.matches("^-p.*", arg)){
                String[] fields = arg.split("p");
                padding = Integer.parseInt(fields[1]);
            }
        }

        log.log(Level.INFO, "Padding: " + padding);
        log.log(Level.INFO, "Minimum Depth: " + minDepth);

        boolean gtfOrGff;
        String line, id;
        HashMap<String, HashSet<GenomicLocation>> targetGenomicLocations = new HashMap<>();
        GatkDepthOfCoverageParser gatkDepthOfCoverageParser = null;

        //determine transcript file format
        String[] fields = args[2].split("\\.");
        if (fields[fields.length - 1].toLowerCase().equals("gtf")){
            gtfOrGff = true;
            log.log(Level.INFO, "Determined transcript file as GTF");
        } else if (fields[fields.length - 1].toLowerCase().equals("gff3")){
            gtfOrGff = false;
            log.log(Level.INFO, "Determined transcript file as GFF3");
        } else {
            throw new IllegalArgumentException("Cannot detemine transcript file format with extension: " + fields[fields.length - 1].toLowerCase());
        }

        //parse gene list
        ListReader listReader = new ListReader(new File(args[1]));

        try {
            listReader.parseListReader();
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read gene list " + e.getMessage());
            System.exit(1);
        }

        //initalise map
        for (String element : listReader.getUniqueElements()){
            targetGenomicLocations.put(element, new HashSet<GenomicLocation>()); //1-based coordinates
        }

        //create target list of bases; read GTF/GFF
        log.log(Level.INFO, "Extracting target regions of interest");
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(args[2])))){
            while ((line = reader.readLine()) != null) {

                //skip headers
                if (!Pattern.matches("^#.*", line))
                {
                    if (gtfOrGff){

                        //split record and add to array
                        GTFRecord gtfRecord = new GTFRecord(line);
                        gtfRecord.parseGTFRecord();

                        // ROI
                        if ((gtfRecord.getFeature().equals("CDS") || gtfRecord.getFeature().equals("stop_codon")) && gtfRecord.getAttributes().get("transcript_biotype").equals("protein_coding")){

                            if (listReader.getUniqueElements().contains(gtfRecord.getAttributes().get("gene_name"))){
                                id = gtfRecord.getAttributes().get("gene_name");
                            } else if (listReader.getUniqueElements().contains(gtfRecord.getAttributes().get("transcript_id"))){
                                id = gtfRecord.getAttributes().get("transcript_id");
                            } else {
                                continue;
                            }

                            //split regions feature into single base coordinates
                            for (int j = gtfRecord.getGenomicLocation().getStartPosition() - padding; j < (gtfRecord.getGenomicLocation().getEndPosition()+1) + padding; ++j) {
                                targetGenomicLocations.get(id).add(new GenomicLocation(gtfRecord.getGenomicLocation().getContig(), j)); //1-based
                            }

                        }

                    } else {

                        //split record and add to array
                        GFF3Record gff3Record = new GFF3Record(line);
                        gff3Record.parseGFFRecord();

                        // ROI
                        if (gff3Record.getType().equals("CDS")){

                            if (listReader.getUniqueElements().contains(gff3Record.getAttributes().get("gene"))){
                                id = gff3Record.getAttributes().get("gene");
                            } else if (listReader.getUniqueElements().contains(gff3Record.getAttributes().get("transcript_id"))){
                                id = gff3Record.getAttributes().get("transcript_id");
                            } else {
                                continue;
                            }

                            //split regions feature into single base coordinates
                            for (int j = gff3Record.getGenomicLocation().getStartPosition() - padding; j < (gff3Record.getGenomicLocation().getEndPosition()+1) + padding; ++j) {
                                targetGenomicLocations.get(id).add(new GenomicLocation(gff3Record.getGenomicLocation().getContig(), j)); //1-based
                            }

                        }

                    }
                }

            }

            reader.close();

        } catch (IOException e) {
            log.log(Level.SEVERE, "Could not read GTF/GFF file: " + e.getMessage());
            System.exit(1);
        }

        //print missing ids
        for (String element : listReader.getUniqueElements()){
            if (targetGenomicLocations.get(element).size() == 0){
                log.log(Level.WARNING, "Could not find: " + element);
            }
        }

        /*print ROI for debugging
        try (PrintWriter printWriter = new PrintWriter("debug.bed")) {
            for (GenomicLocation location : GenomicLocation.mergeBases(new ArrayList<>(targetGenomicLocations.get("BRCA1")))){
                printWriter.println(location.getContig() + "\t" + (location.getStartPosition() - 1) + "\t" + location.getEndPosition());
            }
        } catch (IOException e){
            log.log(Level.SEVERE, "could not write debug file");
        }*/

        /*try (PrintWriter printWriter = new PrintWriter("debug.bed")) {
            for (Map.Entry<String, HashSet<GenomicLocation>> it : targetGenomicLocations.entrySet()){
                for (GenomicLocation genomicLocation : it.getValue()){
                    printWriter.println(genomicLocation.getConcatenatedLocation());
                }
            }
        } catch (IOException e){
            log.log(Level.SEVERE, "could not write debug file");
            System.exit(1);
        }*/


        //read GATK depth of coverage file
        try {
            gatkDepthOfCoverageParser = new GatkDepthOfCoverageParser(new File(args[0])); //1-based
            gatkDepthOfCoverageParser.parseGatkDepthOfCoverageFile();
        } catch (IOException e){
            log.log(Level.SEVERE, "could not read coverage file: " + e.getMessage());
            System.exit(1);
        }

        //loop over samples and report coverage metrics
        for (int n = 0; n < gatkDepthOfCoverageParser.getSampleIds().size(); ++n){
            log.log(Level.INFO, "Generating coverage metrics for sample: " + gatkDepthOfCoverageParser.getSampleIds().get(n));

            try (PrintWriter printWriter = new PrintWriter(gatkDepthOfCoverageParser.getSampleIds().get(n) + "_gaps.bed")) {

                for (Map.Entry<String, HashSet<GenomicLocation>> target : targetGenomicLocations.entrySet()){
                    log.log(Level.FINE, "Generating coverage metrics for target: " + target.getKey());

                    //calculate coverage from per base coverage file
                    Coverage coverage = new Coverage(target.getValue(), gatkDepthOfCoverageParser.getDepthOfCoverage().get(n), minDepth);
                    coverage.populateCoverageMetics();

                    System.out.println(gatkDepthOfCoverageParser.getSampleIds().get(n) + "\t" + target.getKey() + "\t" + String.format("%.2f", coverage.getCoveragePercentage() * 100) + "%");

                    for (GenomicLocation location : coverage.getWindowMissingBases()){
                        printWriter.println(location.getContig() + "\t" + (location.getStartPosition() - 1) + "\t" + location.getEndPosition()); //print to bed; convert to 0-based on the fly
                    }

                }

            } catch (IOException e){
                log.log(Level.SEVERE, "Could not write gaps for: " + gatkDepthOfCoverageParser.getSampleIds().get(n));
                System.exit(1);
            }

        }

    }
}
