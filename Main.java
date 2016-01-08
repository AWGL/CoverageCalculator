package nhs.genetics.cardiff;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        int minDepth = 30;
        int padding = 20;

        if (args.length != 3){
            System.err.println("Usage: <PerBaseCoverage> <Gene/ENSTList> <GTF>");
            System.err.println("-d Minimum depth [" + minDepth + "]");
            System.err.println("-p Interval padding [" + padding + "]");
            System.exit(1);
        }

        String line, id;
        HashMap<String, HashSet<GenomicLocation>> targetGenomicLocations = new HashMap<>();

        //parse gene list
        ListReader listReader = new ListReader(new File(args[1]));
        listReader.parseListReader();

        //initalise map
        for (String element : listReader.getUniqueElements()){
            targetGenomicLocations.put(element, new HashSet<GenomicLocation>());
        }

        //create target list of bases; read GTF
        log.log(Level.INFO, "Extracting target regions of interest");
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(args[2])))){
            while ((line = reader.readLine()) != null) {

                //skip headers
                if (!Pattern.matches("^#.*", line))
                {
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
                        for (int j = gtfRecord.getGenomicLocation().getStartPosition() - padding; j < gtfRecord.getGenomicLocation().getEndPosition() + padding; ++j) {
                            targetGenomicLocations.get(id).add(new GenomicLocation(gtfRecord.getGenomicLocation().getContig(), j));
                        }

                    }

                }

            }

            reader.close();

        } catch (IOException e) {
            log.log(Level.SEVERE, "Could not read GTF file: " + e.getMessage());
        }

        //print missing ids
        for (String element : listReader.getUniqueElements()){
            if (targetGenomicLocations.get(element).size() == 0){
                log.log(Level.WARNING, "Could not find: " + element);
            }
        }

        //print ROI for debugging
        try (PrintWriter printWriter = new PrintWriter("debug.bed")) {
            for (GenomicLocation location : GenomicLocation.mergeBases(new ArrayList<>(targetGenomicLocations.get("BRCA1")))){
                printWriter.println(location.getContig() + "\t" + (location.getStartPosition() - 1) + "\t" + location.getEndPosition());
            }
        } catch (IOException e){
            log.log(Level.SEVERE, "could not write debug file");
        }

        //read GATK depth of coverage file
        GatkDepthOfCoverageParser gatkDepthOfCoverageParser = new GatkDepthOfCoverageParser(new File(args[0]));
        gatkDepthOfCoverageParser.parseGatkDepthOfCoverageFile();

        //loop over samples and report coverage metrics
        for (int n = 0; n < gatkDepthOfCoverageParser.getSampleIds().size(); ++n){
            log.log(Level.INFO, "Generating coverage metrics for sample: " + gatkDepthOfCoverageParser.getSampleIds().get(n));

            try (PrintWriter printWriter = new PrintWriter(gatkDepthOfCoverageParser.getSampleIds().get(n) + "_gaps.bed")) {

                for (Map.Entry<String, HashSet<GenomicLocation>> target : targetGenomicLocations.entrySet()){
                    log.log(Level.INFO, "Generating coverage metrics for target: " + target.getKey());

                    //calculate coverage from per base coverage file
                    Coverage coverage = new Coverage(target.getValue(), gatkDepthOfCoverageParser.getDepthOfCoverage().get(n), minDepth);
                    coverage.populateCoverageMetics();

                    System.out.println(gatkDepthOfCoverageParser.getSampleIds().get(n) + "\t" + target.getKey() + "\t" + coverage.getCoveragePercentage());

                    for (GenomicLocation location : coverage.getWindowMissingBases()){
                        printWriter.println(location.getContig() + "\t" + location.getStartPosition() + "\t" + location.getEndPosition());
                    }

                }

            } catch (IOException e){
                log.log(Level.SEVERE, "Could not write gaps for: " + gatkDepthOfCoverageParser.getSampleIds().get(n));
            }

        }

    }
}
