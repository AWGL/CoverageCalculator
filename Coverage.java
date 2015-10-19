package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by matt on 10/03/15.
 */

public class Coverage {

    private static final Logger log = Logger.getLogger(Coverage.class.getName());

    private int threshold;
    private File coverageFile, bedFile;

    private ArrayList<String> sampleIds = new ArrayList<>();

    private HashMap<String, HashSet<GenomicLocation>> targetBasesByGene = new HashMap<>(); //Gene:[single-base pos]
    private HashSet<GenomicLocation> totalTargetBases = new HashSet<>(); //[all single base pos]

    private ArrayList<HashSet<GenomicLocation>> totalPassingBasesBySample = new ArrayList<>();
    private ArrayList<ArrayList<GenomicLocation>> totalMissingBasesBySample = new ArrayList<>();
    private ArrayList<HashMap<String, ArrayList<GenomicLocation>>> geneMissingBasesBySample = new ArrayList<>();

    private ArrayList<ArrayList<GenomicLocation>> totalMissingRegionsBySample = new ArrayList<>();

    public Coverage(File coverageFile, File bedFile, int threshold){
        this.coverageFile = coverageFile;
        this.bedFile = bedFile;
        this.threshold = threshold;
    }
    public void populateTargetBases() {

        log.log(Level.INFO, "Extracting target bases");

        //read target bases into memory
        try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(bedFile.toString(), new BEDCodec(BEDCodec.StartOffset.ONE), false)) {
            Iterable<BEDFeature> iter = reader.iterator();

            for (BEDFeature feature : iter) {

                String[] fields = feature.getName().split(":");

                if (!targetBasesByGene.containsKey(fields[0])) targetBasesByGene.put(fields[0], new HashSet<GenomicLocation>());

                //split bed feature into single base coordinates
                for (int j = feature.getStart(); j < feature.getEnd() + 1; ++j) { //convert to single base 1-pos coordinates
                    GenomicLocation genomicLocation = new GenomicLocation(feature.getContig(), j);

                    //store unique bases from a single region
                    targetBasesByGene.get(fields[0]).add(genomicLocation);

                    //store all target bases
                    totalTargetBases.add(genomicLocation);
                }

            }

            reader.close();
        } catch (IOException e) {
            log.log(Level.SEVERE, e.toString());
        }

    }
    public void extractBasesPassingMinDepth(){

        log.log(Level.INFO, "Extracting bases passing minimum depth");

        boolean passHeaderLine = false;
        String line;

        //read gatk per-base depth file; store bases passing QC
        try (BufferedReader reader = new BufferedReader(new FileReader(coverageFile))){

            while ((line = reader.readLine()) != null) {

                if (!line.equals("")) {
                    String[] fields = line.split("\t");

                    //store SampleIDs
                    if (!passHeaderLine){

                        for (int n = 3 ; n < fields.length; ++n){

                            //store sampleIDs
                            String[] subFields = fields[n].split("_");
                            sampleIds.add(subFields[2]);

                            //initalise array
                            totalPassingBasesBySample.add(new HashSet<GenomicLocation>());
                            totalMissingBasesBySample.add(new ArrayList<GenomicLocation>());
                            geneMissingBasesBySample.add(new HashMap<String, ArrayList<GenomicLocation>>());

                        }

                        passHeaderLine = true;

                    } else {

                        String[] subFields = fields[0].split(":");
                        GenomicLocation genomicLocation = new GenomicLocation(subFields[0], Integer.parseInt(subFields[1]));

                        for (int n = 3; n < fields.length; ++n) {

                            //base passes minDP
                            if (Integer.parseInt(fields[n]) >= threshold) {
                                totalPassingBasesBySample.get(n - 3).add(genomicLocation);
                            }
                        }

                    }

                }

            }

            reader.close();
        } catch (IOException e) {
            log.log(Level.SEVERE, "Problem reading gatk depth file: " + e.getMessage());
        }

    }
    public void extractBasesFailingMinDepth(){ //todo optimise

        log.log(Level.INFO, "Extracting bases failing minimum depth");

        //Store bases failing QC wrt to targets
        for (int n = 0; n < sampleIds.size(); ++n){

            //total targets
            for (GenomicLocation genomicLocation : totalTargetBases){
                if (!totalPassingBasesBySample.get(n).contains(genomicLocation)){ //check if base passes QC
                    totalMissingBasesBySample.get(n).add(genomicLocation);
                }
            }

            //gene targets
            for (Map.Entry<String, HashSet<GenomicLocation>> geneSet : targetBasesByGene.entrySet()){

                //initalise
                geneMissingBasesBySample.get(n).put(geneSet.getKey(), new ArrayList<GenomicLocation>());

                for (GenomicLocation genomicLocation : geneSet.getValue()){

                    if (!totalPassingBasesBySample.get(n).contains(genomicLocation)){ //check if base passes QC
                        geneMissingBasesBySample.get(n).get(geneSet.getKey()).add(genomicLocation);
                    }

                }
            }

        }

    }
    public void windowMissingBases() {

        log.log(Level.INFO, "Converting single base coordinates to regions");

        //Sort failed bases by chromosomal coordinate and condense into regions
        for (int n = 0; n < sampleIds.size(); ++n) {
            totalMissingRegionsBySample.add(GenomicLocation.mergeBases(totalMissingBasesBySample.get(n)));
        }

    }
    public void printCoverageMetrics(){

        log.log(Level.INFO, "Printing coverage metrics");

        for (int n = 0; n < sampleIds.size(); ++n){

            //print gaps
            try (PrintWriter writer = new PrintWriter(sampleIds.get(n) + "_Gaps.bed")){

                for (GenomicLocation genomicLocation : totalMissingRegionsBySample.get(n)){
                    writer.println(genomicLocation.getContig() + "\t" + (genomicLocation.getStartPosition() - 1) + "\t" + genomicLocation.getEndPosition()); //convert to BED 0-based
                }

                writer.close();
            } catch (IOException e){
                log.log(Level.SEVERE, "Could not write to file: " + e.getMessage());
            }

        }

        //print percentage coverage
        try (PrintWriter writer = new PrintWriter("PercentageCoverage.txt")){

            //print headers
            writer.print("SampleId");

            for (Map.Entry<String, HashSet<GenomicLocation>> geneSet : targetBasesByGene.entrySet()){
                writer.print("\t" + geneSet.getKey());
            }

            writer.println();

            //print % coverage by Gene
            for (int n = 0; n < sampleIds.size(); ++n){

                //print metrics
                writer.print(sampleIds.get(n));

                for (Map.Entry<String, HashSet<GenomicLocation>> geneSet : targetBasesByGene.entrySet()){
                    writer.print("\t" + String.format("%.4g", ((double) (geneSet.getValue().size() - geneMissingBasesBySample.get(n).get(geneSet.getKey()).size()) / geneSet.getValue().size()) * 100));
                }

                writer.println();
            }

            writer.close();
        } catch (IOException e) {
            log.log(Level.SEVERE, "Could not write to percentage coverage file: " + e.getMessage());
        }

    }

}
