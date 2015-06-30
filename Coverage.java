package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by matt on 10/03/15.
 */

public class Coverage {

    private static final Logger log = Logger.getLogger(Coverage.class.getName());

    private int threshold;
    private File file;
    private ArrayList<String> sampleIDs = new ArrayList<>();
    private ListReader bedFilePathList;
    private ArrayList<HashSet<GenomicLocation>> basesPassingQC = new ArrayList<HashSet<GenomicLocation>>();

    private ArrayList<HashSet<GenomicLocation>> targetBasesByRegion = new ArrayList<HashSet<GenomicLocation>>();
    private ArrayList<ArrayList<ArrayList<GenomicLocation>>> missingBasesByRegion = new ArrayList<ArrayList<ArrayList<GenomicLocation>>>();

    private HashSet<GenomicLocation> totalTargetBases = new HashSet<>();
    private ArrayList<ArrayList<GenomicLocation>> totalMissingBases = new ArrayList<ArrayList<GenomicLocation>>();
    private ArrayList<ArrayList<GenomicLocation>> totalMissingRegions = new ArrayList<ArrayList<GenomicLocation>>();

    public Coverage(File file, ListReader bedFilePathList, int threshold){
        this.file = file;
        this.bedFilePathList = bedFilePathList;
        this.threshold = threshold;
    }

    public void populateTargetBases(){

        BEDCodec codec = new BEDCodec();

        //loop over BED files and store unique list of target bases (in total and per region)
        for (int n = 0; n < bedFilePathList.getElements().size(); ++n) {

            targetBasesByRegion.add(new HashSet<GenomicLocation>());

            //read target bases into memory
            try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(bedFilePathList.getElements().get(n), codec, false)){

                Iterable<BEDFeature> iter = reader.iterator();

                for (BEDFeature feature : iter) {
                    for (int j = feature.getStart(); j < feature.getEnd() + 1; ++j){ //TODO: check DOC : works for VCF

                        GenomicLocation temp = new GenomicLocation(feature.getContig(), j);

                        //store unique bases from a single region
                        if (!targetBasesByRegion.get(n).contains(temp)){
                            targetBasesByRegion.get(n).add(temp);
                        }

                        //store unique bases from all regions
                        if (!totalTargetBases.contains(temp)) {
                            totalTargetBases.add(temp);
                        }
                    }
                }

            } catch (IOException e){
                log.log(Level.SEVERE, e.toString());
            }

            log.log(Level.INFO, "Number of unique target bases for " + bedFilePathList.getElements().get(n) + ": " + targetBasesByRegion.get(n).size() + " bp");
        }

        log.log(Level.INFO, "Total number of unique target bases: " + totalTargetBases.size() + " bp");
    }

    public void extractBasesPassingMinThresholdFromVCF(){

        VCFFileReader vcfFileReader = new VCFFileReader(file, new File(file + ".idx"));
        Iterator<VariantContext> it = vcfFileReader.iterator();

        //store sampleIDs
        for (String sampleID : vcfFileReader.getFileHeader().getGenotypeSamples()){
            sampleIDs.add(sampleID);
            basesPassingQC.add(new HashSet<GenomicLocation>());
        }

        //read VCF file and fail bases if we are not sure on the genotype
        while(it.hasNext()) {

            VariantContext variantContext = it.next();

            //loop over sampleIDs and wildtype genotypes above minGQ
            for (int n = 0; n < vcfFileReader.getFileHeader().getNGenotypeSamples(); ++n) {

                if (variantContext.getGenotype(n).hasAnyAttribute("RGQ")){
                    if (Integer.parseInt((String) variantContext.getGenotype(n).getAnyAttribute("RGQ")) >= threshold) {
                        basesPassingQC.get(n).add(new GenomicLocation(variantContext.getContig(), variantContext.getStart()));
                    }
                } else if (variantContext.getGenotype(n).hasAnyAttribute("GQ")){
                    if ((Integer) variantContext.getGenotype(n).getAnyAttribute("GQ") >= threshold) {
                        basesPassingQC.get(n).add(new GenomicLocation(variantContext.getContig(), variantContext.getStart()));
                    }
                }

            }

        }
    }

    public void extractBasesPassingMinThresholdFromDepthOfCoverage(){

        String line;
        int lineNumber = 0;

        //read gatk per-base depth file; store bases passing QC
        try (BufferedReader depthReader = new BufferedReader(new FileReader(file))){

            while ((line = depthReader.readLine()) != null) {

                if (!line.equals("")) {
                    lineNumber++;

                    String[] fields = line.split("\t");

                    //store SampleIDs
                    if (lineNumber == 1){

                        for (int n = 3 ; n < fields.length; ++n){

                            //store sampleIDs
                            String[] subFields = fields[n].split("_");
                            sampleIDs.add(subFields[2]);

                            //initalise array
                            basesPassingQC.add(new HashSet<GenomicLocation>());

                        }

                    } else {

                        String[] subFields = fields[0].split(":");
                        GenomicLocation location = new GenomicLocation(subFields[0], Integer.parseInt(subFields[1]));

                        for (int n = 3; n < fields.length; ++n) {

                            //base passes minDP
                            if (Integer.parseInt(fields[n]) >= threshold) {
                                basesPassingQC.get(n - 3).add(location);
                            }
                        }

                    }

                }

            }

            depthReader.close();

        } catch (IOException e) {
            log.log(Level.SEVERE, "Problem reading gatk depth file: " + e.getMessage());
        }

    }

    public void extractBasesFailingMinThreshold(){

        //initalise arrays
        for (int n = 0; n < sampleIDs.size(); ++n){
            totalMissingBases.add(new ArrayList<GenomicLocation>());
            missingBasesByRegion.add(new ArrayList<ArrayList<GenomicLocation>>());

            for (int j = 0; j < bedFilePathList.getElements().size(); ++j){
                missingBasesByRegion.get(n).add(new ArrayList<GenomicLocation>());
            }

        }

        //Store bases failing QC wrt to target bed
        for (int n = 0; n < sampleIDs.size(); ++n){

            //loop over total target bases
            for (GenomicLocation loc : totalTargetBases){
                if (!basesPassingQC.get(n).contains(loc)){ //check if base passes QC
                    totalMissingBases.get(n).add(loc);
                }
            }

            //loop over target bases by target
            for (int j = 0; j < targetBasesByRegion.size(); ++j){ //loop over regions
                for (GenomicLocation loc : targetBasesByRegion.get(j)){
                    if (!basesPassingQC.get(n).contains(loc)){ //check if base passes QC
                        missingBasesByRegion.get(n).get(j).add(loc);
                    }
                }
            }

        }

    }

    public void windowMissingBases() {

        //Sort failed bases by chromosomal coordinate and condense into regions
        for (int n = 0; n < sampleIDs.size(); ++n) { //loop over SampleIDs
            totalMissingRegions.add(GenomicLocation.mergeBases(totalMissingBases.get(n)));
        }

    }

    public void printCoverageMetrics(){

        String setName = bedFilePathList.getFileName().substring(0, bedFilePathList.getFileName().lastIndexOf("."));

        //loop over SampleIDs
        for (int n = 0; n < sampleIDs.size(); ++n){

            //print gaps
            try (PrintWriter writer = new PrintWriter(sampleIDs.get(n) + "_" + setName + "_Gaps.bed")){

                for (GenomicLocation location : totalMissingRegions.get(n)){
                    writer.println(location.getChromosome() + "\t" + (location.getStartPosition() - 1) + "\t" + location.getEndPosition()); //convert to BED 0-based
                }

                writer.close();
            } catch (IOException e){
                log.log(Level.SEVERE, "Could not write to file: " + e.getMessage());
            }

            //print percentage coverage
            try (PrintWriter writer = new PrintWriter(sampleIDs.get(n) + "_" + setName + "_PercentageCoverage.txt")){

                //print headers
                writer.print("SampleID\t" + setName); //set

                for (int j = 0; j < bedFilePathList.getElements().size(); ++j){ //loop over beds belonging to this set

                    File targetPath = new File(bedFilePathList.getElements().get(j));
                    String targetName = targetPath.getName().substring(0, targetPath.getName().lastIndexOf("."));

                    writer.print("\t" + targetName);
                }
                writer.println();

                //print metrics
                writer.write(sampleIDs.get(n));
                writer.print("\t" + String.format("%.4g", (double) (totalTargetBases.size() - totalMissingBases.get(n).size()) / totalTargetBases.size())); //set

                //print metrics for each target
                for (int k = 0; k < bedFilePathList.getElements().size(); ++k){ // target
                    writer.print("\t" + String.format("%.4g", (double) (targetBasesByRegion.get(k).size() - missingBasesByRegion.get(n).get(k).size()) / targetBasesByRegion.get(k).size()));
                }

                writer.println();

                writer.close();
            } catch (IOException e) {
                log.log(Level.SEVERE, "Could not write to file: " + e.getMessage());
            }
        }

    }

}
