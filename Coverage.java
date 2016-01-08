package nhs.genetics.cardiff;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by matt on 10/03/15.
 */

public class Coverage {

    private int minDepth;
    private HashSet<GenomicLocation> targetGenomicLocations;
    private HashMap<GenomicLocation, Integer> depthOfCoverage;
    private ArrayList<GenomicLocation> missingBases = new ArrayList<>();
    private int passedBases = 0;

    public Coverage(HashSet<GenomicLocation> targetGenomicLocations, HashMap<GenomicLocation, Integer> depthOfCoverage, int minDepth){
        this.minDepth = minDepth;
        this.targetGenomicLocations = targetGenomicLocations;
        this.depthOfCoverage = depthOfCoverage;
    }

    public void populateCoverageMetics(){
        for (GenomicLocation genomicLocation : targetGenomicLocations){
            if (depthOfCoverage.containsKey(genomicLocation) && depthOfCoverage.get(genomicLocation) >= minDepth){
                passedBases++;
            } else {
                missingBases.add(genomicLocation);
            }
        }
    }

    public ArrayList<GenomicLocation> getWindowMissingBases() {
        return GenomicLocation.mergeBases(missingBases);
    }
    public int getPassedBases() {
        return passedBases;
    }
    public double getCoveragePercentage(){
        return (double) passedBases / targetGenomicLocations.size();
    }

}
