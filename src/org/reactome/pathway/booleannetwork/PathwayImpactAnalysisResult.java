/*
 * Created on Sep 14, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.List;

/**
 * A model class to store impact analysis results based on Boolean network simulation.
 * @author wug
 *
 */
public class PathwayImpactAnalysisResult {
    
    private Long dbId;
    private String pathwayName;
    private Double sum;
    private Double average;
    private List<String> targetGenes;
    
    public PathwayImpactAnalysisResult() {
    }

    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public String getPathwayName() {
        return pathwayName;
    }

    public void setPathwayName(String pathwayName) {
        this.pathwayName = pathwayName;
    }

    public Double getSum() {
        return sum;
    }

    public void setSum(Double sum) {
        this.sum = sum;
    }

    public Double getAverage() {
        return average;
    }

    public void setAverage(Double average) {
        this.average = average;
    }

    public List<String> getTargetGenes() {
        return targetGenes;
    }

    public void setTargetGenes(List<String> targetGenes) {
        this.targetGenes = targetGenes;
    }
    
    @Override
    public String toString() {
        return dbId + "\t" + pathwayName + "\t" + sum + "\t" + average + "\t" + String.join(",", targetGenes);
    }

}
