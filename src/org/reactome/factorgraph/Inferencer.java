/*
 * Created on Oct 27, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.Map;


/**
 * An interface that defines what methods should be implemtented for a class doing inference.
 * @author gwu
 *
 */
public interface Inferencer {
    
    public void setFactorGraph(FactorGraph fg);
    
    public FactorGraph getFactorGraph();
    
    public <T extends Number> void setObservation(Observation<T> observation);
    
    public <T extends Number> void setObservation(Map<Variable, T> varToAssignment);
    
    public void clearObservation();
    
    public Observation<? extends Number> getObservation();
    
    public void runInference() throws InferenceCannotConvergeException;
    
    public double calculateLogZ();
    
}
