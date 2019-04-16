/*
 * Created on Oct 23, 2014
 *
 */
package org.reactome.factorgraph;

/**
 * If LoopyBeliefPropagation cannot converge, this Exception will be thrown.
 * @author gwu
 *
 */
public class InferenceCannotConvergeException extends Exception {
    
    /**
     * Default constructor.
     */
    public InferenceCannotConvergeException() {
    }
 
    public InferenceCannotConvergeException(String message) {
        super(message);
    }
    
}
