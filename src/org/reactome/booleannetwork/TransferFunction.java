/*
 * Created on Sep 6, 2018
 *
 */
package org.reactome.booleannetwork;

/**
 * Transfer function as defined in the original logic model simulation.
 * @author wug
 *
 */
public interface TransferFunction {

    public Double transfer(Number value);
    
}
