/*
 * Created on Dec 2, 2014
 *
 */
package org.reactome.fi.pgm;

import java.util.Map;
import java.util.Set;

/**
 * This interface is used to determine the default factor function for a pair-wise
 * Markov Random Field.
 * @author gwu
 *
 */
public interface FIPGMFactorValuesAssigner {
    
    /**
     * Define the way to convert a FI into a factor.
     * @param fi
     * @param nameToPartners
     * @return
     */
    public double[] getValuesForFIFactor(String fi,
                                        Map<String, Set<String>> nameToPartners);
    
    /**
     * Define the way to convert a single gene in a FI network into a factor.
     * @param gene
     * @return
     */
    public double[] getValuesForGeneFactor();
    
    /**
     * Return true if a gene factor should be created. If true is returned in this
     * method, getValuesForGeneFactor() should return a non-null value. 
     * @return
     */
    public boolean isGeneFactorValueSupported();
    
}
