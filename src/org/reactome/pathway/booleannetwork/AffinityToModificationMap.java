/*
 * Created on Sep 7, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

/**
 * This interface is used to map from affinity (e.g. binding constants) to modification strength
 * used for Boolean network simulation.
 * @author wug
 *
 */
public interface AffinityToModificationMap {

    /**
     * Map from an affinity value to a modification strength.
     * @param value
     * @return
     */
    public double getModificationStrenth(double value);
    
}
