/*
 * Created on Sep 14, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * This interface is used to map a drug to its targets as a map from genes to inhibition strengths.
 * @author wug
 *
 */
public interface DrugToTargetsMapper {
    
    public Map<String, Double> getGeneToInhibition(Collection<String> genes, 
                                                   String drug) throws Exception;
    
    public Map<String, Double> getGeneToActivation(Collection<String> genes,
                                                   String drug) throws Exception;
    
    public Set<String> getDrugTargets(String drug) throws Exception;

}
