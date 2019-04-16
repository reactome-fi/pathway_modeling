/*
 * Created on Apr 7, 2011
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * Results from MCL clustering results.
 * @author wgm
 */
@XmlRootElement
public class MCLClusteringResult extends ProcessCallResult {
    // Have to specify a concrete class type for JAXB purpose.
    private ArrayList<String> clusters;
    
    public MCLClusteringResult() {
        
    }
    
    public void setClusters(ArrayList<String> clusters) {
        this.clusters = clusters;
    }
    
    
    public ArrayList<String> getClusters() {
        return this.clusters;
    }
    
}
