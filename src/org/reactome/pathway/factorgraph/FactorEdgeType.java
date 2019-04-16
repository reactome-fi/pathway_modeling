/*
 * Created on May 23, 2012
 *
 */
package org.reactome.pathway.factorgraph;


/**
 * This enum lists types of edges between a factor node and a variable node in the converted
 * factor graph from Reactome pathway.
 * @author gwu
 *
 */
public enum FactorEdgeType {
    // types for reactions
    INPUT,
    CATALYST,
    ACTIVATOR,
    INHIBITOR,
    OUTPUT,
    // Complex
    COMPLEX,
    // EntitySet
    MEMBER,
    // Observation
    OBSERVATION;
    
    static enum FactorEdgeLabel {
        MIN,
        MAX
    }
    
    public static FactorEdgeLabel mapTypeToLabel(FactorEdgeType type) {
        switch (type) {    
            case INPUT :
            case CATALYST :
            case ACTIVATOR :
            case COMPLEX :
                return FactorEdgeLabel.MIN;
            case MEMBER :
                return FactorEdgeLabel.MAX;
        }
        return null;
    }
    
    public static int getTypeWeight(FactorEdgeType type) {
        switch (type) {
            // Assign different weights for edge types
//            case INPUT : return 1;
//            case CATALYST : return 3;
//            case ACTIVATOR : return 2;
//            case INHIBITOR : return 2;
//            default : return 1;
            // Assign same weight for edge types
            case INPUT : return 1;
            case CATALYST : return 1;
            case ACTIVATOR : return 1;
            case INHIBITOR : return 1;
            default : return 1;
        }
    }
    
    /**
     * If there are multiple inputs, actiavators or inhibitors, a group node will be created.
     * This method is used to provide the edge types for edges between individual entities and
     * these group nodes.
     * @param type
     * @return
     */
    public static FactorEdgeType mapGroupTypeToSingleType(FactorEdgeType type) {
        if (type == ACTIVATOR || type == INHIBITOR)
            return MEMBER;
        return type;
    }
}
