/*
 * Created on Oct 14, 2014
 *
 */
package org.reactome.pathway.factorgraph;

/**
 * @author gwu
 *
 */
public class PathwayFGConstants {
    
    // This is used to control if a complex or entityset should be split into
    // multiple factor nodes. For example, a complex has 9 subunits, it will generate
    // a double array containing 59,049 numbers. If we split it into 3 groups (randomly),
    // we will get 324 numbers only (4 * 81).
    public static final int MAXIMUM_AUGUMENT_NODE = 3;
    // The states for each Var. Three states right now: normal, down, up
    public static int NUMBER_OF_STATES = 3;
    // Two states for variable: 0 for normal and 1 for perturbed
//    public static int NUMBER_OF_STATES = 2;
    // Used as epison value: see the original PARADIM implementation
    public static final double EPSILON_VALUE = 0.001;
    
}
