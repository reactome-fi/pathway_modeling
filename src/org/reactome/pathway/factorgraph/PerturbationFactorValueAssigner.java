/*
 * Created on Jan 25, 2017
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.List;

import org.gk.util.StringUtils;
import org.junit.Test;
import org.reactome.factorgraph.Variable;

/**
 * This customized FactorValueAssigner is used to assign values for factors using two-state
 * variables: 0 for normal and 1 for perturbed.
 * @author gwu
 *
 */
public class PerturbationFactorValueAssigner extends FactorValueAssigner {
    // This value is 1/10 of minumum value, excluding 0, from the factor graph converted from 
    // PIP3 activates AKT signaling when beta = 5.0d
//    private final static double MINIMUM_VALUE = 0.0005d;
//    private final static double BETA = 5.0d;
    
    // Using the same pathway and beta = 3, 1/10 of the minimum value of factor functions
//    private final static double MINIMUM_VALUE = 0.005d;
//    private final double BETA = 3.0d;
    
    // Using beta = 10
    private final static double MINIMUM_VALUE = 4.5e-5 / 10.0d;
    private final double BETA = 10.0d;
    
    private final static double MAXIMUM_VALUE = 1.0d - MINIMUM_VALUE;
    
    /**
     * Default constructor.
     */
    public PerturbationFactorValueAssigner() {
    }
    
    @Test
    public void testTypicalReactionFactorValues() {
        DEBUG = true;
        List<Variable> variables = new ArrayList<Variable>();
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        // Check reaction factor
        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.CATALYST, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.ACTIVATOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.INHIBITOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        
        // Test a binding or a complex formation
//        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
//        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
//        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
//        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        
        // Test an EntitySet
//        testGenerateVariable(FactorEdgeType.MEMBER, variables, edgeTypes);
//        testGenerateVariable(FactorEdgeType.MEMBER, variables, edgeTypes);
//        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        
        generateFactorValues(variables, edgeTypes);
    }
    
    @Override
    public List<Double> generateFactorValues(List<Variable> variables,
                                             List<FactorEdgeType> edgeTypes) {
        List<Double> values = new ArrayList<Double>();
        List<Integer> states = new ArrayList<Integer>();
        for (Variable var : variables)
            states.add(0); // State with state == 0
        boolean isChanged = true;
        int count = 0;
        if (DEBUG)
            System.out.println(StringUtils.join("\t", variables));
        while (isChanged) {
            double expectedState = calculateExpectedOutputState(states, edgeTypes);
            int actualOutputState = getState(states, 
                                             edgeTypes,
                                             FactorEdgeType.OUTPUT);
            
            double value = 0.0d;
            // Test code to make sure the output is used as a directed probability
            if (actualOutputState == 1) {
                value = 1.0d - Math.exp(-expectedState * BETA);
                if (value < MINIMUM_VALUE)
                    value = MINIMUM_VALUE;
                else if (value > MAXIMUM_VALUE)
                    value = MAXIMUM_VALUE;
            }
            else {
                double diff = Math.abs(expectedState - actualOutputState);
                // Use an exponential function
                value = Math.exp(-diff * BETA);
                if (value > MAXIMUM_VALUE)
                    value = MAXIMUM_VALUE;
                else if (value < MINIMUM_VALUE)
                    value = MINIMUM_VALUE;
            }
            
            values.add(value);
            if (DEBUG) {
                System.out.println(StringUtils.join("\t", states) + "\t" + value);
            }
            // Change the state
            isChanged = false;
            for (int i = 0; i < states.size(); i++) {
                int state = states.get(i);
                Variable var = variables.get(i);
                if (state < var.getStates() - 1) {
                    states.set(i, ++state);
                    // Need to set previous states back
                    for (int j = i - 1; j >= 0; j--) {
                        states.set(j, 0); 
                    }
                    isChanged = true;
                    break;
                }
            }
        }
        ensureValues(values);
        return values;
    }
    
    @Override
    protected double calculateExpectedOutputState(List<Integer> states,
                                                  List<FactorEdgeType> edgeTypes) {
        if (isMemberEdge(edgeTypes))
            return calculateExpectedOutputStateForMembers(states, edgeTypes);
//        return getMaximumOutputState(states, edgeTypes);
        int total = 0;
        int totalWeight = 0;
        for (int i = 0; i < edgeTypes.size(); i++) {
            FactorEdgeType edgeType = edgeTypes.get(i);
            if (edgeType == FactorEdgeType.OUTPUT) {
                continue;
            }
            int state = states.get(i);
            int weight = FactorEdgeType.getTypeWeight(edgeType);
            total += state * weight; 
            totalWeight += weight;
        }
        // Arithematic mean
        double expectedState = (double) total / totalWeight;
        return expectedState;
    }
    
    private double getMaximumOutputState(List<Integer> states,
                                         List<FactorEdgeType> edgeTypes) {
        for (int i = 0; i < edgeTypes.size(); i++) {
            FactorEdgeType edgeType = edgeTypes.get(i);
            if (edgeType == FactorEdgeType.OUTPUT) {
                continue;
            }
            int state = states.get(i);
            if (state > 0)
                return state;
        }
        return 0.0d;
    }
    
    /**
     * Search if there is a normal state in the member list. If yes, return it.
     * @param states
     * @param edgeTypes
     * @return
     */
    private double calculateExpectedOutputStateForMembers(List<Integer> states,
                                                          List<FactorEdgeType> edgeTypes) {
        int rtn = 1; // Perturbed state
        for (int i = 0; i < edgeTypes.size(); i++) {
            FactorEdgeType edgeType = edgeTypes.get(i);
            if (edgeType == FactorEdgeType.OUTPUT)
                continue;
            int state = states.get(i);
            if (state < rtn)
                rtn = state;
        }
        return rtn;
    }
    
    private boolean isMemberEdge(List<FactorEdgeType> types) {
        boolean isMemberType = false;
        boolean hasOtherType = false;
        // Make sure only one input type: Member
        for (FactorEdgeType type : types) {
            if (type == FactorEdgeType.OUTPUT)
                continue;
            if (type == FactorEdgeType.MEMBER)
                isMemberType = true;
            else {
                hasOtherType = true;
            }
        }
        if (isMemberType && hasOtherType)
            throw new IllegalArgumentException("Member edge type should not be mixed with other type!");
        return isMemberType;
    }
    
}
