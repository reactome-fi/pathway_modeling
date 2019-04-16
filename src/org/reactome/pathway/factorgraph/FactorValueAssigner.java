/*
 * Created on May 22, 2012
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.gk.util.StringUtils;
import org.junit.Test;
import org.reactome.factorgraph.Variable;
import org.reactome.pathway.factorgraph.FactorEdgeType.FactorEdgeLabel;

/**
 * This interface is used to assign values for a factor. See the documentation to class PathwayToFactorGraphConverter
 * for details about the implementation of this class.
 * @author gwu
 *
 */
public class FactorValueAssigner {
    protected double major;
    protected double minor;
    protected static boolean DEBUG = false;
    
    public FactorValueAssigner() {
        major = 1.0d - PathwayFGConstants.EPSILON_VALUE;
        minor = PathwayFGConstants.EPSILON_VALUE / 2.0d;
    }
    
    public List<Double> generateFactorValues(List<Variable> variables,
                                             List<FactorEdgeType> edgeTypes) {
        // In our implementation, there should be only one edge label
        Set<FactorEdgeLabel> edgeLabels = new HashSet<FactorEdgeLabel>();
        for (FactorEdgeType type : edgeTypes) {
            FactorEdgeLabel label = FactorEdgeType.mapTypeToLabel(type);
            if (label != null)
                edgeLabels.add(label);
        }
        if (edgeLabels.size() > 1) {
            throw new IllegalArgumentException("More than one edge labels in the passed EdgeTypes!");
        }
        List<Double> values = new ArrayList<Double>();
        List<Integer> states = new ArrayList<Integer>();
        for (Variable var : variables)
            states.add(0); // State with state == 0
        boolean isChanged = true;
        int count = 0;
        while (isChanged) { // Need to go through each combination
            double value = calculateFactorValue(states, edgeTypes);
//            double value = calculateFactorValueViaDistribution(states, edgeTypes);
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
    
    /**
     * Make sure there is no NaN in the value
     * @param values
     */
    protected void ensureValues(List<Double> values) {
        for (Double value : values)
            if (Double.isNaN(value))
                throw new IllegalStateException("Generated values has NaN: " + values);
    }
    
    /**
     * A new way to calculate a factor value using normal distribution by checking the difference
     * between expected state and the actual state of the output. The expected value will be calculated
     * based on a weighted average from the left-hand-side of a reaction, including input, catalyst,
     * inhibitor, and activator. Weights are pre-assigned and hard-coded in class FactorEdgeType.
     * @return
     */
    private double calculateFactorValueViaDistribution(List<Integer> states,
                                                       List<FactorEdgeType> edgeTypes) {
        // Calculate expected state of the output
        double expectedState = calculateExpectedOutputState(states, edgeTypes);
        int actualOutputState = getState(states, 
                                         edgeTypes,
                                         FactorEdgeType.OUTPUT);
        // Use normal distribution to calculate p-value on how confidence we can
        // assume the correct of the expected state
        NormalDistribution normal = new NormalDistribution();
        double diff = Math.abs(expectedState - actualOutputState);
        // Use two-sides p-value
        double value = (1.0 - normal.cumulativeProbability(diff)) * 2;
        // Give it some buffer and don't assign 1.0 for factor value
        if (value == 1.0d)
            value = major;
        return value;
    }
    
    protected double calculateExpectedOutputState(List<Integer> states,
                                                  List<FactorEdgeType> edgeTypes) {
        int totalWeight = 0;
        int total = 0;
        for (int i = 0; i < edgeTypes.size(); i++) {
            FactorEdgeType edgeType = edgeTypes.get(i);
            int state = states.get(i);
            if (edgeType == FactorEdgeType.OUTPUT) {
                continue;
            }
            if (edgeType == FactorEdgeType.INHIBITOR)
                state = PathwayFGConstants.NUMBER_OF_STATES - 1 - state;
            // Different pre-assigned weights for different edge types
            int weight = FactorEdgeType.getTypeWeight(edgeType);
            total += state * weight;
            totalWeight += weight;
            
            // Geometric mean using the same weight
//            total *= (state + 1);
//            totalWeight ++;
            
        }
        // Arithematic mean
        double expectedState = (double) total / totalWeight;
        // Geometric mean
//        double expectedState = Math.pow(total, 1.0d / totalWeight) - 1.0d;
        return expectedState;
    }
    
    /**
     * Calculate a factor value for one setting, which is defined by the states.
     * @param states
     * @param edgeTypes
     * @return
     */
    private double calculateFactorValue(List<Integer> states,
                                        List<FactorEdgeType> edgeTypes) {
        Integer mergedState = mergeStateForSameLabel(states, edgeTypes);
        // For our reaction based factor, we have only one inhibitor edge
        Integer inhibitorState = getState(states, edgeTypes, FactorEdgeType.INHIBITOR);
        // Calculate expected output state
        Integer expectedOutput = null;
        if (inhibitorState == null) {
            expectedOutput = mergedState;
        }
        else {
            int stateFromInhibitor = PathwayFGConstants.NUMBER_OF_STATES - 1 - inhibitorState;
            expectedOutput = Math.min(mergedState, stateFromInhibitor);
        }
        Integer outputState = getState(states, edgeTypes, FactorEdgeType.OUTPUT);
        return expectedOutput == outputState ? major : minor;
    }
    
    protected Integer getState(List<Integer> states,
                               List<FactorEdgeType> edgeTypes,
                               FactorEdgeType requiredType) {
        Integer state = null;
        for (int i = 0; i < edgeTypes.size(); i++) {
            FactorEdgeType type = edgeTypes.get(i);
            if (type == requiredType) {
                state = states.get(i);
                break;
            }
        }
        return state;
    }
    
    /**
     * Since we have only one edge label, we can use this simple implementation.
     * @param states
     * @param edgeTypes
     * @return
     */
    private Integer mergeStateForSameLabel(List<Integer> states,
                                           List<FactorEdgeType> edgeTypes) {
        Integer currentState = null;
        for (int i = 0; i < edgeTypes.size(); i++) {
            Integer state = states.get(i);
            FactorEdgeType type = edgeTypes.get(i);
            FactorEdgeLabel label = FactorEdgeType.mapTypeToLabel(type);
            if (label == null)
                continue;
            if (currentState == null)
                currentState = state;
            if (label == FactorEdgeLabel.MIN) {
                if (currentState > state)
                    currentState = state;
            }
            else if (label == FactorEdgeLabel.MAX) {
                if (currentState < state)
                    currentState = state;
            }
        }
        if (currentState == null)
            currentState = 1; // Normal state
        return currentState;
    }
    
    @Test
    public void testTypicalReactionFactorValues() {
        List<Variable> variables = new ArrayList<Variable>();
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        // Check reaction factor
        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.CATALYST, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.ACTIVATOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.INHIBITOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        
        checkFactorValuesForDebug(variables, edgeTypes);
    }
    
    /**
     * This method is used to analyze the factor values. The logic confirmed should be
     * implemented in method generateFacrorValues(List<Variable>, List<FactorEdgeType>).
     * @param variables
     * @param edgeTypes
     */
    private void checkFactorValuesForDebug(List<Variable> variables,
                                           List<FactorEdgeType> edgeTypes) {
        // The following output is for debugging
        StringBuilder builder = new StringBuilder();
        for (Variable var : variables)
            builder.append(var.getName()).append("\t");
        builder.append("Expected\tDiff\tProb(beta=1)\tProb(beta=5)\tProb(normal)");
        System.out.println(builder.toString());
        // How many values
        Set<Double> totalValues = new HashSet<Double>();
        Set<Double> totalDiffs = new HashSet<Double>();
        NormalDistribution normalDist = new NormalDistribution();
        Set<Double> normalValues = new HashSet<Double>();
        
        List<Integer> states = new ArrayList<Integer>();
        for (Variable var : variables)
            states.add(0); // State with state == 0
        boolean isChanged = true;
        int count = 0;
        while (isChanged) { // Need to go through each combination
            // Calculate expected state of the output
            double expectedState = calculateExpectedOutputState(states, edgeTypes);
            int actualOutputState = getState(states, 
                                             edgeTypes,
                                             FactorEdgeType.OUTPUT);
            double diff = Math.abs(expectedState - actualOutputState);
            // Use an exponential function
            double value1 = Math.exp(-diff * 1);
            double value2 = Math.exp(-diff * 5);
            // Use two-sides p-value
            double value3 = (1.0 - normalDist.cumulativeProbability(diff)) * 2;
            normalValues.add(value3);
            
            totalValues.add(value2);
            totalDiffs.add(diff);
            
            System.out.println(StringUtils.join("\t", states) + "\t" + 
                               expectedState + "\t" + 
                               diff + "\t" + 
                               value1 + "\t" +
                               value2 + "\t" + 
                               value3);
            
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
        
        // Check types of values and diffs
        List<Double> valuesList = new ArrayList<Double>(totalValues);
        Collections.sort(valuesList);
        List<Double> diffsList = new ArrayList<Double>(totalDiffs);
        Collections.sort(diffsList);
        List<Double> normalValueList = new ArrayList<Double>(normalValues);
        Collections.sort(normalValueList);
        
        System.out.println("\nDiff\tValue(beta=5)\tValue(normal)");
        for (int i = 0; i < diffsList.size(); i++) {
            int valueIndex = valuesList.size() - i - 1;
            System.out.println(diffsList.get(i) + "\t" + 
                               valuesList.get(valueIndex) + "\t" + 
                               normalValueList.get(valueIndex));
        }
    }
    
    @Test
    public void testGenerateFactorValues() {
        DEBUG = true;
        List<Variable> variables = new ArrayList<Variable>();
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        // Check reaction factor
        System.out.println("Test Reaction:");
        testGenerateVariable(FactorEdgeType.INPUT, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.CATALYST, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.ACTIVATOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.INHIBITOR, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        List<Double> values = generateFactorValues(variables, edgeTypes);
        System.out.println("Values: " + values.size());
        // Check complex
        variables.clear();
        edgeTypes.clear();
        System.out.println("\nTest Complex:");
        testGenerateVariable(FactorEdgeType.COMPLEX, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.COMPLEX, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.COMPLEX, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        values = generateFactorValues(variables, edgeTypes);
        System.out.println("Values: " + values.size());
        // Check EntitySet
        System.out.println("\nTest EntitySet:");
        variables.clear();
        edgeTypes.clear();
        testGenerateVariable(FactorEdgeType.MEMBER, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.MEMBER, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.MEMBER, variables, edgeTypes);
        testGenerateVariable(FactorEdgeType.OUTPUT, variables, edgeTypes);
        values = generateFactorValues(variables, edgeTypes);
        System.out.println("Values: " + values.size());
    }
    
    protected void testGenerateVariable(FactorEdgeType type, 
                                      List<Variable> variables,
                                      List<FactorEdgeType> edgeTypes) {
        Variable variable = new Variable(PathwayFGConstants.NUMBER_OF_STATES);
        variable.setName(type.toString());
        variables.add(variable);
        edgeTypes.add(type);
    }
}
