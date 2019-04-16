/*
 * Created on Aug 30, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;

/**
 * This class is used to assign values to factors by iterating all combinations of states. It needs to be 
 * set a Factor object first. The client code should call method iterate() to get all state combinations.
 * @author gwu
 *
 */
public class FactorValueAssignment {
    private Factor factor;
    
    /**
     * Default constructor.
     */
    public FactorValueAssignment() {
    }

    public Factor getFactor() {
        return factor;
    }

    public void setFactor(Factor factor) {
        this.factor = factor;
    }
    
    /**
     * Get the wrapped factor value specified by an index.
     * @param index
     * @return
     */
    public double getFactorValue(int index) {
        return factor.getValues()[index];
    }

    /**
     * The client should call this method to iterate all combinations of states for the assigned
     * factor.
     */
    public List<Map<Variable, Integer>> iterate() {
        // Make sure both maxState and factor have been assigned.
        if (factor == null)
            throw new IllegalStateException("Assign a factor for this FactorValueAssignment object.");
        // Keep track states of variables in the factor. 
        List<Integer> states = new ArrayList<Integer>();
        List<Variable> variables = factor.getVariables();
        for (Variable var : variables) {
            states.add(0);
        }
        // Start iterating
        List<Map<Variable, Integer>> assignments = new ArrayList<Map<Variable,Integer>>();
        // The first assignment has not been covered by the following while loop.
        Map<Variable, Integer> assignment = createAssignment(states, variables);
        assignments.add(assignment);
        while (!isDone(states)) {
            for (int i = 0; i < states.size(); i++) {
                Integer state = states.get(i);
                Variable var = variables.get(i);
                if (state < var.getStates() - 1) {
                    states.set(i, ++state);
                    for (int j = i - 1; j >= 0; j--) {
                        states.set(j, 0); // Just reset
                    }
                    break;
                }
            }
            assignment = createAssignment(states, variables);
            assignments.add(assignment);
        }
        return assignments;
    }

    /**
     * Create an assignment from the states List based on the original Variable list.
     * @param states
     * @param variables
     * @return
     */
    private Map<Variable, Integer> createAssignment(List<Integer> states,
                                                    List<Variable> variables) {
        // Assign states for return
        Map<Variable, Integer> assignment = new HashMap<Variable, Integer>();
        for (int i = 0; i < variables.size(); i++) {
            Variable var = variables.get(i);
            Integer state = states.get(i);
            assignment.put(var, state);
        }
        return assignment;
    }
    
    /**
     * Check if all states have been iterated through. In other words, make sure all variables
     * have got their maximum states. This is a lazy implementation. There is a a quicker way to
     * check this.
     * @param states
     * @return
     */
    private boolean isDone(List<Integer> states) {
        List<Variable> variables = factor.getVariables();
        for (int i = 0; i < states.size(); i++) {
            Variable var = variables.get(i);
            Integer state = states.get(i);
            if (state < var.getStates() - 1) // State starts from 0.
                return false; // If the current state is still less than the total state, 
                              // the iteration is not done yet.
        }
        return true;
    }
    
    /**
     * Test the iterate method.
     */
    @Test
    public void testIterate() {
        Factor factor = new Factor();
        List<Variable> variables = new ArrayList<Variable>();
        for (int i = 0; i < 3; i++) {
            Variable var = new Variable();
            var.setName("Var" + i);
            var.setStates(3);
            variables.add(var);
        }
        factor.setVariables(variables);
        setFactor(factor);
        List<Map<Variable, Integer>> assignments = iterate();
        System.out.println("Total assignments: " + assignments.size());
        for (Map<Variable, Integer> assignment : assignments) {
            for (Variable var : variables) {
                System.out.print(var + ": " + assignment.get(var) + "\t");
            }
            System.out.println();
        }
    }
    
}
