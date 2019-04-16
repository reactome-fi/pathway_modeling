/*
 * Created on Oct 15, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;

/**
 * This class is used to help to handle a set of variables during converting.
 * @author gwu
 *
 */
public class VariableSetHandler {
    // Other helper objects
    private PathwayVariableManager variableManager;
    private FactorValueAssigner valueAssigner;
    
    /**
     * Default constructor.
     */
    public VariableSetHandler() {
    }
    
    public void setVariableManager(PathwayVariableManager variableManager) {
        this.variableManager = variableManager;
    }

    public void setValueAssigner(FactorValueAssigner valueAssigner) {
        this.valueAssigner = valueAssigner;
    }

    /**
     * Use this method to handle a set of variables, e.g., a set of variables for complex's subunits
     * or entityset's members.
     * @param variables
     * @param targetVar
     * @param edgeType
     * @param factors
     */
    public void handleSetOfVariables(List<Variable> variables,
                                     Variable targetVar,
                                     FactorEdgeType edgeType,
                                     Set<Factor> factors) {
        // Create a factor with this size of variables
        if (variables.size() <= PathwayFGConstants.MAXIMUM_AUGUMENT_NODE) {
            Factor factor = createFactor(variables,
                                         targetVar,
                                         edgeType);
            factors.add(factor);
            return;
        }
        // Want to make sure the results are always the same, so a soring is needed.
        Collections.sort(variables, new Comparator<Variable>() {
            public int compare(Variable var1, Variable var2) {
                if (var1.getName() == null || var2.getName() == null)
                    return 0; // Do nothing
                return var1.getName().compareTo(var2.getName());
            }
        });
        // To be used for next method call in the recursive calling
        List<Variable> nextList = new ArrayList<Variable>();
        // For variable naming
        for (int i = 0; i < variables.size(); i += PathwayFGConstants.MAXIMUM_AUGUMENT_NODE) {
            // If the left variables is less than the maximum node, stop here
            if (i + PathwayFGConstants.MAXIMUM_AUGUMENT_NODE > variables.size()) {
                nextList.addAll(variables.subList(i, variables.size()));
                break;
            }
            List<Variable> subList = variables.subList(i, i + PathwayFGConstants.MAXIMUM_AUGUMENT_NODE);
            // Create an accessory variable node
            // We will not use cached variable. So always call this method to create a new variable
            Variable accessoryVar = variableManager.createVar(targetVar.getName() + "_accessory",
                                                              PathwayPGMConfiguration.getConfig().getNumberOfStates());
            Factor factor = createFactor(subList, 
                                         accessoryVar,
                                         edgeType);
            factors.add(factor);
            nextList.add(accessoryVar);
        }
        handleSetOfVariables(nextList,
                             targetVar, 
                             edgeType,
                             factors);
    }
    
    private Factor createFactor(List<Variable> variables, 
                                Variable targetVar,
                                FactorEdgeType edgeType) {
        // Copy the original list to be used for generate values
        List<Variable> varList = new ArrayList<Variable>(variables);
        List<FactorEdgeType> typeList = new ArrayList<FactorEdgeType>();
        for (int i = 0; i < varList.size(); i++)
            typeList.add(edgeType);
        // Add targetVar
        varList.add(targetVar);
        typeList.add(FactorEdgeType.OUTPUT);
        Factor factor = variableManager.createFactor(varList, typeList, valueAssigner, targetVar.getName());
        // Copy DB_ID
        String dbId = targetVar.getProperty(ReactomeJavaConstants.DB_ID);
        factor.setProperty(ReactomeJavaConstants.DB_ID, dbId);
        return factor;
    }
    
}
