/*
 * Created on Jul 7, 2015
 *
 */
package org.reactome.factorgraph.common;

import org.reactome.factorgraph.EMFactor;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;


/**
 * This customized ObservationFactorHandler is used to create a discrete Variable for an observed data item.
 * @author gwu
 *
 */
public class DiscreteObservationFactorhandler extends ObservationFactorHandler {
    
    /**
     * Default constructor.
     */
    public DiscreteObservationFactorhandler() {
    }

    @Override
    public VariableAssignment<Number> parseValue(Double value, 
                                                 DataType dataType,
                                                 Variable variable) {
        double[] threshold = configuration.getThreshold(dataType);
        Integer assignment = getAssignment(value,
                                           threshold,
                                           variable.getStates()); // We use only one value here.
        VariableAssignment<Number> varAssgn = new VariableAssignment<Number>();
        varAssgn.setVariable(variable);
        varAssgn.setAssignment(assignment);
        return varAssgn;
    }
    
    private Integer getAssignment(Double value,
                                  double[] threshold,
                                  int states) {
        Integer assignment = null;
        // Two cases: 
        // If there is only one value
        if (threshold.length == 1) {
            if (value < threshold[0])
                assignment = 0;
            else
                assignment = 1;
            return assignment;
        }
        for (int i = 0; i < threshold.length; i++) {
            if (value < threshold[i]) {
                assignment = i;
                break;
            }
        }
        // The largest state
        if (assignment == null)
            assignment = threshold.length;
        if (states == threshold.length + 1)
            return assignment;
        if (states == 2 && threshold.length == 2) {
            if (assignment == 1)
                return 0;
            return 1;
        }
        throw new IllegalArgumentException("Cannot support the lenght of threshold and the number of states: " + threshold.length + ", " + states);
    }

    @Override
    public Variable createObservationVar(String geneName, 
                                         DataType dataType,
                                         VariableManager variableManager) {
        String varName = geneName + "_" + dataType;
        return variableManager.getVarForName(varName, 
                                             configuration.getNumberOfStates());
    }

    @Override
    public Factor createObservationFactor(Variable anchorVar,
                                          Variable obsVar, 
                                          DataType dataType) {
        double[] factorValues = configuration.getDataTypeValues(dataType);
        // Support EM learning
        // The following order is important: if two nodes has the following causal direction
        // A -> B. B should be added first into the VecVar, and the values for the vector should
        // be this order: a0 b0, a0 b1, a1 b0, a1 b1 (the sums of the first two and the last
        // two should be 1.0 respectively.)
        EMFactor factor = new EMFactor(anchorVar, obsVar, factorValues);
        return factor;
    }
    
}
