/*
 * Created on Jul 9, 2015
 *
 */
package org.reactome.factorgraph.common;

import org.reactome.factorgraph.ContinuousVariable;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;

/**
 * This abstract class generalizes some common methods that may be used for continuous variables and factors.
 * @author gwu
 *
 */
public abstract class ContinuousObservationFactorHandler extends ObservationFactorHandler {
    
    @Override
    public VariableAssignment<Number> parseValue(Double value,
                                                 DataType dataType,
                                                 Variable obsVar) {
        if (!(obsVar instanceof ContinuousVariable))
            throw new IllegalArgumentException(obsVar.getName() + " is not a ContinuousVariable object!");
        VariableAssignment<Number> varAssgn = new VariableAssignment<Number>();
        varAssgn.setVariable(obsVar);
        varAssgn.setAssignment(value);
        return varAssgn;
    }
    
    @Override
    public Variable createObservationVar(String geneName, 
                                         DataType dataType,
                                         VariableManager variableManager) {
        String varName = geneName + "_" + dataType;
        ContinuousVariable variable = new ContinuousVariable();
        variable.setName(varName);;
        variableManager.register(variable);
        return variable;
    }
    
}
