/*
 * Created on Jul 2, 2015
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.List;

/**
 * A customized Factor class containing one ContinuousVariable and one Variable, which
 * can be used as a leaf factor node and converted into a discrete Factor object based
 * on the value in the ContinuousVariable using some pre-specified distribution.
 * Note: This class is an abstract class. A subclass should provide its own implementation
 * to method marginalizeForDiscrete(VariableAssignment<? extends Number>).
 * @author gwu
 *
 */
public abstract class ContinuousFactor extends Factor {
    // The wrapped ContinuousVariable object
    protected ContinuousVariable continuousVariable;
    // The wrapped discrete Variable object.
    protected Variable discreteVariable;

    public Variable getDiscreteVariable() {
        return discreteVariable;
    }

    public void setDiscreteVariable(Variable discreteVariable) {
        this.discreteVariable = discreteVariable;
    }

    public ContinuousVariable getContinuousVariable() {
        return continuousVariable;
    }

    public void setContinuousVariable(ContinuousVariable variable) {
        this.continuousVariable = variable;
    }

    @Override
    protected double[] sendMessage(FGNode target, InferenceType inferenceType, boolean logSpace) {
        throw new IllegalStateException("This method is not supported in class " + getClass().getName());
    }

    @Override
    protected void updateBelief(boolean logSpace) {
        throw new IllegalStateException("This method is not supported in class CLGFactor " + getClass().getName());
    }

    @Override
    public List<Variable> getVariables() {
        List<Variable> list = new ArrayList<Variable>();
        list.add(discreteVariable);
        list.add(continuousVariable);
        return list;
    }
    
    @Override
    protected double[] marginalize(double[] factorValues,
                                   Variable target) {
        throw new IllegalStateException("This method is not supported in class CLGFactor " + getClass().getName());
    }
    
    /**
     * Calculate the marginal for the wrapped discrete Variable object without the observed value of the
     * wrapped ContinuousVariable.
     * @return
     */
    protected abstract double[] marginalizeForDiscrete();
    
    /**
     * Calculate the marginal of this ContiounsFactor object based on the assignment of the wrapped
     * ContinuousVariable object.
     * @param target
     * @return
     */
    public double[] marginalizeForDiscrete(VariableAssignment<? extends Number> varAssign) {
        if (varAssign == null)
            return marginalizeForDiscrete();
        if (varAssign.getAssignment() == null)
            throw new IllegalArgumentException("The assignment is empty in the varAssign object!");
        if (varAssign.getVariable() != continuousVariable)
            throw new IllegalArgumentException("The wrapped Variable in the varAssign object "
                    + "should be the ContinuousVariable in this ContinuousFactor object.");
        return _marginalizeForDiscrete(varAssign);
    }

    /**
     * The actual method to marginalize this ContinuousFactor object. All required values related to 
     * this ContinousFactor object should be checked before this method is called. However, required
     * values by a subclass should be checked by the subclass itself.
     * @param varAssign
     * @return
     */
    protected abstract double[] _marginalizeForDiscrete(VariableAssignment<? extends Number> varAssign);
    
}
