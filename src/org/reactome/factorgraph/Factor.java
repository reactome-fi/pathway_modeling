/*
 * Created on Feb 26, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlIDREF;
import javax.xml.bind.annotation.XmlTransient;

/**
 * A simple factor.
 * @author gwu
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class Factor extends FGNode {
    // Variables contained by this Factor.
    @XmlElement(name="variable")
    @XmlIDREF
    protected List<Variable> variables;
    // Stride is defined as in Page 358 of the PGM book by Holler and Friedman
    protected List<Integer> strides;
    // Cached strides in order to assign values
    // Use double array should be faster than using ArrayList<Double>
    @XmlElement(name="values")
    protected double[] values;
    // Kept a set of double array for re-using
    @XmlTransient
    private Map<Integer, double[]> statesToMessage;
    
    /**
     * Default constructor.
     */
    public Factor() {
        statesToMessage = new HashMap<Integer, double[]>();
    }
    
    public Factor(Variable from,
                  Variable to,
                  double[] values) {
        this();
        // The values for the vector should (here a = to, from = b)
        // be this order: a0 b0, a0 b1, a1 b0, a1 b1 (the sums of the first two and the last
        // two should be 1.0 respectively.)
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(to);
        variables.add(from);
        setName(from.getName() + "->" + to.getName());
        setVariables(variables);
        setValues(values);
    }
    
    /**
     * Set a value in the factor function.
     */
    public void setValue(Double value,
                         Map<Variable, Integer> variableToState) {
        validVariableToStateArgument(variableToState);
        int index = getIndexForAssignment(variableToState);
        values[index] = value;
    }
    
    public void setValues(List<Double> values) {
        double[] tmp = new double[values.size()];
        for (int i = 0; i < values.size(); i++)
            tmp[i] = values.get(i);
        setValues(tmp);
    }
    
    public void setValues(double[] values) {
        if (getTotalStates() != values.length)
            throw new IllegalArgumentException("The passed values has a size that is not consisitent "
                    + "with the total states of variables contained by this Factor object.");
        this.values = values;
    }
    
    private int getTotalStates() {
        int total = 1;
        if (variables != null) {
            for (Variable var : variables)
                total *= var.getStates();
        }
        return total;
    }

    /**
     * Get the index for a specific assignment to this Factor object.
     * @param variableToState
     * @return
     */
    public int getIndexForAssignment(Map<Variable, Integer> variableToState) {
        int index = 0;
        for (int i = 0; i < variables.size(); i++) {
            Variable var = variables.get(i);
            Integer state = variableToState.get(var);
            if (state == null) {
                throw new IllegalArgumentException("The passed argument doesn't contain state for variable " + var);
            }
            index += state * strides.get(i);
        }
        return index;
    }

    private void validVariableToStateArgument(Map<Variable, Integer> variableToState) {
        // Make sure variableStates has the same length as variable
        if (variables == null || values == null)
            throw new IllegalStateException("Variables have not been assigned to this Factor object.");
        // The following check can be disabled since the argument variableToState can be a bigger map having more
        // than variables for this Factor as long as it has states for all variables.
//        if (variables.size() != variableToState.size())
//            throw new IllegalArgumentException("The passed argument has a different length to the contained variables.");
    }
    
    /**
     * Get the value for a factor state.
     * @param variableToState
     * @return
     */
    public Double getValue(Map<Variable, Integer> variableToState) {
        validVariableToStateArgument(variableToState);
        int index = getIndexForAssignment(variableToState);
        return values[index];
    }
    
    /**
     * Send message from this factor to a target variable node.
     * @param variable
     * @param values
     * @return
     */
    @Override
    protected double[] sendMessage(FGNode target,
                                   InferenceType inferenceType,
                                   boolean logSpace) {
        // The following two checks have been truned off for performance reason
        // Have to make sure the links are correct in the client code.
//        if (!(target instanceof Variable))
//            throw new IllegalArgumentException("The passed target should be a Variable object.");
//        // Make sure the target variable is contained by this factor
//        if (!variables.contains(target))
//            throw new IllegalArgumentException("The passed target, " + target + ", is not contained by this factor object.");
        if (message == null)
            message = new double[values.length];
        // Reset message
        for (int i = 0; i < values.length; i++)
            message[i] = (logSpace ? Math.log(values[i]) : values[i]);
        multiple(target, message, logSpace);
        // Marginalize variables except target
        if (logSpace)
            convertLogToProb(message);
        double[] rtn = null;
        if (inferenceType == InferenceType.MAX_PRODUCT) {
            rtn = maximize(message, (Variable)target);
        }
        else {// Default inference type should be SUM_PRODUCT
            rtn = marginalize(message, (Variable)target);
        }
        normalize(rtn, false, logSpace);
//        // The following code should not be used in the production environment
//        for (int i = 0; i < rtn.length; i++) {
//            if (Double.isNaN(rtn[i]))
//                throw new IllegalStateException("A Message contains NaN: a possible numerical underflow occurs. Probably the log-space should be used for computation.");
////            if (Double.isInfinite(rtn[i])) // This should be allowed 
////                throw new IllegalStateException("A Message contains Infinity: a possible numerical overflow occurs.");
//        }
        return rtn;
    }
    
    private double[] multiple(FGNode target,
                              double[] factorValues,
                              boolean logSpace) {
        for (Edge edge : getInEdges()) {
            if (edge.getFromNode() == target)
                continue;
            multiple(factorValues,
                     edge.getMessage(),
                     (Variable)edge.getFromNode(),
                     logSpace);
        }
        return factorValues;
    }
    
    @Override
    protected void updateBelief(boolean logSpace) {  
        if (belief == null)
            belief = new double[values.length];
        System.arraycopy(values, 0, belief, 0, values.length);
        if (logSpace) {
            for (int i = 0; i < belief.length; i++)
                belief[i] = Math.log(belief[i]);
        }
        for (Edge edge : getInEdges()) {
            Variable var = (Variable) edge.getFromNode();
            multiple(belief,
                     edge.getMessage(),
                     var,
                     logSpace);
        }
        normalize(belief, logSpace, false);
    }
    
    protected double[] getDoubleArray(Integer states) {
        double[] array = statesToMessage.get(states);
        if (array == null) {
            array = new double[states];
            statesToMessage.put(states, array);
        }
        else {
            // Reset
            for (int i = 0; i < array.length; i++)
                array[i] = 0.0d;
        }
        return array;
    }
    
    /**
     * Marginalize the passed factorValues for the passed target.
     * @param factorValues
     * @param target
     * @return
     */
    protected double[] marginalize(double[] factorValues, 
                                   Variable target) {
        // Get index of the target variable so that we can get its assignment
        int targetIndex = variables.indexOf(target);
        int targetStride = strides.get(targetIndex);
        double[] rtn = getDoubleArray(target.getStates());
        for (int i = 0; i < factorValues.length; i++) {
            int targetState = (int) Math.floor(i / (double) targetStride) % target.getStates();
            rtn[targetState] += factorValues[i];
        }
        return rtn;
    }
    
    /**
     * This method is used to find the maximum message for MAX_PRODUCT inference.
     * Don't mix this method with methods used for EM learning. The implementation 
     * of this method is similar to method marginalize().
     * @param factorValues
     * @param target
     * @return
     */
    private double[] maximize(double[] factorValues,
                              Variable target) {
        // Get index of the target variable so that we can get its assignment
        int targetIndex = variables.indexOf(target);
        int targetStride = strides.get(targetIndex);
        double[] rtn = getDoubleArray(target.getStates());
        for (int i = 0; i < factorValues.length; i++) {
            int targetState = (int) Math.floor(i / (double) targetStride) % target.getStates();
            if (factorValues[i] > rtn[targetState])
                rtn[targetState] = factorValues[i];
        }
        return rtn;
    }
    
    /**
     * Get the assignment for variables contained by this Factor object for the specified
     * index in the values.
     * @param valueIndex
     * @return
     */
    public Map<Variable, Integer> getAssignment(int valueIndex) {
        Map<Variable, Integer> varToAssign = new HashMap<Variable, Integer>();
        for (int i = 0; i < variables.size(); i++) {
            Variable var = variables.get(i);
            int stride = strides.get(i);
            int assign = (int) Math.floor(valueIndex / (double) stride) % var.getStates();
            varToAssign.put(var, assign);
        }
        return varToAssign;
    }
    
    /**
     * The following implementation is based on Algorithm 10.A.1
     * in the PGM book by Koller and Friedman. The implementation is
     * simplified since this is a factor graph and the passed variable
     * should be contained by this Factor object.
     * @param factorValues
     * @param inMessage
     * @param variable
     */
    private void multiple(double[] factorValues,
                          double[] inMessage,
                          Variable variable,
                          boolean logSpace) {
        int j = 0; // index for factorValues
        int k = 0; // index for inMessage
        int[] assignment = new int[variables.size()]; // picked state for each variable in this Factor object.
        for (int i = 0; i < factorValues.length; i++) { // Cycle through each value
            if (logSpace)
                factorValues[i] += inMessage[k];
            else
                factorValues[i] *= inMessage[k]; // Multiple two values from two lists.
            // The following for-loop is used to reset j and k based on each variable's assignment
            for (int v = 0; v < variables.size(); v++) {
                Variable tmpVar = variables.get(v);
                // The stride is defined as 0 if tmpVar is not the same as variable.
                // Since inMessage is for one variable only, its stride is 1 always.
                int inMessageStride = (tmpVar == variable ? 1 : 0);
                assignment[v] ++;
                if (assignment[v] == tmpVar.getStates()) { // Need to reset assignment for tmpVar
                    assignment[v] = 0; // Reset back to the first state
                    j -= (tmpVar.getStates() - 1) * strides.get(v);
                    k -= (tmpVar.getStates() - 1) * inMessageStride; 
                }
                else {
                    j += strides.get(v);
                    k += inMessageStride;
                    break;
                }
            }
        }
    }
    
    public double[] getValues() {
        return this.values;
    }

    public List<Variable> getVariables() {
        return variables;
    }

    public void setVariables(List<Variable> variables) {
        if (variables == null || variables.size() == 0)
            throw new IllegalArgumentException("Variables is empty!");
        this.variables = variables;
        int totalState = getTotalStates();
        this.values = new double[totalState];
        // Strides
        strides = new ArrayList<Integer>(variables.size());
        int stride = 1;
        for (int i = 0; i < variables.size(); i++) {
            strides.add(stride); // The first stride is 1 always
            stride *= variables.get(i).getStates();
        }
    }
    
    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("Variables: " + getVariables() + ", values:");
        if (values != null) {
            for (double value : values)
                builder.append(" ").append(value);
        }
        return builder.toString();
    }

}
