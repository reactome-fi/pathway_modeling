/*
 * Created on Feb 26, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlIDREF;

/**
 * Variable copied from org.libdai.Var.
 * @author gwu
 *
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class Variable extends FGNode {
    // How many states this variable has. The variable state is counted from
    // 0 to states - 1 (e.g. for a three-state variable, the states are 0, 1, 2).
    @XmlAttribute
    private int states;
    // Usually there is no need to maintain this. However,
    // sometimes it is handy to avoid a repetitive query.
    // Factors contains this Variable object.
    @XmlElement(name="factor")
    @XmlIDREF
    protected List<Factor> factors;

    /**
     * Default constructor.
     */
    public Variable() {
    }
    
    public Variable(int states) {
        this.states = states;
    }

    public int getStates() {
        return states;
    }

    public void setStates(int states) {
        if (this.states != 0)
            throw new IllegalStateException("The state of a Variable object cannot be changed!");
        this.states = states;
    }
    
    public List<Factor> getFactors() {
        return factors;
    }

    public void setFactors(List<Factor> factors) {
        // Check to make sure factors contain this variable
        for (Factor factor : factors) {
            if (!factor.getVariables().contains(this))
                throw new IllegalArgumentException(factor + " doesn't contain the variable: " + this);
        }
        this.factors = factors;
    }
    
    public void addFactor(Factor factor) {
        if (!factor.getVariables().contains(this))
            throw new IllegalArgumentException(factor + " doesn't contain the variable: " + this);
        if (factors == null)
            factors = new ArrayList<Factor>();
        if (factors.contains(factor))
            return;
        factors.add(factor);
    }
    
    public void removeFactor(Factor factor) {
        if (factors == null)
            return;
        factors.remove(factor);
    }

    /**
     * Send message from this Variable to the passed Factor object. The passed InferenceType
     * object doesn't play any role in this method since this is a Variable.
     * @param target
     * @param type 
     * @return
     */
    protected double[] sendMessage(FGNode target,
                                   InferenceType type,
                                   boolean logSpace) {
        // The following checks have been turned off for performance reason.
//        if (!(target instanceof Factor))
//            throw new IllegalArgumentException("The passed target should be a Factor object.");
        // Make sure the passed factor contains this variable
//        Factor targetFactor = (Factor) target;
        // The following check if rather expensive. Turn it off!
//        if (!targetFactor.getVariables().contains(this))
//            throw new IllegalArgumentException("The passed Factor object doesn't contain this variable, " + this + ".");
        // It is much faster to use a double array instead of a ArrayList<Double>
        if (message == null)
            message = new double[states];
        for (int i = 0; i < message.length; i++)
            message[i] = logSpace ? 0.0d : 1.0d;
//        for (Edge edge : inEdges) {
        // Use the following for loop is faster than the above simplified for loop.
        for (int i = 0; i < inEdges.size(); i++) {
            Edge edge = inEdges.get(i);
            if (edge.getFromNode() == target)
                continue; // Escape the target factor
            double[] message1 = edge.getMessage();
            multiple(message, message1, logSpace);
        }
        // We need to normalize these messages so that the logZ calculation is correct in 
        // class LoopyBeiefPropagation.
        normalize(message, logSpace, logSpace);
        // The following code should not be used in the production environment
//        for (int i = 0; i < totalMessage.length; i++) {
//            if (Double.isNaN(totalMessage[i]))
//                throw new IllegalStateException("A Message contains NaN: a possible numerical underflow occurs. Probably the log-space should be used for computation.");
////            if (Double.isInfinite(rtn[i])) // This should be allowed 
////                throw new IllegalStateException("A Message contains Infinity: a possible numerical overflow occurs.");
//        }
        return message;
    }
    
    @Override
    protected void updateBelief(boolean logSpace) {
        if (belief == null)
            belief = new double[states];
        // Initialize 
        for (int i = 0; i < states; i++)
            belief[i] = 1.0d;
        for (Edge edge : getInEdges()) {
            double[] message = edge.getMessage();
            multiple(belief, edge.getMessage(), logSpace);
        }
        normalize(belief, logSpace, false);
    }
    
    /**
     * Do a simple multiplication for paired elements in these two list.
     * @param belief
     * @param message
     */
    private void multiple(double[] belief, 
                          double[] message,
                          boolean logSpace) {
        // Make sure the passed two lists have the same length
        if (belief.length != message.length)
            throw new IllegalArgumentException("The passed two list of Double have different lengths!");
        for (int i = 0; i < belief.length; i++) {
            if (logSpace)
                belief[i] += message[i];
            else
                belief[i] *= message[i];
        }
    }
    
}
