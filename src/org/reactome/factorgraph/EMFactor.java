/*
 * Created on Jun 10, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.List;

/**
 * This class is implemented as a subclass to Factor, which is subject to EM learning.
 * The current EM implementation can be used for conditional probability learning. In 
 * other words, parameters related to two variables, one parent and one child.
 * @author gwu
 *
 */
public class EMFactor extends Factor {
    protected Variable parent;
    protected Variable child;
    // Used to hold temp values for EM
    protected double[] counts;
    
    /**
     * Default constructor.
     */
    public EMFactor() {
    }
    
    public EMFactor(Variable from,
                    Variable to,
                    double[] values) {
        super(from, to, values);
        setParent(from);
        setChild(to);
    }

    public Variable getParent() {
        return parent;
    }
    
    public void setCounts(double[] counts) {
        this.counts = counts;
    }
    
    public double[] getCounts() {
        return this.counts;
    }
    
    public void initCounts() {
        if (counts == null)
            counts = new double[getValues().length];
        // Assigning 1 as initial counts can avoid 0 in cases where no values
        // are available. This is different from the original algorithm described
        // in the PGM book. Using 1.0 is followed the implementation in libdai. Using
        // 1.0 may result a little increase of loglikelihood of the trained model.
        // This is also called "Laplace's correction" (Page 735, the PGM book).
        for (int i = 0; i < counts.length; i++)
            counts[i] = 1.0d;
    }
    
    /**
     * This is the actual expectation step in EM.
     */
    public void updateCounts() {
        for (int i = 0; i < belief.length; i++)
            counts[i] += belief[i];
    }
    
    /**
     * This is the actual maximization step in EM
     */
    public void updateFactorValues() {
        double[] parentCounts = marginalize(counts, parent);
        int parentIndex = variables.indexOf(parent);
        int parentStride = strides.get(parentIndex);
        for (int i = 0; i < counts.length; i++) {
            int parentState = (int) Math.floor(i / (double) parentStride) % parent.getStates();
            values[i] = counts[i] / parentCounts[parentState];
        }
    }
    
    public void randomFactorValues() {
        // Random counts
        if (counts == null)
            counts = new double[getValues().length];
        for (int i = 0; i < counts.length; i++)
            counts[i] = Math.random();
        updateFactorValues();
    }

    public void setParent(Variable parent) {
        if (variables == null || !variables.contains(parent))
            throw new IllegalArgumentException(parent + " is not listed as a variable in this LearningFactor object!");
        this.parent = parent;
    }

    public Variable getChild() {
        return child;
    }

    public void setChild(Variable child) {
        if (variables == null || !variables.contains(child))
            throw new IllegalArgumentException(child + " is not listed as a variable in this LearningFactor object.");
        this.child = child;
    }

    @Override
    public void setVariables(List<Variable> variables) {
        if (variables == null || variables.size() != 2)
            throw new IllegalArgumentException("Only two variables can be assigned to a LearningFactor object.");
        super.setVariables(variables);
    }
    
}
