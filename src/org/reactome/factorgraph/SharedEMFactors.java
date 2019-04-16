/*
 * Created on Jun 12, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.List;

/**
 * This class is used to describe a set of EMFactors that should have the same set of parameters. 
 * In other words, this is used for performing shared parameter learning.
 * @author gwu
 *
 */
public class SharedEMFactors extends EMFactor {
    // A list of EMFactors sharing parameters.
    private List<EMFactor> sharedFactors;
    
    /**
     * Default constructor.
     */
    public SharedEMFactors() {
    }

    public List<EMFactor> getSharedFactors() {
        return sharedFactors;
    }

    public void setSharedFactors(List<EMFactor> sharedFactors) {
        validateFactors(sharedFactors);
        this.sharedFactors = sharedFactors;
    }
    
    public void addSharedFactor(EMFactor factor) {
        if (sharedFactors == null)
            sharedFactors = new ArrayList<EMFactor>();
        if (sharedFactors.size() == 0) {
            sharedFactors.add(factor);
            // There is no need to validate
            return;
        }
        EMFactor factor0 = sharedFactors.get(0);
        validateFactors(factor, factor0);
        // If there is no exception thrown in the above method, it is safe to
        // add the passed factor into the list.
        sharedFactors.add(factor); 
    }
    
    /**
     * Validate if factors in the shared list are compatible for learning.
     * All factors in the list should have the same cardinalities in both
     * parent and child nodes. Parent and child should have the same indices
     * in the variable list.
     * @param sharedFactors
     * @return
     */
    private void validateFactors(List<EMFactor> sharedFactors) {
        // Pair-wise comparison
        for (int i = 0; i < sharedFactors.size() - 1; i++) {
            EMFactor factor1 = sharedFactors.get(i);
            for (int j = i + 1; j < sharedFactors.size(); j++) {
                EMFactor factor2 = sharedFactors.get(j);
                validateFactors(factor1, factor2);
            }
        }
    }
    
    private void validateFactors(EMFactor factor1, EMFactor factor2) {
        if (factor1.getParent().getStates() != factor2.getParent().getStates())
            throw new IllegalArgumentException("Parent variables have different cardinalities: " + factor1 + ", " + factor2);
        if (factor1.getChild().getStates() != factor1.getChild().getStates())
            throw new IllegalArgumentException("Child variables different cardinalities: " + factor1 + ", " + factor2);
        int parentIndex1 = factor1.getVariables().indexOf(factor1.getParent());
        int parentIndex2 = factor2.getVariables().indexOf(factor2.getParent());
        if (parentIndex1 != parentIndex2)
            throw new IllegalArgumentException("Orders of parent and child variables are not the same: " + factor1 + ", " + factor2);
    }
    
    
    
    @Override
    public void updateCounts() {
        for (EMFactor factor : sharedFactors) {
            double[] belief = factor.getBelief();
            for (int i = 0; i < belief.length; i++)
                counts[i] += belief[i];
        }
    }

    /**
     * Get the values from the first EMFactor object contained by this object.
     */
    @Override
    public double[] getValues() {
        if (sharedFactors == null || sharedFactors.size() == 0)
            return null;
        return sharedFactors.get(0).getValues();
    }
    
    /**
     * Override this method so that the passed values can be set for all EMFactor objects 
     * managed by this SharedEMFactors object. Note: The passed values object will be
     * shared by all EMFactor objects. So changing one in one EMFactor will be propagated
     * into others.
     */
    @Override
    public void setValues(double[] values) {
        if (sharedFactors == null || sharedFactors.size() == 0)
            return;
        for (EMFactor factor : sharedFactors)
            factor.setValues(values);
    }

    @Override
    public void updateFactorValues() {
        // Get the total parent counts so that we can get conditional probabilities
        // Just need to query one EMFactor in order to parent's stride since we have
        // validate all parents should have the same index and therefore stride
        EMFactor factor0 = sharedFactors.get(0);
        Variable parent = factor0.getParent();
        int parentIndex = factor0.getVariables().indexOf(parent);
        int parentStride = factor0.strides.get(parentIndex);
        // Get all counts for calculate conditional probabilities
        double[] parentCounts = new double[parent.getStates()];
        for (int i = 0; i < counts.length; i++) {
            int targetState = (int) Math.floor(i / (double) parentStride) % parent.getStates();
            parentCounts[targetState] += counts[i];
        }
        if (values == null)
            values = new double[factor0.getValues().length];
        for (int i = 0; i < counts.length; i++) {
            int parentState = (int) Math.floor(i / (double) parentStride) % parent.getStates();
            values[i] = counts[i] / parentCounts[parentState];
        }
        
        // Now copy the learned values to each EMFactor
        // Note: all EM factors contained by this object will share the same values object. This
        // may be error-prone.
        // This method may be called only once for each factor. 
        if (factor0.getValues() != values) {
            for (EMFactor factor : sharedFactors) {
                factor.setValues(values);
            }
        }
    }
    
    @Override
    public String toString() {
        return getName() + ": " + super.toString();
    }
    
}
