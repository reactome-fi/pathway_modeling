/*
 * Created on Jul 2, 2015
 *
 */
package org.reactome.factorgraph;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.reactome.factorgraph.ContinuousVariable.DistributionType;

/**
 * This Factor object is used to describe a relationship between a ContinuousVariable, whose distribution is
 * modeled by EmpiricalDistribution, and a discrete Variable.
 * @author gwu
 *
 */
public class EmpiricalFactor extends ContinuousFactor {
    
    /**
     * Default constructor.
     */
    public EmpiricalFactor() {
    }
    
    @Override
    protected double[] marginalizeForDiscrete() {
        double[] rtn = getDoubleArray(discreteVariable.getStates());
        // Assign the same values for all elements
        for (int i = 0; i < rtn.length; i++) {
            rtn[i] = 1.0d;
        }
        return rtn;
    }
    
    /* (non-Javadoc)
     * @see org.reactome.factorgraph.ContinuousFactor#marginalizeForDiscrete(org.reactome.factorgraph.VariableAssignment)
     */
    @Override
    protected double[] _marginalizeForDiscrete(VariableAssignment<? extends Number> varAssign) {
        if (varAssign.getDistribution() == null)
            throw new IllegalArgumentException("EmpiricalDistribution is not set in the varAssign object!");
        double[] rtn = null;
        // The following should be safe since it should be checked by the supclass method.
        ContinuousVariable cVar = (ContinuousVariable) varAssign.getVariable();
        EmpiricalDistribution distribution = varAssign.getDistribution();
        if (cVar.getDistributionType() == DistributionType.TWO_SIDED) {
            if (discreteVariable.getStates() == 2) {
                // There will be two states: normal and abnormal
                double cpd = distribution.cumulativeProbability(varAssign.getAssignment().doubleValue());
                rtn = getDoubleArray(2);
                if (cpd > 0.50d) {
                    cpd = 1.0d - cpd; // Switch to the other direction in the distribution
                }
                rtn[0] = cpd * 2.0d; // 0 should be normal
                rtn[1] = 1.0d - rtn[0]; // 1 should be abnormal
            }
            else if (discreteVariable.getStates() == 3) {
                // There will be three states: lower, normal, and higher
                double cpd = distribution.cumulativeProbability(varAssign.getAssignment().doubleValue());
                rtn = getDoubleArray(3);
                if (cpd < 0.50d) { // Should be more likely in the lower state
                    rtn[0] = 1.0d - 2.0d * cpd; // Lower
                    rtn[1] = 2.0d * cpd; // Normal
                    rtn[2] = 0.0d; // Assume there is zero change in the higher state by thinking the mean as a scale.
                }
                else if (cpd > 0.50d) { // Should be more likely in the upper state
                    cpd = 1.0d - cpd;
                    rtn[0] = 0.0d;
                    rtn[1] = 2.0d * cpd;
                    rtn[2] = 1.0d - 2.0d * cpd;
                }
                else { // In the normal state by thinking the scale is on balance
                    rtn[0] = 0.0d;
                    rtn[1] = 2.0d * cpd;
                    rtn[2] = 0.0d;
                }
            }
        }
        else { // ONE_SIDED should be used as the default
            // Only two states are supported for ONE_SIDED distribution
            if (discreteVariable.getStates() == 2) {
                // There will be two states: normal and abnormal
                double cpd = distribution.cumulativeProbability(varAssign.getAssignment().doubleValue());
                rtn = getDoubleArray(2);
                rtn[0] = 1.0d - cpd;
                rtn[1] = cpd;
            }
        }
        if (rtn != null)
            return rtn;
        throw new IllegalStateException("Cannot marginalize a ContinuousFactor: " +
                cVar.getDistributionType() + " ContinuousVariable and " +
                discreteVariable.getStates() + " states discrete Variable.");
    }
    
}
