/*
 * Created on Mar 30, 2015
 *
 */
package org.reactome.factorgraph;

import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * A customized Factor class containing CLGVariable with mixture Gaussian distribution. These factors should
 * not exist in the current implementation of LBP. In this implementation, only one discrete variable and one
 * CLGVariable are supported.
 * @author gwu
 *
 */
public class CLGFactor extends ContinuousFactor {
    // One state may be corresponding to multiple distributions: e.g. abnormal in 
    // a two-state gene expression should be related to both up and down distributions.
    // So we need a map from discrete states of the wrapped discrete variable to
    // the normal distributions. These requirements have been wrapped in this list of
    // CLGFactorDistribution objects.
    private List<CLGFactorDistribution> distributions;
    
    
    /**
     * Default constructor.
     */
    public CLGFactor() {
    }
    
    /**
     * Set the distributions for the wrapped continuous Variables. The size of the passed map
     * should be the same as the total states of wrapped discrete Variable for marginallizing.
     * Before calling this method, the wrapped discrete Variable should be set first. 
     * @param stateToDistributions
     */
    public void setDistributions(List<CLGFactorDistribution> distributions) {
        if (distributions == null)
            throw new IllegalArgumentException("The passed argument should not be null!");
        if (discreteVariable == null)
            throw new IllegalStateException("Set the discrete Variable object first before assigning the distributions.");
        if (discreteVariable.getStates() != distributions.size()) {
            throw new IllegalArgumentException("The size of the passed map and the number of discrete Variable states are not the same!");
        }
        // Make sure each state in the discrete variable should be specified
        for (int i = 0; i < discreteVariable.getStates(); i++) {
            boolean found = false;
            for (CLGFactorDistribution fd : distributions) {
                if (fd.getState() == i) {
                    found = true;
                    break;
                }
            }
            if (!found)
                throw new IllegalArgumentException("State " + i + " has no distribution defined.");
        }
        this.distributions = distributions;
    }

    @Override
    protected double[] marginalizeForDiscrete() {
        double[] rtn = getDoubleArray(discreteVariable.getStates());
        for (CLGFactorDistribution fd : distributions) {
            rtn[fd.getState()] = fd.sumWeights();
        }
        normalize(rtn, false, false); // In case they are not normalized
        return rtn;
    }

    @Override
    protected double[] _marginalizeForDiscrete(VariableAssignment<? extends Number> varAssign) {
        return marginalizeForDiscrete(varAssign.getAssignment().doubleValue());
    }

    /**
     * Calculate the marginal for the sole discrete Variable.
     * @param target
     * @return
     */
    public double[] marginalizeForDiscrete(Double clgVariableValue) {
        if (clgVariableValue == null) {
            return marginalize(null, discreteVariable);
        }
        else {
            double[] rtn = getDoubleArray(discreteVariable.getStates());
            for (CLGFactorDistribution fd : distributions) {
                rtn[fd.getState()] = fd.getDensity(clgVariableValue);
            }
            normalize(rtn, false, false);
            return rtn;
        }
    }
    
    /**
     * The distribution for one state for a CLGFactor.
     * @author gwu
     *
     */
    public static class CLGFactorDistribution {
        private Integer state;
        // Mixture weights
        private double[] weights;
        private NormalDistribution[] distributions;
        
        public CLGFactorDistribution(Integer state,
                                     double weight,
                                     NormalDistribution distribution) {
            this.state = state;
            weights = new double[]{weight};
            distributions = new NormalDistribution[]{distribution};
        }
        
        public CLGFactorDistribution(Integer state,
                                     double[] weights,
                                     NormalDistribution[] distributions) {
            if (weights.length != distributions.length)
                throw new IllegalArgumentException("The lengths of weights and distributions should be the same.");
            this.state = state;
            this.weights = weights;
            this.distributions = distributions;
        }

        public Integer getState() {
            return state;
        }

        public double[] getWeights() {
            return weights;
        }

        public NormalDistribution[] getDistributions() {
            return distributions;
        }
        
        /**
         * Sum all weights together.
         * @return
         */
        public double sumWeights() {
            double sum = 0.0d;
            for (double weight : weights)
                sum += weight;
            return sum;
        }
        
        public double getDensity(double value) {
            double density = 0.0d;
            for (int i = 0; i < weights.length; i++) {
                density += weights[i] * distributions[i].density(value);
            }
            return density;
        }
        
    }

}
