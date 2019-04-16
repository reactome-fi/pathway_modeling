/*
 * Created on Nov 24, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 * Using random restart for EM learning.
 * @author gwu
 *
 */
public class RandomRestartEM extends ExpectationMaximization {
    private final Logger logger = Logger.getLogger(RandomRestartEM.class);
    private int restart = 10;
    
    /**
     * Default constructor.
     */
    public RandomRestartEM() {
    }
    
    public void setRestart(int restart) {
        if (restart <= 0)
            return; // Don't accept these values
        this.restart = restart;
    }
    
    public int getRestart() {
        return restart;
    }

    @Override
    public double learn(FactorGraph factorGraph, List<EMFactor> factors) throws InferenceCannotConvergeException {
        List<LearnResult> results = new ArrayList<LearnResult>();
        for (int i = 0; i < restart; i++) {
            logger.info("Random parameters:");
            for (EMFactor factor : factors) {
                factor.randomFactorValues();
                logger.info(factor.getName() + ": " + Arrays.toString(factor.getValues()));
            }
            double loglikelihood = super.learn(factorGraph, factors);
            logger.info("Parameters learned:");
            for (EMFactor factor : factors) {
                logger.info(factor.getName() + ": " + Arrays.toString(factor.getValues()));
            }
            storeLearnResult(loglikelihood, factors, results);
        }
        // Find the result with the highest loglikelihood value
        double max = Double.NEGATIVE_INFINITY;
        LearnResult maxResult = null;
        for (LearnResult result : results) {
            if (result.logLikelihood > max) {
                maxResult = result;
                max = result.logLikelihood;
            }
        }
        // Copy values
        for (EMFactor factor : factors) {
            double[] values = maxResult.factorToValues.get(factor);
            factor.setValues(values);
        }
        return max;
    }
    
    private void storeLearnResult(double logLikehood,
                                  List<EMFactor> factors,
                                  List<LearnResult> results) {
        Map<EMFactor, double[]> factorToValues = new HashMap<EMFactor, double[]>();
        for (EMFactor factor : factors) {
            double[] values = factor.getValues();
            factorToValues.put(factor, Arrays.copyOf(values, values.length));
        }
        LearnResult result = new LearnResult();
        result.factorToValues = factorToValues;
        result.logLikelihood = logLikehood;
        results.add(result);
    }
    
    private class LearnResult {
        private Map<EMFactor, double[]> factorToValues;
        private double logLikelihood;
    }
    
}
