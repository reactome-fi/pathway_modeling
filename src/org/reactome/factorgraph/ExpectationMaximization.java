/*
 * Created on Jun 10, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * This class implements EM for conditional probability. The implementation has
 * made sure the total probabilities have not changed to avoid the change of the
 * partition function (aka Z), which will breaks the calculation of likelihood
 * calculation. 
 * This class is a direct implementation of Algorithm 19.2 in the PGM book by
 * Koller and Friedman.
 * @author gwu
 *
 */
public class ExpectationMaximization {
    private static final Logger logger = Logger.getLogger(ExpectationMaximization.class);
    private List<Observation<Integer>> evidences;
    // The inference algorithm to be used. Currently
    // only LBP is supported
    private Inferencer inferencer;
    // total iteration
    private int maxIteration = 50; // Default
    private double tolerance = 1.0e-3; // log likelihood tolerance between two iterations. This is a relative change.
    // To control if debug is needed: the debug is controled by log4j.
    // So in order to see debugging information, log4j should be enabled
    private boolean debug;
    
    /**
     * Default constructor.
     */
    public ExpectationMaximization() {
        // Default to use LBP
        inferencer = new LoopyBeliefPropagation();
    }
    
    public void setDebug(boolean debug) {
        this.debug = debug;
    }
    
    public boolean getDebug() {
        return this.debug;
    }
    
    public void setInferenser(Inferencer inferencer) {
        this.inferencer = inferencer;
    }
    
    public Inferencer getInferencer() {
        return this.inferencer;
    }

    public List<Observation<Integer>> getEvidences() {
        return evidences;
    }

    public void setEvidences(List<Observation<Integer>> evidences) {
        this.evidences = evidences;
    }
    
    public int getMaxIteration() {
        return maxIteration;
    }

    public void setMaxIteration(int maxIteration) {
        this.maxIteration = maxIteration;
    }

    public double getTolerance() {
        return tolerance;
    }

    public void setTolerance(double tolerance) {
        this.tolerance = tolerance;
    }

    /**
     * Call this method to learn parameters for the passed list of factors.
     * @param factorGraph
     * @param factors parameters in those factors should be learned
     * @return the loglikelihood value at the end of the learning
     */
    public double learn(FactorGraph factorGraph,
                        List<EMFactor> factors) throws InferenceCannotConvergeException {
        // Need a target to learn
        if (factors == null || factors.size() == 0)
            throw new IllegalArgumentException("No learning factors has been provided.");
        // Have to make sure evidences have been assigned
        if (evidences == null || evidences.size() == 0)
            throw new IllegalStateException("No evidences have been assigned for the EM learning.");
        inferencer.setFactorGraph(factorGraph);
        int iteration = 0;
        // Here the change of logZ is used to measure the log likelihood change since
        // log likelihood can be calculated as logZ - logZ0, here logZ is for a factor 
        // graph having evidence assigned, logZ0 the original factor graph.
        double preLogLikelihood = 0.0d;
        double logLikelihoodDiff = Double.MAX_VALUE;
        while (iteration < maxIteration && logLikelihoodDiff > tolerance) {
            double logLikelihood = expect(factorGraph, factors);
            maximize(factorGraph, factors);
            // At least to run twice in order to compare likelihood
            if (iteration > 0) {
                logLikelihoodDiff = logLikelihood - preLogLikelihood;
                if (logLikelihoodDiff < 0)
                    logger.warn("During learning loglikelihood decrease by: " + logLikelihoodDiff);
                logLikelihoodDiff = Math.abs(logLikelihoodDiff / preLogLikelihood);
            }
            preLogLikelihood = logLikelihood;
            iteration ++;
            if (debug) {
                logger.info("Loglikelihood for iteration " + iteration + ": " + logLikelihood + 
                            " and difference: " + (iteration > 1 ? logLikelihoodDiff : "N/A"));
                logger.info("Learned parameters:");
                for (EMFactor factor : factors)
                    logger.info(factor.getName() + ": " + Arrays.toString(factor.getValues()));
            }
        }
        return preLogLikelihood;
    }
    
    /**
     * The implementation of Procedure Computer-ESS in Algorithm 19.2 in the PGM book.
     * @param fg
     * @param factors
     */
    private double expect(FactorGraph fg, List<EMFactor> factors) throws InferenceCannotConvergeException {
        for (EMFactor factor : factors)
            factor.initCounts();
        // Since we want to get a difference between two runs of EM learning,
        // we don't need to actually calculate the background logZ in order to 
        // calculate loglikelihood to avoid the problem of multiple maximum: two
        // runs may get different values.
        // However, keeping the background logZ can increase the convergence speed.
        // So for the time being, this is kept.
        inferencer.clearObservation(); // Have to reset observation
        inferencer.runInference();
        double logZ = inferencer.calculateLogZ();
        if (debug)
            logger.info("LogZ: " + logZ);
        double loglikelihood = 0.0d;
        for (Observation<Integer> obs : evidences) {
//            logger.info("Sample " + obs.getName());
            Map<Variable, Integer> varToAssignment = obs.getVariableToAssignment();
            inferencer.setObservation(varToAssignment);
            inferencer.runInference();
            for (EMFactor factor : factors) {
                expect(factor);
            }
//            loglikelihood += lbp.calculateLogZ(fg);
            loglikelihood += (inferencer.calculateLogZ() - logZ);
        }
        return loglikelihood;
    }
    
    protected void expect(EMFactor factor) {
        factor.updateCounts();
    }
    
    private void maximize(FactorGraph fg, List<EMFactor> factors) {
        for (EMFactor factor : factors)
            factor.updateFactorValues();
    }
}
