/*
 * Created on Nov 24, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.List;

import org.apache.log4j.Logger;

/**
 * This class implements a algorithm called deterministic annealing variant of the EM
 * algorithm reported by Ueda and Nakano:  
 * http://papers.nips.cc/paper/941-deterministic-annealing-variant-of-the-em-algorithm.pdf.
 * This implementation actually cannot work in our setting and should not use it.
 * @author gwu
 */
@Deprecated
public class AnnealingEM extends ExpectationMaximization {
    private final static Logger logger = Logger.getLogger(AnnealingEM.class);
    private double currentBeta = 1.0d;
    private double[] beta;
    
    /**
     * Default constrctuor.
     */
    public AnnealingEM() {
    }
    
    /**
     * Set up a set of beta values for annealing. The last value in the
     * array should be 1.0.
     * @param values
     */
    public void setBeta(double[] values) {
        if (!(new Double(values[values.length - 1]).equals(1.0d))) {
            throw new IllegalArgumentException("The last value in the values argument should be 1.0.");
        }
        this.beta = values;
    }
    
    public double[] getBeta() {
        return this.beta;
    }

    @Override
    public double learn(FactorGraph factorGraph, List<EMFactor> factors) throws InferenceCannotConvergeException {
        // The following hard-coded a protocol of annealing step
        double[] beta = new double[] {
               0.001d, 0.005d, 0.01d, 0.05d, 0.10d, 0.50d, 0.75d, 1.00d 
        };
        double loglikelihood = Double.NEGATIVE_INFINITY;
        for (double beta1 : beta) {
            currentBeta = beta1;
            logger.info("currentBeta: " + currentBeta);
            loglikelihood = super.learn(factorGraph, factors);
        }
        return loglikelihood;
    }

    /**
     * The expect step for each factor.
     */
    @Override
    protected void expect(EMFactor factor) {
        // This is rather a hack way for updating counts.
        double[] counts = factor.getCounts();
        if (factor instanceof SharedEMFactors) {
            SharedEMFactors sharedFactors = (SharedEMFactors) factor;
            for (EMFactor factor1 : sharedFactors.getSharedFactors()) {
                double[] belief = factor1.getBelief();
                for (int i = 0; i < belief.length; i++)
                    counts[i] += Math.pow(belief[i], currentBeta);
            }
        }
        else {
            double[] belief = factor.getBelief();
            for (int i = 0; i < belief.length; i++)
                counts[i] += Math.pow(belief[i], currentBeta);
        }
    }
    
}
