/*
 * Created on Jul 8, 2015
 *
 */
package org.reactome.factorgraph.common;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.reactome.factorgraph.CLGFactor;
import org.reactome.factorgraph.CLGFactor.CLGFactorDistribution;
import org.reactome.factorgraph.ContinuousVariable;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;

/**
 * @author gwu
 *
 */
public class CLGObservationFactorHandler extends ContinuousObservationFactorHandler {
    
    /**
     * Default constructor.
     */
    public CLGObservationFactorHandler() {
    }
    
    /* (non-Javadoc)
     * @see org.reactome.factorgraph.common.ObservationFactorHandler#createObservationFactor(org.reactome.factorgraph.Variable, org.reactome.factorgraph.Variable, org.reactome.factorgraph.common.DataType)
     */
    @Override
    public Factor createObservationFactor(Variable anchorVar,
                                          Variable obsVar, 
                                          DataType dataType) {
        if (!(obsVar instanceof ContinuousVariable))
            throw new IllegalArgumentException(obsVar.getName() + " is not a ContinuousVariable object!");
        // Support two or three states discrete variable only
        if (anchorVar.getStates() != 2 && anchorVar.getStates() != 3) {
            throw new IllegalArgumentException("Support two or three states CLG observation only.");
        }
        ContinuousVariable clgTo = (ContinuousVariable) obsVar;
        CLGFactor clgFactor = new CLGFactor();
        clgFactor.setContinuousVariable(clgTo);
        clgFactor.setDiscreteVariable(anchorVar);
        // Need to create a list of Gaussian distributions
        List<NormalDistribution> dists = new ArrayList<NormalDistribution>();
        // These distributions are created based on threshold values
        double[] thresholds = configuration.getThreshold(dataType);
        //TODO: This will be changed. Right now, just use this simple scheme.
        for (int i = 0; i < thresholds.length; i++) {
            double value = thresholds[i];
            NormalDistribution dist = new NormalDistribution(value, 1.0d);
            dists.add(dist);
            // We also need the value for the middle one
            if (i < thresholds.length - 1) {
                value = (thresholds[i] + thresholds[i + 1]) / 2.0d;
                dist = new NormalDistribution(value, 1.0d);
                dists.add(dist);
            }
        }
        //TODO: The assignment of weights should be changed too. Right now we just assign 
        // equal weights.
        double[] weights = new double[dists.size()];
        double value = 1.0d / weights.length;
        for (int i = 0; i < weights.length; i++)
            weights[i] = value;
        List<CLGFactorDistribution> factorDistrbutions = generateFactorDistributions(anchorVar,
                                                                                     dists,
                                                                                     weights);
        if (factorDistrbutions.size() == 0)
            throw new IllegalStateException("This type of CLGFactor is not supported: state " + anchorVar.getStates() + " and threshold length " + thresholds.length);
        clgFactor.setDistributions(factorDistrbutions);
        return clgFactor;
    }

    private List<CLGFactorDistribution> generateFactorDistributions(Variable anchorVar,
                                                                    List<NormalDistribution> dists,
                                                                    double[] weights) {
        List<CLGFactorDistribution> factorDistrbutions = new ArrayList<CLGFactorDistribution>();
        if (anchorVar.getStates() == 2) {
            if (dists.size() == 2) {
                for (int i = 0; i < dists.size(); i++) {
                    CLGFactorDistribution cfd = new CLGFactorDistribution(i, weights[i], dists.get(i));
                    factorDistrbutions.add(cfd);
                }
            }
            else if (dists.size() == 3) {
                // Pick up the middle one as the normal state
                CLGFactorDistribution cfd = new CLGFactorDistribution(0, weights[1], dists.get(1));
                factorDistrbutions.add(cfd);
                cfd = new CLGFactorDistribution(1, 
                                                new double[]{weights[0], weights[2]},
                                                new NormalDistribution[]{dists.get(0), dists.get(2)});
                factorDistrbutions.add(cfd);
            }
        }
        else if (anchorVar.getStates() == 3) {
            if (dists.size() == 3) {
                for (int i = 0; i < dists.size(); i++) {
                    CLGFactorDistribution cfd = new CLGFactorDistribution(i, weights[i], dists.get(i));
                    factorDistrbutions.add(cfd);
                }
            }
        }
        return factorDistrbutions;
    }
    
}
