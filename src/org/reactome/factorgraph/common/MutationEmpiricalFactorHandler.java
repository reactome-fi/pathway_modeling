/*
 * Created on Jul 14, 2015
 *
 */
package org.reactome.factorgraph.common;

import java.io.IOException;
import java.util.List;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.reactome.cancer.base.MAFMutationAssessorAnnotator;

/**
 * This class is used to handle mutation related factor generation.
 * @author gwu
 *
 */
public class MutationEmpiricalFactorHandler extends EmpiricalFactorHandler {
    
    /**
     * Default constructor.
     */
    public MutationEmpiricalFactorHandler() {
    }
    
    /**
     * Override this method so that mutation score files can be loaded for distribution.
     * The mutation scores used for distribution are pre-stored in a file and not based on
     * the actual scores in a MAF file.
     */
    @Override
    public EmpiricalDistribution getDistribution(DataType dataType) {
        EmpiricalDistribution distribution = super.getDistribution(DataType.Mutation);
        if (distribution.isLoaded())
            return distribution;
        // Need to load the maScoreFile
        // Try to get from the configuration file
        String maScoreFile = getConfiguration().getProperties().get("maScoreFile");
        if (maScoreFile == null)
            throw new IllegalStateException("maScoreFile has not been set!");
        try {
            List<Double> scores = new MAFMutationAssessorAnnotator().loadMAScores(maScoreFile);
            double[] scoreArray = new double[scores.size()];
            for (int i = 0; i < scores.size(); i++)
                scoreArray[i] = scores.get(i);
            distribution.load(scoreArray);
            return distribution;
        }
        catch(IOException e) {
            throw new IllegalStateException("getDistribution: " + e.getMessage());
        }
    }
    
//    @Test
//    public void test() {
//        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
//        setConfiguration(config);
//        EmpiricalDistribution dist = getDistribution(DataType.Mutation);
//        System.out.println("Mean: " + dist.getNumericalMean());
//        System.out.println("Bin: " + dist.getBinCount());
//        System.out.println("Variance: " + dist.getNumericalVariance());
//        System.out.println("SuuportLowerBound: " + dist.getSupportLowerBound());
//        System.out.println("UpperBound: " + dist.getSupportUpperBound());
//        double cumProb = dist.cumulativeProbability(-5.41d);
//        System.out.println("Prob for 0.0: " + cumProb);
//    }
    
}
