/*
 * Created on May 19, 2016
 *
 */
package org.reactome.factorgraph.common;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.reactome.r3.util.FileUtility;

/**
 * Use this class to create an empirical distribution directly from a MAF file.
 * @author gwu
 *
 */
public class MAFMutationEmpiricalFactorHandler extends EmpiricalFactorHandler {
    private String mafFileName;
    private String mutationScoreColumnName;
    
    /**
     * Default constructor.
     */
    public MAFMutationEmpiricalFactorHandler() {
    }

    public String getMafFileName() {
        return mafFileName;
    }

    public void setMafFileName(String mafFileName) {
        this.mafFileName = mafFileName;
    }

    public String getMutationScoreColumnName() {
        return mutationScoreColumnName;
    }

    public void setMutationScoreColumnName(String mutationScoreColuName) {
        this.mutationScoreColumnName = mutationScoreColuName;
    }

    @Override
    public EmpiricalDistribution getDistribution(DataType dataType) {
        EmpiricalDistribution distribution = super.getDistribution(dataType);
        if (distribution.isLoaded())
            return distribution;
        if (mafFileName == null) {
            // Try to get from the configuration
            mafFileName = getConfiguration().getTypeToEvidenceFile().get(DataType.Mutation);
        }
        if (mafFileName == null || mutationScoreColumnName == null)
            throw new IllegalStateException("Mutation file in the MAF format and/or mutation score "
                    + "column name have not been set!");
        try {
            FileUtility fu = new FileUtility();
            fu.setInput(mafFileName);
            String line = fu.readLine();
            String[] tokens = line.split("\t");
            int scoreColIndex = -1;
            for (int i = 0; i < tokens.length; i++) {
                if (tokens[i].equals(mutationScoreColumnName)) {
                    scoreColIndex = i;
                    break;
                }
            }
            if (scoreColIndex == -1) {
                throw new IllegalStateException("Cannot find a column for " + mutationScoreColumnName + " in the MAF file!");
            }
            // Get the score
            Map<String, Double> keyToScore = new HashMap<String, Double>();
            while ((line = fu.readLine()) != null) {
                tokens = line.split("\t");
                if (tokens.length - 1 < scoreColIndex)
                    continue;
                if (tokens[scoreColIndex].length() == 0)
                    continue;
                // Get the key based on mutation position and type
                // 1: chromosome; 2: start_position; 3: end_position; 6: Reference_allele; 7: alternate_allele
                String key = tokens[1] + ":" + tokens[2] + ":" + tokens[3] + ":" + tokens[6] + ":" + tokens[7];
                if (keyToScore.containsKey(key))
                    continue;
                keyToScore.put(key, new Double(tokens[scoreColIndex]));
            }
            double[] scores = new double[keyToScore.size()];
            int index = 0;
            for (String key : keyToScore.keySet())
                scores[index ++] = keyToScore.get(key);
            distribution.load(scores);
            return distribution;
        }
        catch(IOException e) {
            throw new IllegalStateException("getDistribution: " + e.getMessage());
        }
    }
    
//    @Test
//    public void test() {
//        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
//        config.typeToFileName.put(DataType.Mutation,
//                                  "/Users/gwu/datasets/ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.normalized.simple.maf");
//        setConfiguration(config);
////        setMutationScoreColumnName("MetaLR_RankScore");
//        setMutationScoreColumnName("Log_NormalizedScore");
//        EmpiricalDistribution dist = getDistribution(DataType.Mutation);
//        System.out.println("Mean: " + dist.getNumericalMean());
//        System.out.println("Bin: " + dist.getBinCount());
//        System.out.println("Variance: " + dist.getNumericalVariance());
//        System.out.println("SupportLowerBound: " + dist.getSupportLowerBound());
//        System.out.println("UpperBound: " + dist.getSupportUpperBound());
//        double value = -15.25232219696045d;
//        double cumProb = dist.cumulativeProbability(value);
//        System.out.println("Prob for " + value + ": " + cumProb);
//        value = -19.120713257022654;
//        cumProb = dist.cumulativeProbability(value);
//        System.out.println("Prob for " + value + ": " + cumProb);
//        value = -10.4533472061157;
//        cumProb = dist.cumulativeProbability(value);
//        System.out.println("Prob for " + value + ": " + cumProb);
//    }
    
}
