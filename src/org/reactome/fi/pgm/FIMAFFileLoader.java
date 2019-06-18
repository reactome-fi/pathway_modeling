/*
 * Created on Jun 13, 2016
 *
 */
package org.reactome.fi.pgm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.factorgraph.common.DataType;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.UniProtProteinLengthHelper;

/**
 * A customized MATFileLoader to control the mutation file loading.
 * @author gwu
 *
 */
public class FIMAFFileLoader extends MATFileLoader {
    
    /**
     * Default constructor.
     */
    public FIMAFFileLoader() {
    }
    
    /**
     * In this implementation, we use only missense mutation and average all mutation
     * scores to mitigate the biases from long proteins.
     */
    @Override
    public Map<String, Map<String, Float>> loadSampleToGeneToFIScore(String fileName,
                                                                     String scoreColName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String tokens[] = line.split("\t");
        int maScoreIndex = getMAColumIndex(scoreColName, tokens);
        if (maScoreIndex < 0)
            throw new IllegalArgumentException("There is no MA_FI.score column in the maf file: " + fileName);
        // Get needed index
        Integer variantClassifierIndex = null;
        Integer sampleIndex = null;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(Tumor_Sample_Barcode))
                sampleIndex = i;
            else if (tokens[i].equals(Variant_Classification))
                variantClassifierIndex = i;
        }
        // This map will be returned
        Map<String, Map<String, Float>> sampleToGeneToScore = new HashMap<String, Map<String, Float>>();
        Map<String, Map<String, Integer>> sampleToGeneToMutations = new HashMap<String, Map<String, Integer>>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // If a sample is unknown, exclude this line
            if (tokens[sampleIndex].equals("Unknown"))
                continue;
            // Only load allowed types
            if (!tokens[variantClassifierIndex].equals("Missense_Mutation"))
                continue;
            // Just want to get the code for sample
            String sample = parseSample(tokens[sampleIndex]);
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            Map<String, Integer> geneToMutations = sampleToGeneToMutations.get(sample);
            if (geneToScore == null) {
                geneToScore = new HashMap<String, Float>();
                sampleToGeneToScore.put(sample, geneToScore);
                geneToMutations = new HashMap<String, Integer>();
                sampleToGeneToMutations.put(sample, geneToMutations);
            }
            // Get the score
            Float score = null;
            if (tokens.length - 1 >= maScoreIndex && tokens[maScoreIndex].length() > 0) {
                String token = tokens[maScoreIndex];
                try {
                    score = Float.parseFloat(token);
                }
                catch(NumberFormatException e) {
                    // Just ignore this case
                }
            }
            if (score == null) { // Just in case. This should not occur 
                continue;
            }
            // Store the parsed results
            Float oldScore = geneToScore.get(tokens[0]);
            if (oldScore == null) {
                geneToScore.put(tokens[0], score);
                geneToMutations.put(tokens[0], 1);
            }
//            else if (score > oldScore)
//                geneToScore.put(tokens[0], score);
            else {
                geneToScore.put(tokens[0], oldScore + score);
                int mutations = geneToMutations.get(tokens[0]);
                geneToMutations.put(tokens[0], ++mutations);
            }
        }
        fu.close();
        // Average
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            Map<String, Integer> geneToMutations = sampleToGeneToMutations.get(sample);
            for (String gene : geneToScore.keySet()) {
                Float score = geneToScore.get(gene);
                Integer mutations = geneToMutations.get(gene);
                Float mean = score / mutations;
                geneToScore.put(gene, mean);
            }
        }
        return sampleToGeneToScore;
    }
    
    @Test
    public void checkCorrelationBetweenScoreAndProteinLength() throws Exception {
        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
        config.getTypeToEvidenceFile().put(DataType.Mutation,
                "/Users/gwu/datasets/ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.normalized.simple.maf");
        
        String maFile = config.getTypeToEvidenceFile().get(DataType.Mutation);
        String fiScoreCol = "Log_NormalizedScore";
        setSampleNameLength(null);
        Map<String, Map<String, Float>> sampleToGeneToScore = loadSampleToGeneToFIScore(maFile,
                                                                                        fiScoreCol);
        SummaryStatistics stat = new SummaryStatistics();
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            for (String gene : geneToScore.keySet()) {
                Float score = geneToScore.get(gene);
                stat.addValue(score.doubleValue());
            }
        }
        
        Set<String> fis = config.getFIs();
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        
        System.out.println("Total samples: " + sampleToGeneToScore.size());
        Map<String, Double> geneToScoreSum = new HashMap<String, Double>();
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            for (String gene : fiGenes) {
                Float score = geneToScore.get(gene);
                if (score == null)
                    score = new Float(stat.getMin());
                Double scoreSum = geneToScoreSum.get(gene);
                if (scoreSum == null)
                    geneToScoreSum.put(gene, score.doubleValue());
                else
                    geneToScoreSum.put(gene, score + scoreSum);
            }
        }
        Map<String, Integer> geneToLength = new UniProtProteinLengthHelper().loadGeneToProteinLength();
        List<Double> scoreList = new ArrayList<Double>();
        List<Double> lengthList = new ArrayList<Double>();
        for (String gene : geneToScoreSum.keySet()) {
            Double scoreSum = geneToScoreSum.get(gene);
            Integer length = geneToLength.get(gene);
            if (length == null)
                continue;
            scoreList.add(scoreSum);
            lengthList.add(length.doubleValue());
//            System.out.println(gene + "\t" + scoreSum + "\t" + length);
        }
        PearsonsCorrelation cor = MathUtilities.constructPearsonCorrelation(scoreList, lengthList);
        System.out.println("Pearson correlation: " + cor.getCorrelationMatrix().getEntry(0, 1));
        System.out.println("pvalue: " + cor.getCorrelationPValues().getEntry(0, 1));
    }
    
    @Test
    public void testloadSampleToGeneToFIScore() throws IOException {
        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
        config.getTypeToEvidenceFile().put(DataType.Mutation,
                "/Users/gwu/datasets/ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.normalized.simple.maf");

        String maFile = config.getTypeToEvidenceFile().get(DataType.Mutation);
        String fiScoreCol = config.getProperties().get("FIScoreColumnName");
        setSampleNameLength(null);
        Map<String, Map<String, Float>> sampleToGeneToScore = loadSampleToGeneToFIScore(maFile,
                                                                                        fiScoreCol);
        for (String sample : sampleToGeneToScore.keySet()) {
            // There are 4 mutations in TTN in this sample.
            if (!sample.equals("96dc785c-8417-4813-8d15-c32b22d78b74"))
                continue;
            System.out.println(sample);
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            for (String gene : geneToScore.keySet()) {
                System.out.println(gene + "\t" + geneToScore.get(gene));
            }
            break;
        }
    }
    
}
