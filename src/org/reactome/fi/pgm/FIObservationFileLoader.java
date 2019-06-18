/*
 * Created on Oct 16, 2014
 *
 */
package org.reactome.fi.pgm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.ObservationFactorHandler;
import org.reactome.factorgraph.common.ObservationFileLoader;
import org.reactome.factorgraph.common.VariableManager;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * Customize the ObservationFileLoader class to handle mutation file.
 * @author gwu
 *
 */
public class FIObservationFileLoader extends ObservationFileLoader {
    
    /**
     * Default constructor.
     */
    public FIObservationFileLoader() {
    }
    
    @Override
    protected void loadMutationFile(String fileName,
                                    Collection<Factor> factors,
                                    VariableManager varManager) throws IOException {
        //        parseNoMutationDataFile(fileName, 
        //                                DataType.Mutation, 
        //                                factors, 
        //                                varManager);
        Map<String, Map<String, Float>> sampleToGeneToValue = loadMAFFile(fileName);
        for (String sample : sampleToGeneToValue.keySet()) {
            Map<String, Float> geneToValue = sampleToGeneToValue.get(sample);
            for (String gene : geneToValue.keySet()) {
                ensureCentralDogmaNodes(gene, factors, varManager);
                Float value = geneToValue.get(gene);
                addObservation(value, 
                        gene,
                        sample,
                        DataType.Mutation,
                        varManager,
                        factors);
            }
        }
        // Need to finish the loading by some post-processing steps.
        ObservationFactorHandler factorHandler = getObservationFactorHandler(DataType.Mutation);
        factorHandler.finish();
    }
    
    @Override
    protected Map<String, Map<String, Float>> loadMAFFile(String fileName) throws IOException {
        // If there is no gene mutated in the file, we assume it has 0 state
        // Load all genes used in the FI network
        FileUtility fu = new FileUtility();
        Set<String> fis = FIPGMConfiguration.getConfig().getFIs();
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
                
        //        // Just a simple test for loading gene to score file
        //        fu.setInput(fileName);
        //        Map<String, Float> geneToScore = new HashMap<String, Float>();
        //        String line = fu.readLine();
        //        while ((line = fu.readLine()) != null) {
        //            String[] tokens = line.split("\t");
        //            String gene = tokens[0];
        //            if (!allGenes.contains(gene))
        //                continue;
        //            float score = new Float(tokens[13]);
        ////            System.out.println(gene + "\t" + score);
        //            if (score < 3.66E-16f)
        //                score = 3.66E-16f;
        //            geneToScore.put(gene, new Float(-Math.log(score)));
        //        }
        //        fu.close();
        //        // Perform a randomization
        //        geneToScore = MathUtilities.permutate(geneToScore);
        //        Map<String, Map<String, Float>> sampleToGeneToScore = new HashMap<String, Map<String,Float>>();
        //        sampleToGeneToScore.put("ICGC", geneToScore);
        //        return sampleToGeneToScore;
        
//        MATFileLoader mafFileLoader = new MATFileLoader();
        MATFileLoader mafFileLoader = new FIMAFFileLoader();
        String sampleNameLength = FIPGMConfiguration.getConfig().getProperties().get("sampleNameLength");
        if (sampleNameLength != null) {
            Integer length = new Integer(sampleNameLength);
            if (length == -1)
                mafFileLoader.setSampleNameLength(null);
            else
                mafFileLoader.setSampleNameLength(length);
        }
        String fiScoreColName = FIPGMConfiguration.getConfig().getProperties().get("FIScoreColumnName");
        Map<String, Map<String, Float>> sampleToGeneToMAScore = mafFileLoader.loadSampleToGeneToFIScore(fileName, 
                                                                                                        fiScoreColName);
        // We will assign the minimum value for genes that have not been mentioned in the MAF file
        float minScore = Float.MAX_VALUE;
        // Get the minimum value
        for (Map<String, Float> geneToScore : sampleToGeneToMAScore.values()) {
            for (Float score : geneToScore.values())
                if (score < minScore)
                    minScore = score;
        }
        for (String sample : sampleToGeneToMAScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToMAScore.get(sample);
            // Keep FI genes only
            geneToScore.keySet().retainAll(allGenes);
            Set<String> mutatedGenes = new HashSet<String>(geneToScore.keySet()); // Make a copy to avoid changes in this set
            for (String gene : allGenes) {
                if (!mutatedGenes.contains(gene))
                    geneToScore.put(gene, minScore);
            }
        }
        return sampleToGeneToMAScore;
    }
    
    @Test
    public void testLoadMAFFile() throws Exception {
        String fileName = "test_data/tcga_ov/ov.maf.txt";
        Map<String, Map<String, Float>> sampleToGeneToScore = loadMAFFile(fileName);
        List<Float> values = new ArrayList<Float>();
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            for (Float value : geneToScore.values())
                values.add(value);
        }
        System.out.println("Total values: " + values.size());
        Collections.sort(values);
        FileUtility fu = new FileUtility();
        fu.setOutput("tmp/ov.maf.ma.scores.txt");
        for (float value : values)
            fu.printLine(value + "");
        fu.close();
    }
    
}
