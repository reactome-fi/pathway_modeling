/*
 * Created on Jul 23, 2014
 *
 */
package org.reactome.fi.pgm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.cancer.CancerGeneExpressionCommon;

/**
 * This class is used to process results generated from the PGM_FI model.
 * @author gwu
 *
 */
public class FIPGMResultAnalyzer {
    
    /**
     * Default construcotr.
     */
    public FIPGMResultAnalyzer() {
    }
    
    @Test
    public void analyzeInferenceResults() throws IOException {
        String fileName = FIPGMConfiguration.RESULT_DIR + "PGM_FI_Inference_Results_072214.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(fileName);
        Set<String> samples = new HashSet<String>();
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            samples.addAll(sampleToValue.keySet());
        }
        List<String> sampleList = new ArrayList<String>(samples);
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            int count = 0;
            for (String gene : geneToSampleToValue.keySet()) {
                Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                Double value = sampleToValue.get(sample);
                if (value >= 0.5)
                    count ++;
            }
            System.out.println(sample + "\t" + count);
        }
    }
    
}
