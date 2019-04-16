/*
 * Created on Jan 11, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.reactome.r3.util.FileUtility;

/**
 * This class group a set of methods that are related to gene expression data analysis.
 * @author wgm
 *
 */
public class CancerGeneExpressionCommon {
    private FileUtility fu = new FileUtility();
    
    public CancerGeneExpressionCommon() {
    }
    
    public Map<String, Map<String, Double>> loadGeneExp(String fileName) throws IOException {
        return loadGeneExp(fileName, true);
    }
    
    /**
     * Load gene expression data per gene. The value is stored as gene to sample to value.
     * @return
     * @throws IOException
     */
    private Map<String, Map<String, Double>> loadGeneExp(String fileName,
                                                        boolean escapeRowContainNA) throws IOException {
        Map<String, Map<String, Double>> geneToData = new HashMap<String, Map<String,Double>>();
        fu.setInput(fileName);
        // Sample list
        String line = fu.readLine();
        List<String> sampleList = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (String token : tokens) {
            String sample = token.replaceAll("\"", "");
            //if (sample.length() > 12)
            //    sample = sample.substring(0, 12);
            sampleList.add(sample);
        }
        // Starting parsing
        while ((line = fu.readLine()) != null) {
//            if (line.contains("NA"))
//                continue; // Don't want any genes containing "NA".
//
            int index = line.indexOf("\t");
            if (line.substring(index + 1).contains("NA") && escapeRowContainNA)
                continue; // Don't want any genes with values containing "NA".
//            System.out.println(line);
            tokens = line.split("\t");
            // The first one is gene name
            String gene = tokens[0].replace("\"", "");
            Map<String, Double> sampleToValue = new HashMap<String, Double>();
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA") || tokens[i].length() == 0)
                    continue; // Just escape these values
                String sample = sampleList.get(i); // The first sample has been checked.
                Double value = new Double(tokens[i]);
                sampleToValue.put(sample, value);
            }
            if (geneToData.containsKey(gene)) {
                System.out.println("Duplicated gene: " + gene);
            }
            geneToData.put(gene, sampleToValue);
        }
        fu.close();
        // These two genes having NA values in exp
        //geneToData.remove("C1ORF129");
        //geneToData.remove("LRRC50");
        return geneToData;
    }
}
