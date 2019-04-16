/*
 * Created on Oct 22, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * A class that is used to analyze the inference results using the implementation in this
 * package.
 * @author gwu
 *
 */
public class ResultsAnalyzer {
    
    /**
     * Default constructor.
     */
    public ResultsAnalyzer() {
    }
    
//    @Test
//    public void compareTwoFiles() throws IOException {
////        String file1 = "results/paradigm/twoCases/hnsc/112914/Node3/YAP1-_and_WWTR1__TAZ_-stimulated_gene_expression.txt";
////        String file2 = "results/paradigm/twoCases/hnsc/121814/Node4/YAP1-_and_WWTR1__TAZ_-stimulated_gene_expression.txt";
//        
////        String file1 = "results/paradigm/twoCases/hnsc/112914/Node4/Metabolism_of_nucleotides.txt";
////        String file2 = "results/paradigm/twoCases/hnsc/121814/Node10/Metabolism_of_nucleotides.txt";
//        
//        String file1 = "tmp/Metabolism_of_nucleotides_fixed_var_1.txt";
//        String file2 = "tmp/Metabolism_of_nucleotides_fixed_var_2.txt";
//        
//        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
//        Map<String, Map<String, Double>> sampleToNodeToValue1 = helper.loadGeneExp(file1);
//        Map<String, Map<String, Double>> sampleToNodeToValue2 = helper.loadGeneExp(file2);
//        for (String sample : sampleToNodeToValue1.keySet()) {
//            Map<String, Double> nodeToValue1 = sampleToNodeToValue1.get(sample);
//            Map<String, Double> nodeToValue2 = sampleToNodeToValue2.get(sample);
//            for (String node : nodeToValue1.keySet()) {
//                Double value1 = nodeToValue1.get(node);
//                Double value2 = nodeToValue2.get(node);
//                Double diff = Math.abs(value2 - value1);
//                if (diff > 1.0)
//                    System.out.println(sample + "\t" + node + "\t" + diff);
//            }
//        }
//    }
    
    @Test
    public void averageLearnedParameters() throws IOException {
//        String dirName = "results/paradigm/twoCases/brca/102414/";
        // Results for TCGA HNSC
        String dirName = "results/paradigm/twoCases/hnsc/110214/";
        double diffThreshold = 1.0e-5;
        List<File> loggingFiles = getLoggingFiles(dirName);
        FileUtility fu = new FileUtility();
        String type = "mRNA_EXP";
        List<List<Double>> values = parseLearnedParameters(loggingFiles, 
                                                           type,
                                                           fu);
        System.out.println("Parameters for " + type);
        for (List<Double> list : values) {
            System.out.println(list);
        }
        List<Double> average = average(values, diffThreshold);
        System.out.println("\nAverage: " + average);
        
        type = "CNV";
        values = parseLearnedParameters(loggingFiles, 
                                        type,
                                        fu);
        System.out.println("\nParameters for " + type);
        for (List<Double> list : values) {
            System.out.println(list);
        }
        average = average(values, diffThreshold);
        System.out.println("\nAverage: " + average);
    }
    
    private List<Double> average(List<List<Double>> values,
                                 double diffThreshold) {
        List<Double> rtn = new ArrayList<Double>();
        int count = 0;
        for (List<Double> list : values) {
            if (!isGood(list, diffThreshold)) {
                System.out.println("Not good: " + list);
                continue;
            }
            if (rtn.size() == 0) {
                rtn.addAll(list);
                continue;
            }
            for (int j = 0; j < list.size(); j++) {
                Double v = rtn.get(j);
                rtn.set(j, v + list.get(j));
            }
            count ++;
        }
        System.out.println("Count: " + count + " out of " + values.size());
        // Perform a normalization
        double sum = 0.0d;
        for (int i = 0; i < rtn.size(); i++) {
            sum += rtn.get(i);
            if ((i + 1) % 3 == 0) {
                for (int j = i; j >= i - 2 ; j--) {
                    rtn.set(j, rtn.get(j) / sum);
                }
                sum = 0.0d;
            }
        }
        return rtn;
    }
    
    private boolean isGood(List<Double> list, double diffThreshold) {
        // Just the first values of probabilities (0, 3, 6) and make sure they have enough difference
        int[] indices = new int[]{0, 3, 6};
        for (int i = 0; i < indices.length - 1; i++) {
            double v1 = list.get(i);
            for (int j = i + 1; j < indices.length; j++) {
                double v2 = list.get(j);
                if (Math.abs(v2 - v1) < diffThreshold)
                    return false;
            }
        }
        // Make sure the first value is larger than the last
        return list.get(0).compareTo(list.get(6)) > 0;
    }

    private List<List<Double>> parseLearnedParameters(List<File> loggingFiles, 
                                                      String type,
                                                      FileUtility fu) throws IOException {
        List<List<Double>> values = new ArrayList<List<Double>>();
        for (File file : loggingFiles) {
            fu.setInput(file.getAbsolutePath());
            String line = null;
            while ((line = fu.readLine()) != null) {
                if (line.contains("Learned parameters: " + type)) {
                    int index = line.lastIndexOf(":");
                    List<Double> list = parseLearnedParameters(line.substring(index + 1).trim());
                    values.add(list);
                }
            }
            fu.close();
        }
        return values;
    }
    
    private List<Double> parseLearnedParameters(String text) {
        String[] tokens = text.split(", |\\]|\\[");
        List<Double> list = new ArrayList<Double>();
        for (String token : tokens) {
            if (token.length() == 0)
                continue;
            list.add(new Double(token));
        }
        return list;
    }
    
    @Test
    public void checkLoggings() throws IOException {
//        String dirName = "results/paradigm/twoCases/brca/102214/";
//        String dirName = "results/paradigm/twoCases/brca/102414/";
        // Results for TCGA HNSC
        String dirName = "results/paradigm/twoCases/hnsc/110214/";
        // Running results using parameters from 110214
        dirName = "results/paradigm/twoCases/hnsc/110614/";
        // Running results on Nov 10, 2014 by using the log-space
        dirName = "results/paradigm/twoCases/hnsc/111014/";
        // Several bug fixes and performance tuning (e.g. use shared double[])
        dirName = "results/paradigm/twoCases/hnsc/111414/";
        // Use fixed hard-coded parameters without learning
        dirName = "results/paradigm/twoCases/hnsc/111714/";
        // Re-run after finding a bug in the Java code
//        dirName = "results/paradigm/twoCases/hnsc/112614/";
        // After many optimizations: running without throwing any exception
        dirName = "results/paradigm/twoCases/hnsc/112914/";
        // Only gene expression data is used in this run
        dirName = "results/paradigm/twoCases/hnsc/121014/";
        // repeat 111914's analysis
        dirName = "results/paradigm/twoCases/hnsc/121814/";
        // Results using MAX_PROUDUCT LBP
        dirName = "results/paradigm/twoCases/hnsc/121914/";
        // Second run with MAX_PRODUCT LBP
        dirName = "results/paradigm/twoCases/hnsc/121914_1/";
        // sum product with gibbs as the back up
        dirName = "results/paradigm/twoCases/hnsc/020215/";
        List<File> loggingFiles = getLoggingFiles(dirName);
        // Check total pathway processed
        int totalProcessed = 0;
        Map<String, String> pathwayToTime = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        int index1, index2;
        for (File logging : loggingFiles) {
            fu.setInput(logging.getAbsolutePath());
            String line = null;
            while ((line = fu.readLine()) != null) {
                if (line.contains("Converting pathway [Pathway")) {
                    totalProcessed ++;
                    index1 = line.lastIndexOf("[");
                    index2 = line.lastIndexOf("...");
                    String pathway = line.substring(index1, index2);
                    pathwayToTime.put(pathway, null);
                }
                else if (line.contains("Time for")) {
                    index1 = line.lastIndexOf("[");
                    index2 = line.lastIndexOf(":");
                    String pathway = line.substring(index1, index2);
                    String time = line.substring(index2 + 1).trim();
                    if (!pathwayToTime.containsKey(pathway))
                        throw new IllegalStateException("Cannot find starting point for pathway " + pathway + 
                                                        " in file " + logging.getAbsolutePath());
                    pathwayToTime.put(pathway, time);
                }
            }
            fu.close();
        }
        System.out.println("Total processed pathway: " + totalProcessed);
        System.out.println("Size of pathwayToTime: " + pathwayToTime.size());
        System.out.println("\nPathways that have not finished:");
        int total = 0;
        for (String pathway : pathwayToTime.keySet()) {
            String time = pathwayToTime.get(pathway);
            if (time == null) {
                System.out.println(pathway);
                total ++;
            }
        }
        System.out.println("In total: " + total);
        
        System.out.println("\nPathways running time (minutes):");
        final Map<String, Double> pathwayToMins = new HashMap<String, Double>();
        for (String pathway : pathwayToTime.keySet()) {
            String time = pathwayToTime.get(pathway);
            if (time == null)
                continue;
            // Need a little parsing
            int index = time.indexOf(" ");
            pathwayToMins.put(pathway, new Double(time.substring(0, index)) / 60.0d);
        }
        // Do a sorting
        List<String> sortedPathways = new ArrayList<String>(pathwayToMins.keySet());
        Collections.sort(sortedPathways, new Comparator<String>() {
            public int compare(String pathway1, String pathway2) {
                Double min1 = pathwayToMins.get(pathway1);
                Double min2 = pathwayToMins.get(pathway2);
                return min2.compareTo(min1);
            }
        });
        for (String pathway : sortedPathways) {
            System.out.println(pathway + "\t" + pathwayToMins.get(pathway));
        }
        System.out.println("In total: " + sortedPathways.size());
    }

    private List<File> getLoggingFiles(String dirName) {
        File dir = new File(dirName);
        List<File> loggingFiles = new ArrayList<File>();
        for (File file : dir.listFiles()) {
            if (file.getName().startsWith("Node")) {
                File loggingFile = new File(file, "logging.txt");
                if (loggingFile.exists())
                    loggingFiles.add(loggingFile);
            }
        }
        System.out.println("Total logging files: " + loggingFiles.size());
        return loggingFiles;
    }
    
    
    
}
