package org.reactome.cancer;
/*
 * Created on May 27, 2015
 *
 */


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to attach mutation assessor scores for mutations in a MAF file.
 * @author gwu
 *
 */
public class MAFMutationAssessorAnnotator {
    private String maFIScoreDirName;
    
//    @Test
//    public void generateInputFordbNSFP() throws IOException {
//        // TCGA Pancancer MAF
//        String dir = "datasets/TCGA/TCGA_PanCancer/somatic_mafs_cleaned_filtered/";
//        String mafFile = dir + "pancan12_cleaned_filtered.maf";
//        String outFile = dir + "pancan12_cleaned_filtered.input.txt";
//        DbNSFPAnalyzer analyzer = new DbNSFPAnalyzer();
//        analyzer.generateInputFordbNSFPFromMAF(mafFile,
//                                               outFile);
//        // The above generated output should be used for search_dbNSFP32a as following:
//        // java -Xmx8G search_dbNSFP32a -i ../../ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.input.txt 
//        //                              -o Barcelona_consensus.filter.genes.dbNSFP.output.txt 
//        //                              -v hg19
//    }
//    
//    @Test
//    public void annotateMAFWithdbNSFP() throws IOException {
//        String dir = "datasets/TCGA/TCGA_PanCancer/somatic_mafs_cleaned_filtered/";
//        String mafFile = dir + "pancan12_cleaned_filtered.maf";
//        String outFile = dir + "pancan12_cleaned_filtered.dbNSFP.maf";
//        String scoreFile = dir + "pancan12_cleaned_filtered.output.txt";
//        DbNSFPAnalyzer analyzer = new DbNSFPAnalyzer();
//        analyzer.annotateMAFFileWithdbNSFP(mafFile,
//                                           outFile,
//                                           scoreFile);
//    }
    
    /**
     * Default constructor.
     */
    public MAFMutationAssessorAnnotator() {
    }
    
    /**
     * Set the directory which contains all pre-calculated MA FI scores, and load
     * these scores into memory.
     * @param dirName
     * @throws IOException
     */
    public void setMAFIScoreDirectory(String dirName) throws IOException {
        this.maFIScoreDirName = dirName;
    }

    private Map<String, Double> loadMAFIScore(File file) throws IOException {
        String line;
        String[] tokens;
//        System.out.println("Loading " + file.getName() + "...");
        FileUtility fu = new FileUtility();
        fu.setInput(file.getAbsolutePath());
        line = fu.readLine(); // Escape the first line, which is header
        // Get the score index
        tokens = line.split("\t");
        int scoreIndex = tokens.length - 1;
        Map<String, Double> keyToValue = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
//                System.out.println(line);
            tokens = line.split("\t");
            Double value = null;
            if (tokens.length > scoreIndex) {
                value = new Double(tokens[scoreIndex]);
                keyToValue.put(tokens[0], value);
            }
        }
        fu.close();
        return keyToValue;
    }
    
    /**
     * Annotate a MAF file using MutationAssessor score.
     * @throws IOException
     */
    public void annotateMAFWithMAFIScore(String srcFile,
                                         String targetFile) throws IOException {
        if (maFIScoreDirName == null)
            throw new IllegalStateException("Specify the directory containing MA FI scores first before calling this method.");
        FileUtility fu = new FileUtility();
        fu.setInput(srcFile);
        fu.setOutput(targetFile);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int chrIndex = 0, startPosIndex = 0, endPosIndex = 0, refIndex = 0, tumorAll1Index = 0, tumorAll2Index = 0;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("Chromosome"))
                chrIndex = i;
            if (tokens[i].equals("Start_position"))
                startPosIndex = i;
            if (tokens[i].equals("End_position"))
                endPosIndex = i;
            if (tokens[i].equals("Reference_Allele"))
                refIndex = i;
            if (tokens[i].equals("Tumor_Seq_Allele1"))
                tumorAll1Index = i;
            if (tokens[i].equals("Tumor_Seq_Allele2"))
                tumorAll2Index = i;
        }
        // Make a copy
        fu.printLine(line + "\tMA_FI.score"); // Follow the name used by firehose.
        while ((line = fu.readLine()) != null)
            fu.printLine(line);
        fu.close();
        // We want to annotate by chromosomes to control the usage of memory
        // Load MA FI score
        File maDir = new File(maFIScoreDirName);
        for (File maFile : maDir.listFiles()) {
            Map<String, Double> keyToValue = loadMAFIScore(maFile);
            String chr = getChromosome(maFile);
            // Need to use the target as the source
            File file = new File(targetFile);
            File tmp = new File(targetFile + ".tmp");
            file.renameTo(tmp);
            fu.setInput(tmp.getAbsolutePath());
            fu.setOutput(targetFile);
            line = fu.readLine();
            fu.printLine(line);
            int count = 0;
            while ((line = fu.readLine()) != null) {
                tokens = line.split("\t");
                if (!tokens[chrIndex].equals(chr)) {
                    fu.printLine(line);
                    continue;
                }
                Double score = getScore(chr,
                                        tokens[startPosIndex],
                                        tokens[endPosIndex],
                                        tokens[refIndex],
                                        tokens[tumorAll1Index],
                                        tokens[tumorAll2Index],
                                        keyToValue);
                if (score == null)
                    fu.printLine(line + "\t");
                else {
                    fu.printLine(line + "\t" + score);
                    count ++;
                }
            }
            System.out.println("Total annotated for " + chr + ": " + count);
            fu.close();
        }
    }
    
    private String getChromosome(File maFile) {
        String fileName = maFile.getName();
        int index1 = fileName.indexOf(".");
        int index2 = fileName.indexOf(".", index1 + 1);
        String chr = fileName.substring(index1 + 1, index2);
        chr = chr.substring(3);
        if (chr.startsWith("0"))
            return chr.substring(1);
        else
            return chr;
    }
    
    private Double getScore(String chromosome,
                            String start,
                            String end,
                            String ref,
                            String tumorAll1,
                            String tumorAll2,
                            Map<String, Double> keyToValue) {
        // Some test
        if (!start.endsWith(end))
            return null; // Not supported
        if (ref.equals(tumorAll1) && !ref.equals(tumorAll2)) { // Most common case
            String key = "hg19," + chromosome + "," + start + "," + ref + "," + tumorAll2;
            return keyToValue.get(key);
        }
        // A homogenous case
        if (tumorAll1.equals(tumorAll2) && !ref.equals(tumorAll1)) {
            String key = "hg19," + chromosome + "," + start + "," + ref + "," + tumorAll2;
            return keyToValue.get(key);
        }
        // All other cases are not supported
        return null;
    }
    
    /**
     * Load pre-generated MA scores to build a distribution.
     * @param scoreFileName
     * @return
     * @throws IOException
     */
    public List<Double> loadMAScores(String scoreFileName) throws IOException {
        List<Double> scores = new ArrayList<Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(scoreFileName);
        String line = null;
        while ((line = fu.readLine()) != null)
            scores.add(new Double(line));
        fu.close();
        return scores;
    }
    
    @Test
    public void testPreloadedMAScores() throws IOException {
        long time1 = System.currentTimeMillis();
        String maFIScoreDir = "/Users/gwu/ProgramFiles/mutation_assessor_2014/MA.hg19";
        String allScoreFile = maFIScoreDir + "/AllScores_10.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(allScoreFile);
        String line = null;
        List<Double> allScores = new ArrayList<Double>();
        while ((line = fu.readLine()) != null) {
            allScores.add(new Double(line));
        }
        fu.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for loading: " + (time2 - time1));
//        if (true)
//            return;
        // Just want to take 1/10 sampling
        RandomDataGenerator randomizer = new RandomDataGenerator();
        Object[] allScores1 = randomizer.nextSample(allScores, allScores.size() / 10);
        String scoreFile = maFIScoreDir + "/AllScores_10.txt";
        fu.setOutput(scoreFile);
        for (Object obj : allScores1)
            fu.printLine(obj + "");
        fu.close();
    }
    
    @Test
    public void testLoadMAScores() throws IOException {
        long time1 = System.currentTimeMillis();
        String maFIScoreDir = "/Users/gwu/ProgramFiles/mutation_assessor_2014/MA.hg19";
        setMAFIScoreDirectory(maFIScoreDir);
        // Load MA FI score
        File maDir = new File(maFIScoreDirName);
        List<Double> allScores = new ArrayList<Double>();
        int total = 0;
        for (File maFile : maDir.listFiles()) {
            Map<String, Double> keyToValue = loadMAFIScore(maFile);
            System.out.println(maFile.getName() + ": " + keyToValue.size());
            allScores.addAll(keyToValue.values());
            total += keyToValue.size();
        }
        System.out.println("Total: " + total);
        System.out.println("Size of allScores: " + allScores.size());
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
        String scoreFile = maFIScoreDir + "/AllScores.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(scoreFile);
        for (Double score : allScores)
            fu.printLine(score + "");
        fu.close();
        time1 = System.currentTimeMillis();
        fu.setInput(scoreFile);
        allScores.clear();
        String line = null;
        while ((line = fu.readLine()) != null)
            allScores.add(new Double(line));
        fu.close();
        time2 = System.currentTimeMillis();
        System.out.println("Time for loading scores back: " + (time2 - time1));
    }
    
    @Test
    public void testAnnotateMAFWithMAFFieScore() throws IOException {
        String dirName = "datasets/ICGC/Santa_Cruz_Pilot/AnnotatedMAFs/";
        String srcFileName = dirName + "Sanger.PCAWG_train2.annotated.snv_mnv.filter.genes.maf";
        String targetFileName = dirName + "Sanger.PCAWG_train2.annotated.snv_mnv.filter.genes.MAFI.maf";
        String maFIScoreDir = "/Users/gwu/ProgramFiles/mutation_assessor_2014/MA.hg19";
        setMAFIScoreDirectory(maFIScoreDir);
        annotateMAFWithMAFIScore(srcFileName, targetFileName);
    }
    
}
