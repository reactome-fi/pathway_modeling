/*
 * Created on Jun 7, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This file is used to load MAF (mutation annotation format) based on this URL: 
 * https://wiki.nci.nih.gov/display/TCGA/2.4+Sequence-based+Data.
 * @author wgm
 *
 */
public class MATFileLoader {
    // The following is a list of pre-defined column names
    protected static final String Tumor_Sample_Barcode = "Tumor_Sample_Barcode";
    protected static final String Variant_Classification = "Variant_Classification";
    private static final String amino_acid_change = "amino_acid_change";
    private static final String Reference_Allele = "Reference_Allele";
    // Specify the length of sample name for parsing
    private Integer sampleNameLength = 12; // Default used for TCGA barcode. If this value is null, the whole barcode will be used.
    
    public MATFileLoader() {
    }
    
    public void setSampleNameLength(Integer length) {
        this.sampleNameLength = length;
    }
    
    public Integer getSampleNameLength() {
        return this.sampleNameLength;
    }
    
    /**
     * Load the annotated mutation after checking with the passed gene lengths.
     * @param fileName
     * @param geneToLength
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadSampleToGenes(String fileName,
                                                      String aaChangeColumn,
                                                      Integer sampleNameLength,
                                                      Map<String, Integer> geneToLength) throws IOException {
        setSampleNameLength(sampleNameLength);
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        Set<String> allowedTypes = getAllowedMutationTypes();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Integer aaChangeIndex = null;
        Integer variantClassifierIndex = null;
        Integer sampleIndex = null;
        String[] tokens = line.split("\t");
        if (aaChangeColumn == null)
            aaChangeColumn = amino_acid_change; // Use the default
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(Tumor_Sample_Barcode))
                sampleIndex = i;
            if (tokens[i].equals(Variant_Classification))
                variantClassifierIndex = i;
            if (tokens[i].equals(aaChangeColumn)) {
                aaChangeIndex = i;
            }
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (tokens[sampleIndex].equals("Unknown"))
                continue;
            if (!allowedTypes.contains(tokens[variantClassifierIndex].toLowerCase()))
                continue;
            Integer geneLength = null;
            if (geneToLength != null) 
                geneLength = geneToLength.get(tokens[0]);
            if (geneLength != null && aaChangeIndex != null) {
                Integer aaPos = getAAPosition(tokens[aaChangeIndex]);
                if (aaPos != null && aaPos > geneLength)
                    continue; // The annotation is not the same as we have. Ignore this mutation.
            }
            // Just want to get the code for sample
            String sample = parseSample(tokens[sampleIndex]);
            InteractionUtilities.addElementToSet(sampleToGenes, 
                                                 sample,
                                                 tokens[0]);
        }
        fu.close();
        return sampleToGenes;
    }
    
    protected String parseSample(String token) {
        if (sampleNameLength == null)
            return token;
        else
            return token.substring(0, sampleNameLength);
    }
    
    /**
     * Load the map from sample to gene to maximum mutation score.
     * @param fileName
     * @param maScoreColName
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, Float>> loadSampleToGeneToFIScore(String fileName,
                                                                     String maScoreColName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String tokens[] = line.split("\t");
        int maScoreIndex = getMAColumIndex(maScoreColName, tokens);
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
        // Variant types to be checked with
        Set<String> allowedTypes = getAllowedMutationTypes();
        Set<String> maxScoreTypes = getMutationTypesForMaxMAFI();
        // This map will be returned
        Map<String, Map<String, Float>> sampleToGeneToMAScore = new HashMap<String, Map<String, Float>>();
        Float maxScore = Float.NEGATIVE_INFINITY;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // If a sample is unknown, exclude this line
            if (tokens[sampleIndex].equals("Unknown"))
                continue;
            // Only load allowed types
            if (!allowedTypes.contains(tokens[variantClassifierIndex].toLowerCase()))
                continue;
            // Just want to get the code for sample
            String sample = parseSample(tokens[sampleIndex]);
            Map<String, Float> geneToScore = sampleToGeneToMAScore.get(sample);
            if (geneToScore == null) {
                geneToScore = new HashMap<String, Float>();
                sampleToGeneToMAScore.put(sample, geneToScore);
            }
            // Get the score
            Float score = null;
            if (maxScoreTypes.contains(tokens[variantClassifierIndex].toLowerCase())) {
                score = Float.MAX_VALUE;
            }
            else if (tokens.length - 1 >= maScoreIndex && tokens[maScoreIndex].length() > 0) {
                String token = tokens[maScoreIndex];
                try {
                    score = Float.parseFloat(token);
                    if (score > maxScore) // Want to find the real maximum score
                        maxScore = score;
                }
                catch(NumberFormatException e) {
                    // Just ignore this case
                }
            }
            if (score == null) { // Just in case. This should not occur 
//                System.err.println(new IllegalArgumentException("Error in line: " + line));
                continue;
            }
            // Store the parsed results
            Float oldScore = geneToScore.get(tokens[0]);
            if (oldScore == null || oldScore < score)
                geneToScore.put(tokens[0], score);
        }
        fu.close();
        // Replace Float.MAX_Score with the maximum score
        for (String sample : sampleToGeneToMAScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToMAScore.get(sample);
            Set<String> maxGenes = new HashSet<String>();
            for (String gene : geneToScore.keySet()) {
                Float score = geneToScore.get(gene);
                if (score.equals(Float.MAX_VALUE))
                    maxGenes.add(gene);
            }
            for (String maxGene : maxGenes)
                geneToScore.put(maxGene, maxScore);
        }
        return sampleToGeneToMAScore;
    }
    
    /**
     * This method is used to load sample to genes.
     * @param fileName
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadSampleToGenes(String fileName,
                                                      boolean useHomo) throws IOException {
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        Set<String> allowedTypes = getAllowedMutationTypes();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Integer variantClassifierIndex = null;
        Integer sampleIndex = null;
        Integer referenceAlleleIndex = null;
        String[] tokens = line.split("\t");
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(Tumor_Sample_Barcode))
                sampleIndex = i;
            else if (tokens[i].equals(Variant_Classification))
                variantClassifierIndex = i;
            else if (tokens[i].equals(Reference_Allele))
                referenceAlleleIndex = i;
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // If a sample is unknown, exclude this line
            if (tokens[sampleIndex].equals("Unknown"))
                continue;
            if (!allowedTypes.contains(tokens[variantClassifierIndex].toLowerCase()))
                continue;
            if (useHomo) {
                String ref = tokens[referenceAlleleIndex];
                String allel1 = tokens[referenceAlleleIndex + 1];
                String allel2 = tokens[referenceAlleleIndex + 2];
                if (ref.equals(allel1) ||
                    ref.equals(allel2))
                    continue;
            }
            // Just want to get the code for sample
            String sample = parseSample(tokens[sampleIndex]);
            InteractionUtilities.addElementToSet(sampleToGenes, 
                                                 sample,
                                                 tokens[0]);
        }
        fu.close();
        return sampleToGenes;
    }
    
    /**
     * Load the mutation frequency at each position for each gene from a MAF file.
     * @param fileName
     * @param geneToProteinLength
     * @return
     * @throws IOException
     */
    public Map<String, int[]> loadGeneToSiteMutationFrequency(String fileName,
                                                              String aaChangeColumn,
                                                              Set<String> excludedSamples,
                                                              Map<String, Integer> geneToProteinLength) throws IOException {
        Map<String, int[]> geneToSiteFrequency = new HashMap<String, int[]>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // Headers
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int variantIndex = 0, aaChangeIndex = 0;
        if (aaChangeColumn == null)
            aaChangeColumn = "amino_acid_change";
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("Variant_Classification"))
                variantIndex = i;
            else if (tokens[i].startsWith(aaChangeColumn)) // This most likely works
                aaChangeIndex = i;
        }
        if (variantIndex == 0)
            throw new IllegalArgumentException("Cannot find column for Variant_Classification");
        if (aaChangeIndex == 0)
            throw new IllegalArgumentException("Cannot find column for " + aaChangeColumn);
        Set<String> allowedTypes = getAllowedMutationTypes();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String sample = tokens[15].substring(0, 12);
            if (excludedSamples != null && excludedSamples.contains(sample))
                continue;
            // We will check missense_mutation only
            if (!allowedTypes.contains(tokens[variantIndex].toLowerCase()))
                continue;
            // amino_acid_change_WU should have a format like this: p.G184A
            // We want to extract the position part here
            String aaChange = tokens[aaChangeIndex];
            if (!aaChange.startsWith("p."))
                continue; // Want to check with aa change only
            String gene = tokens[0];
            int[] siteFrequency = geneToSiteFrequency.get(gene);
            if (siteFrequency == null) {
                Integer length = geneToProteinLength.get(gene);
                if (length == null) {
//                    System.err.println(new IllegalArgumentException(gene + " has no length!"));
                    continue;
                }
                siteFrequency = new int[length];
                geneToSiteFrequency.put(gene, siteFrequency);
            }

            Integer position = getAAPosition(aaChange);
            if (position == null)
                continue;
            if (position > siteFrequency.length) {
//                System.err.println(new IllegalArgumentException("The aa change position is too big in line: " +
//                        aaChange + " in gene " + gene + " with length " + siteFrequency.length));
                continue;
            }
            if (position == 0)
                position = 1; // Maybe frameshift at the start of the chain. Use it as one.
            siteFrequency[position - 1] ++; // Assuming the first position is 1.
        }
        fu.close();
        return geneToSiteFrequency;
    }
    
    private Integer getAAPosition(String token) {
        if (!token.startsWith("p."))
            return null;
        Pattern pattern = Pattern.compile("\\d+");
        Matcher matcher = pattern.matcher(token);
        if (matcher.find()) {
            String group = matcher.group();
            return new Integer(group);
        }
        return null;
    }
    
    /**
     * Load gene to MA_FI.score, which is generated by using MSKCC mutation assessor.
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadGeneToMAFIScore(String fileName,
                                                   String fiImpactScoreColumnName,
                                                   Set<String> excludedSamples) throws IOException {
        Map<String, Map<String, Float>> sampleToGeneToScore = loadSampleToGeneToFIScore(fileName,
                                                                                        fiImpactScoreColumnName);
        if (excludedSamples != null && excludedSamples.size() > 0)
            sampleToGeneToScore.keySet().removeAll(excludedSamples);
        Map<String, List<Double>> geneToScores = new HashMap<String, List<Double>>();
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            for (String gene : geneToScore.keySet()) {
                Float score = geneToScore.get(gene);
                List<Double> scores = geneToScores.get(gene);
                if (scores == null) {
                    scores = new ArrayList<Double>();
                    geneToScores.put(gene, scores);
                }
                scores.add(score.doubleValue());
            }
        }
        // Get the maximum score for each gene
        // Use a stat object for easy to switch to other stats
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (String gene : geneToScores.keySet()) {
            List<Double> scores = geneToScores.get(gene);
            stat.clear();
            for (Double score : scores) {
                stat.addValue(score);
            }
            geneToScore.put(gene, stat.getMax());
//            geneToScore.put(gene, stat.getMean());
//            geneToScore.put(gene, stat.getPercentile(50.0d));
        }
        return geneToScore;
    }

    protected int getMAColumIndex(String columnName, String[] tokens) {
        // Find the index of MA_FI.score
        int scoreIndex = -1;
        if (columnName == null)
            columnName = "MA_FI.score";
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(columnName)) {
                scoreIndex = i;
                break;
            }
        }
        return scoreIndex;
    }
    
    public Set<String> getAllowedMutationTypes() {
        String[] arrays = new String[] {"Frame_Shift_Del", "In_Frame_Del", 
                "In_Frame_Ins", "Frame_Shift_Ins", "Nonsense_Mutation", 
                "Nonstop_Mutation", "Missense_Mutation", "Splice_Site_Ins", 
                "Splice_Site_SNP", "Splice_site_SNP", "Splice_Site_Del", "StopCodon_DNP", 
                "Splice_Site", "Stop_Codon_Del", "Init_Met_Del"};
        Set<String> set = new HashSet<String>();
        for (String type : arrays)
            set.add(type.toLowerCase());
        return set;
//        return new HashSet<String>(Arrays.asList(arrays));
    }
    
    private Set<String> getMutationTypesForMaxMAFI() {
        String[] arrays = new String[] {"Frame_Shift_Del", "In_Frame_Del", 
                "In_Frame_Ins", "Frame_Shift_Ins", "Nonsense_Mutation", 
                "Nonstop_Mutation", "Splice_Site_Ins", "Splice_Site_SNP",
                "Splice_Site_Del", "StopCodon_DNP", 
                "Stop_Codon_Del"};
        Set<String> set = new HashSet<String>();
        for (String type : arrays)
            set.add(type.toLowerCase());
        return set;
    }
}
