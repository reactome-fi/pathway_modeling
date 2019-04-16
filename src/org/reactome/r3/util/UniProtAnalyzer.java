/*
 * Created on Jun 20, 2006
 *
 */
package org.reactome.r3.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.persistence.MySQLAdaptor;

public class UniProtAnalyzer {
    //public static final String UNI_DIR_NAME = "D:\\documents\\Stein_lab\\Reactome\\Data\\";
    //private final String UNI_SPROT_FILE_NAME =  UNI_DIR_NAME + "uniprot_sprot_human.dat";
    //private final String UNI_TREMBL_FILE_NAME = UNI_DIR_NAME + "uniprot_trembl_human.dat";
    private final String UNI_SPROT_FILE_NAME = R3Constants.UNIPROT_DIR + "uniprot_sprot_human.dat";
    private final String UNI_TREMBL_FILE_NAME = R3Constants.UNIPROT_DIR + "uniprot_trembl_human.dat";
    private MySQLAdaptor reactomeDba;
    private FileUtility fu = new FileUtility();
    
    /**
     * Load the map from gene to the maximum length of its protein products. A gene may have multiple
     * protein products having different lengths. This method returns the longest length. If you want
     * to get the map from a gene to its all protein lengths, use loadGeneToProteinLengths().
     * @return
     * @throws IOException
     */
    public Map<String, Integer> loadGeneToProteinLength() throws IOException {
        Map<String, List<Integer>> geneToLengths = loadGeneToProteinLengths();
        Map<String, Integer> geneToLength = new HashMap<String, Integer>();
        for (String gene : geneToLengths.keySet()) {
            List<Integer> list = geneToLengths.get(gene);
            if (list.size() == 1)
                geneToLength.put(gene, list.get(0));
            else {
                // Find the longest
                Integer length = Integer.MIN_VALUE;
                for (Integer value : list) {
                    if (value > length)
                        length = value;
                }
                geneToLength.put(gene, length);
            }
        }
        return geneToLength;
    }
    
    private String extractGeneName(String line) {
        int index1;
        int index2;
        String geneName;
        // Get gene name
        index1 = line.indexOf("Name=");
        index2 = line.indexOf(";");
        // As of December, 2014, there are some lines as following:
//                    GN   Name=QRSL1 {ECO:0000255|HAMAP-Rule:MF_03150};
//                    GN   Name=GATB {ECO:0000255|HAMAP-Rule:MF_03147,
        if (index2 < 0) // For case: GN   Name=GATB {ECO:0000255|HAMAP-Rule:MF_03147,
            index2 = line.indexOf(" ", index1);
        geneName = line.substring(index1 + "Name=".length(), index2);
        if (geneName.contains("{")) { // For case: Name=QRSL1 {ECO:0000255|HAMAP-Rule:MF_03150}
            index2 = geneName.indexOf("{");
            geneName = geneName.substring(0, index2).trim();
        }
        return geneName;
    }
    
    public Map<String, List<Integer>> loadGeneToProteinLengths() throws IOException {
        String fileName = UNI_SPROT_FILE_NAME;
        Map<String, List<Integer>> geneToLengths = new HashMap<String, List<Integer>>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        Set<String> genes = new HashSet<String>();
        int index1, index2;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("GN   Name=")) { // Gene name
                String gene = extractGeneName(line);
                genes.add(gene);
                // A hard-coded fix
                if (gene.equals("STK26"))
                    genes.add("MST4");
                // Check Synonyms
                index1 = line.indexOf("Synonyms=");
                if (index1 > 0) {
                    index2 = line.indexOf(";", index1);
                    if (index2 > 0) { // Some synonyms may be broken to more than one line. Just ignore them. This occurs in the trembl file only!
                        String tmp = line.substring(index1 + "Synonyms=".length(), index2);
                        String[] tokens = tmp.split(", ");
                        for (String token : tokens) {
                            genes.add(token);
                        }
                    }
                }
            }
            else if (line.startsWith("SQ   SEQUENCE")) { // Sequence length
                if (genes.size() > 0) {
                    index1 = line.indexOf("AA");
                    String text = line.substring("SQ   SEQUENCE".length() + 1, index1).trim();
                    Integer length = new Integer(text);
                    for (String gene : genes) {
                        List<Integer> list = geneToLengths.get(gene);
                        if (list == null) {
                            list = new ArrayList<Integer>();
                            geneToLengths.put(gene, list);
                        }
                        list.add(length);
                    }
                }
            }
            if (line.startsWith("//"))
                genes.clear();
        }
        return geneToLengths;
    }
    
    
}
