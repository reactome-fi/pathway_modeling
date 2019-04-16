/*
 * Created on May 20, 2014
 *
 */
package org.reactome.fi.pgm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.ObservationFileLoader;
import org.reactome.factorgraph.common.VariableManager;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to construct a PGM from a FI network. Since some
 * of variables are cached in this class, the object of this class should
 * be used only once.
 * @author gwu
 *
 */
public class FIPGMConstructor {
    private boolean enableDebug = false;
    private VariableManager variableManager;
    // Data files
    private Map<DataType, String> evidenceFiles;
    // Load evidence files
    private ObservationFileLoader observationLoader;
    
    /**
     * Default constructor.
     */
    public FIPGMConstructor() {
        variableManager = new VariableManager();
        observationLoader = new FIObservationFileLoader();
        observationLoader.setPGMConfiguration(FIPGMConfiguration.getConfig());
    }
    
    public boolean isEnableDebug() {
        return enableDebug;
    }

    public void setEnableDebug(boolean enableDebug) {
        this.enableDebug = enableDebug;
    }

    public Map<String, Variable> getNameToVar() {
        return variableManager.getNameToVar();
    }
    
    public ObservationFileLoader getObservationLoader() {
        return this.observationLoader;
    }
    
    public void setObservationFileLoader(ObservationFileLoader loader) {
        this.observationLoader = loader;
    }

    private void debugFg(FactorGraph fg, Map<String, Set<String>> nameToPartners) throws IOException {
        String out = "tmp/FIPGM_info.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(out);
        fu.printLine("Total factors: " + fg.getFactors().size() + ", and " + fg.getVariables().size());
        fu.printLine("Name\tVar_Id\tPartners");
        Map<String, Variable> nameToVar = getNameToVar();
        for (String name : nameToPartners.keySet()) {
            fu.printLine(name + "\t" + 
                        nameToVar.get(name).getName() + "\t" +  
                        nameToPartners.get(name).size());
        }
        fu.close();
    }
    
    @Test
    public void testConstructFactorGraph() throws Exception {
        FactorGraph fg = constructFactorGraph(PGMType.PairwiseMRF);
        System.out.println("Total factors: " + fg.getFactors().size());
        System.out.println("Total variables: " + fg.getVariables().size());
        long totalMemory = Runtime.getRuntime().totalMemory();
        long freeMeory = Runtime.getRuntime().freeMemory();
        System.out.println("Used memory: " + (totalMemory - freeMeory) / (1024 * 1024.0d) + " Mb");
    }
    
    /**
     * Use this method to create a factor graph from the FI network.
     * The FI network should be the latest version and configured by the 
     * configured file.
     * @return
     */
    public FactorGraph constructFactorGraph(PGMType pgmType) throws IOException {
        FIPGMFactorValuesAssigner valuesAssigner = PGMType.getValuesAssigner(pgmType);
        return constructFactorGraph(valuesAssigner);
    }

    /**
     * Construct a FactorGraph by using the passed FIPGMFactorValuesAssigner object. By using this
     * method, the client can provide customized implementation of Markov Random Field.
     * @param valuesAssigner
     * @return
     * @throws IOException
     */
    public FactorGraph constructFactorGraph(FIPGMFactorValuesAssigner valuesAssigner) throws IOException {
        Set<String> fis = FIPGMConfiguration.getConfig().getFIs();
        return constructFactorGraph(valuesAssigner, fis);
    }

    public FactorGraph constructFactorGraph(FIPGMFactorValuesAssigner valuesAssigner,
                                            Set<String> fis) throws IOException {
        Map<String, Set<String>> nameToPartners = InteractionUtilities.generateProteinToPartners(fis);
        
        List<Factor> factors = new ArrayList<Factor>();
        for (String fi : fis) {
            Factor factor = convertFIToFactor(fi, 
                                              nameToPartners,
                                              valuesAssigner);
            factors.add(factor);
        }
        // Check if factors are needed for genes
        if (valuesAssigner.isGeneFactorValueSupported()) {
            for (String gene : nameToPartners.keySet()) {
                Factor factor = convertGeneToFactor(gene,
                                                    valuesAssigner);
                factors.add(factor);
            }
        }
        
        if (evidenceFiles != null && evidenceFiles.size() > 0) {
            for (DataType type : evidenceFiles.keySet()) {
                String fileName = evidenceFiles.get(type);
                observationLoader.processDataFile(fileName,
                                               type, 
                                               factors,
                                               variableManager);
            }
        }
        
        FactorGraph fg = new FactorGraph();
        fg.setFactors(new HashSet<Factor>(factors));
        fg.validatVariables();
        if (enableDebug)
            debugFg(fg, nameToPartners);
        return fg;
    }
    
    private Factor convertGeneToFactor(String gene,
                                       FIPGMFactorValuesAssigner valueAssigner) {
        List<Variable> variables = new ArrayList<Variable>();
        Variable variable = variableManager.getVarForName(gene, 2);
        variables.add(variable);
        double[] values = valueAssigner.getValuesForGeneFactor();
        Factor factor = new Factor();
        // For easy to look by assign the gene name as the factor's name
        factor.setName(gene);
        factor.setVariables(variables);
        factor.setValues(values);
        return factor;
    }
    
    private Factor convertFIToFactor(String fi,
                                     Map<String, Set<String>> nameToPartners,
                                     FIPGMFactorValuesAssigner valuesAssigner) {
        List<Variable> variables = new ArrayList<Variable>();
        int index = fi.indexOf("\t");
        String name1 = fi.substring(0, index);
        Variable var1 = variableManager.getVarForName(name1, 2);
        String name2 = fi.substring(index + 1);
        Variable var2 = variableManager.getVarForName(name2, 2);
        variables.add(var1);
        variables.add(var2);
        double[] values = valuesAssigner.getValuesForFIFactor(fi, nameToPartners);
        
        Factor factor = new Factor();
        factor.setName(var1 + " - " + var2);
        factor.setVariables(variables);
        factor.setValues(values);
        return factor;
    }
    
    public void setEvidenceFiles(Map<DataType, String> typeToFile) {
        this.evidenceFiles = typeToFile;
    }
    
    /**
     * Pre-specified types of PGMs that the FI network can be converted.
     */
    public static enum PGMType {
        PairwiseMRF,
        NearestNeighborGibbs,
        Ising,
        DegreeBasedModel; // Probabilities are are assigned based on connection degree.
        
        public static FIPGMFactorValuesAssigner getValuesAssigner(PGMType type) {
            switch (type) {
                case PairwiseMRF : return new PairwiseMRFValuesAssigner();
                case Ising : return new IsingFactorValueAssigner();
                case DegreeBasedModel : return new DegreeModelValuesAssigner();
                case NearestNeighborGibbs : return new NearestNeighborGibbsValuesAssigner();
            }
            throw new IllegalArgumentException(type + " is not supported!");
        }
    }
    
    private static class DegreeModelValuesAssigner implements FIPGMFactorValuesAssigner {
        
        public DegreeModelValuesAssigner() {
        }

        @Override
        public double[] getValuesForFIFactor(String fi,
                                             Map<String, Set<String>> nameToPartners) {
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index);
            String name2 = fi.substring(index + 1);
            
            // See the document slides (PGMForFINetwork.pptx) for how to generate
            // factor functions.
            List<Double> valueList = new ArrayList<Double>();
            // The values are sorted as following for a FI A - B
            // 0 0; 1 0; 0 1; 1 1 (the first for A and the second for B)
            double gamma1 = 1.0d - FIPGMConfiguration.GAMMA;
            
            // Based on node degrees 
            valueList.add(2.0d * gamma1);
            // name1 == A
            valueList.add(1.0d + FIPGMConfiguration.GAMMA - gamma1 / nameToPartners.get(name1).size());
            // name2 == B
            valueList.add(1.0d + FIPGMConfiguration.GAMMA - gamma1 / nameToPartners.get(name2).size());
            valueList.add(gamma1 * (1.0d / nameToPartners.get(name1).size() + 1.0d / nameToPartners.get(name2).size()));
            double[] values = new double[valueList.size()];
            for (int i = 0; i < valueList.size(); i++)
                values[i] = valueList.get(i);
            return values;
        }

        @Override
        public double[] getValuesForGeneFactor() {
            return null;
        }
        
        @Override
        public boolean isGeneFactorValueSupported() {
            return false;
        }
        
    }
    
    private static class PairwiseMRFValuesAssigner implements FIPGMFactorValuesAssigner {
        
        public PairwiseMRFValuesAssigner() {
        }

        @Override
        public double[] getValuesForFIFactor(String fi,
                                             Map<String, Set<String>> nameToPartners) {
            return new double[]{
//                    1.0d - FIPGMConfiguration.THETA,
//                    FIPGMConfiguration.THETA,
//                    FIPGMConfiguration.THETA,
//                    1.0d - FIPGMConfiguration.THETA
//                    1.05, 1, 1, 1.05
                    // These values are based on Nearest Neighbor Gibbs assuming average degree = 30
//                    1.02d, 1.0d, 1.0d, 1.02d // The ratio defines how much the effects can be propagated.
//                    1.002d, 1.0d, 1.0d, 1.002d
//                    100, 1, 1, 100
                    // Just for test
//                    0.26, 0.24, 0.24, 0.26
                    0.31, 0.19, 0.19, 0.31 // Parameters for FIsInGene_031516_BigComp_NoUBC.txt
                    // Just for test
//                    1.0d, 1.02d, 1.02d, 1.0d
            };
        }

        @Override
        public double[] getValuesForGeneFactor() {
            return new double[]{
//                    1.0d - FIPGMConfiguration.ALPHA, 
//                    FIPGMConfiguration.ALPHA
//                    100, 1
                    // These values are based on Nearest Neighbor Gibbs
//                    1.0d, 0.67d
//                    1.0d, 1.0d
                    // Just for test
//                    1.0d, 1.0d
//                    0.6d, 0.4d
                    0.975, 0.025
            };
        }
        
        @Override
        public boolean isGeneFactorValueSupported() {
            return true;
        }
    }
    
    /**
     * The following implementation implements a nearest neighbor Gibbs measure reported
     * by Chen et al: PLoS Genetics 7(4): e1001353 (2011).
     */
    private static class NearestNeighborGibbsValuesAssigner implements FIPGMFactorValuesAssigner {
        // Three parameters that should be tuned. The default values below
        // are copied from the original paper.
//        private double h = -4.0d;
        private double h = -4.6d; // 0.01 in pair-wise MRF
//        private double tao1 = 0.0d;
//        private double tao0 = 0.0d;
        private double tao1 = 0.03d;
//        private double tao0 = 0.002d;
        private double tao0 = 0.21d; // 10 with average degree 30 in pair-wise MRF
        
        
        public NearestNeighborGibbsValuesAssigner() {
        }
        
        @Override
        public double[] getValuesForFIFactor(String fi,
                                             Map<String, Set<String>> nameToPartners) {
            double[] values = new double[4];
            int index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            int degree1 = nameToPartners.get(gene1).size();
            double weight1 = Math.sqrt(degree1);
            int degree2 = nameToPartners.get(gene2).size();
            double weight2 = Math.sqrt(degree2);
            double weightSum = weight1 + weight2;
            values[0] = Math.exp(tao0 * weightSum);
            values[1] = 1.0d;
            values[2] = 1.0d;
            values[3] = Math.exp(tao0 * weightSum);
            return values;
        }
        
        @Override
        public double[] getValuesForGeneFactor() {
            double[] values = new double[] {
                    1.0d, // If a gene is not a driver
                    Math.exp(h) // If a gene is a driver
//                    1.0d - FIPGMConfiguration.ALPHA,
//                    FIPGMConfiguration.ALPHA
            };
            return values;
        }

        @Override
        public boolean isGeneFactorValueSupported() {
            return true;
        }
    }
    
    private static class IsingFactorValueAssigner implements FIPGMFactorValuesAssigner {
        
        public IsingFactorValueAssigner() {
        }

        @Override
        public double[] getValuesForFIFactor(String fi,
                                             Map<String, Set<String>> nameToPartners) {
            double gamma1 = 1.0d - FIPGMConfiguration.GAMMA;
            double[] values = new double[] {
                    gamma1,
                    FIPGMConfiguration.GAMMA,
                    FIPGMConfiguration.GAMMA,
                    gamma1
            };
            return values;
        }

        @Override
        public double[] getValuesForGeneFactor() {
            return null;
        }

        @Override
        public boolean isGeneFactorValueSupported() {
            return false;
        }
        
    }
}
