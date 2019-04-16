/*
 * Created on Nov 10, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.random.RandomDataImpl;
import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.factorgraph.FactorValueAssignment;
import org.reactome.factorgraph.GibbsSampling;
import org.reactome.factorgraph.InferenceCannotConvergeException;
import org.reactome.factorgraph.LoopyBeliefPropagation;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.DataType;
import org.reactome.pathway.factorgraph.ReactomePathwayFGRunner.ConvertedFactorGraph;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to generate a random set of samples using Gibbs sampling for performance test and 
 * other related work. The implementation is based on class org.reactome.factorgraph.GibbsSampling, which
 * is used for PGM inference.
 * @author gwu
 *
 */
public class GibbsSampler extends GibbsSampling {
    private static final String DIR_NAME = "results/paradigm/twoCases/GibbsSampling/"; 
    // Used file for generated samples
    public static final String SAMPLE_FILE_NAME = DIR_NAME + "Transcription_Regulation_Pluripoten_Samples_AssignedParams.txt";
    private static Logger logger = Logger.getLogger(GibbsSampler.class);
    // Used to detect the minimum probability
    private final double MINIMUM_PROB = 1.0e-12; // Almost arbitrary
    
    /**
     * Default constructor. 
     */
    public GibbsSampler() {
    }
    
    @Test
    public void learnParametersFromGeneratedSamples() throws Exception {
        FileUtility.initializeLogging();
        // Used to work with a specific pathway
        Long dbId = 2032785L; // YAP1 and WWTR1-stimulated gene expression
        // Another even simplier pathway
        dbId = 452723L; // Transcriptional regulation of pluripotent stem cells
        // This is a simple test pathway
//        dbId = -4L;
        
        GKInstance pathway = getPathway(dbId);
        ReactomePathwayFGRunner runner = new ReactomePathwayFGRunner();
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        ConvertedFactorGraph cfg = runner.convertPathway(pathway, converter);
        logger.info("Converted factors: " + cfg.fg.getFactors().size() + " factors, and " + cfg.fg.getVariables().size() + " variables.");
//        for (Factor factor : cfg.fg.getFactors()) {
//            if (factor instanceof EMFactor)
//                System.out.println("EMFactor: " + factor.getName());
//        }
//        cfg.fg.exportFG(System.out);
//        Map<DataType, SharedEMFactors> typeToSharedFactors = cfg.typeToFactors;
//        for (DataType type : typeToSharedFactors.keySet()) {
//            SharedEMFactors factors = typeToSharedFactors.get(type);
//            logger.info(type + ": " + factors.getSharedFactors().size() + " EMFactors.");
//            // We want to use a random initial parameters for a tet
//            factors.randomFactorValues();
//            logger.info("Random parameters: " + factors.getName() + ": " + Arrays.toString(factors.getValues()));
////            for (EMFactor factor : factors.getSharedFactors())
////                System.out.println(factor);
//        }
//        if (true)
//            return;
        
//        typeToSharedFactors.remove(DataType.CNV);

//        String fileName = "Top_TestSampling_AssignedParam_geneExpRxt.txt";
//        String fileName = "Top_TestSampling_new_geneExpRxt.txt";
//        String fileName = "TestSampling_new_geneExpRxt.txt";
        Collection<Observation<Integer>> samples = loadSamples(cfg.fg, SAMPLE_FILE_NAME);
        samples = MathUtilities.randomSampling(samples,
                                               1000, 
                                               new RandomDataImpl()); // Random pick up 1000 samples
        
        filterSamples(samples, cfg.fg);
        checkSamples(samples);
        cfg.observations = new ArrayList(samples);
        runner.performLearning(pathway, converter, cfg);
        runner.performInference(pathway, converter, cfg);
    }
    
    private void checkSamples(Collection<Observation<Integer>> observations) {
        Map<String, List<Integer>> typeToCounts = new HashMap<String, List<Integer>>();
        for (Observation obs : observations) {
            Map<Variable, Integer> varToState = obs.getVariableToAssignment();
            
        }
    }
    
    /**
     * Perform a test inference to remove those samples that cannot be converged.
     */
    private void filterSamples(Collection<Observation<Integer>> obs, FactorGraph fg) {
        LoopyBeliefPropagation lbp = PathwayPGMConfiguration.getConfig().getLBP();
        lbp.setFactorGraph(fg);
        for (Iterator<Observation<Integer>> it = obs.iterator(); it.hasNext();) {
            Observation sample = it.next();
            lbp.setObservation(sample.getVariableToAssignment());
            try {
                lbp.runInference();
            }
            catch(InferenceCannotConvergeException e) {
                logger.error(sample.getName() + " cannot converge!");
                it.remove();
            }
        }
        logger.info("Total samples aftering filtering: " + obs.size());
    }
    
    /**
     * Load the generated samples.
     * @param fileName
     * @return
     * @throws IOException
     */
    private List<Observation<Integer>> loadSamples(FactorGraph fg,
                                                 String fileName) throws IOException {
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (Variable var : fg.getVariables())
            nameToVar.put(var.getName(), var);
        
        List<Observation<Integer>> samples = new ArrayList<Observation<Integer>>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] names = line.split("\t");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Observation<Integer> sample = new Observation<Integer>();
            Map<Variable, Integer> varToState = new HashMap<Variable, Integer>();
            sample.setVariableToAssignment(varToState);
            sample.setName(tokens[0]);
            for (int i = 2; i < tokens.length; i++) {
                String name = names[i];
                if (name.endsWith("_mRNA") || name.endsWith("_CNV")) {
                    Variable var = nameToVar.get(name);
                    if (var == null)
                        continue;
                    varToState.put(var, new Integer(tokens[i]));
                }
            }
            samples.add(sample);
        }
        fu.close();
        return samples;
    }
    
    @Test
    public void checkLoglikelihood() throws Exception {
//        FileUtility fu = new FileUtility();
//        fu.setInput(SAMPLE_FILE_NAME);
//        String line = fu.readLine();
//        double min = Double.POSITIVE_INFINITY;
//        double max = Double.NEGATIVE_INFINITY;
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            double value = new Double(tokens[1]);
//            if (value > max)
//                max = value;
//            if (value < min)
//                min = value;
//        }
//        fu.close();
//        System.out.println("Min: " + min);
//        System.out.println("Max: " + max);
//        if (true)
//            return;
        
        FileUtility.initializeLogging();
        // Used to work with a specific pathway
        Long dbId = 452723L; // Transcriptional regulation of pluripotent stem cells
        
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        GKInstance pathway = getPathway(dbId);
        FactorGraph fg = converter.convertPathway(pathway);
        logger.info("Converted factors: " + fg.getFactors().size() + " factors, and " + fg.getVariables().size() + " variables.");
        
        setFactorGraph(fg);
        setBurnin(100000);
        List<Observation<Integer>> samples = generateSamples(100000);
        
        // We want to fix these variables
        Set<String> outputNames = getOutputEntityNames();
        Set<Variable> outputVars = new HashSet<Variable>();
        for (Variable var : fg.getVariables()) {
            if (outputNames.contains(var.getName()))
                outputVars.add(var);
        }
        
        double max = Double.NEGATIVE_INFINITY;
        String fileName = DIR_NAME + "Loglikelihood_Test.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        for (Observation obs : samples) {
            Map<Variable, Integer> varToState = obs.getVariableToAssignment();
            double originalValue = fg.getLogLikelihood(varToState);
            for (Variable outputVar : outputVars)
                varToState.put(outputVar, 2);
            double value = fg.getLogLikelihood(varToState);
            if (value > max)
                max = value;
            fu.printLine(obs.getName() + "\t" + originalValue + "\t" + value);
        }
        fu.close();
        System.out.println("Maximum loglikelihood after clamping the output variables to 2: " + max);
    }
    
    /**
     * This method is used to run the sample generation.
     */
    @Test
    public void runGenerateSamples() throws Exception {
        FileUtility.initializeLogging();
        // Used to work with a specific pathway
        Long dbId = 2032785L; // YAP1 and WWTR1-stimulated gene expression
        dbId = 452723L; // Transcriptional regulation of pluripotent stem cells
        
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        GKInstance pathway = getPathway(dbId);
        FactorGraph fg = converter.convertPathway(pathway);
        logger.info("Converted factors: " + fg.getFactors().size() + " factors, and " + fg.getVariables().size() + " variables.");
        setFactorGraph(fg);
        setBurnin(100000);
        List<Observation<Integer>> samples = generateSamples(100000);
        //String fileName = "tmp/TestSampling_AssignedParam_geneExpRxt.txt";
        String fileName = SAMPLE_FILE_NAME;
                
        outputSamples(samples, 
                      pathway,
                      fg,
                      converter.getInstToVarMap(),
                      fileName);
    }
    
    private Set<Variable> getObservationVariables(FactorGraph fg,
                                                  DataType... dataTypes) {
        Set<Variable> variables = new HashSet<Variable>();
        for (Variable var : fg.getVariables()) {
            String name = var.getName();
            for (DataType dataType : dataTypes) {
                if (name.endsWith(dataType.toString())) {
                    variables.add(var);
                    break;
                }
            }
        }
        return variables;
    }
    
    private Set<Variable> getGeneVariables(FactorGraph fg,
                                           Set<Variable> observationVars) {
        Set<String> geneNames = new HashSet<String>();
        for (Variable obsVar : observationVars) {
            String name = obsVar.getName();
            int index = name.indexOf("_");
            geneNames.add(name.substring(0, index));
        }
        Set<Variable> rtn = new HashSet<Variable>();
        for (Variable var : fg.getVariables()) {
            String name = var.getName();
            int index = name.indexOf("_");
            if (index <= 0)
                continue;
            if (geneNames.contains(name.substring(0, index)))
                rtn.add(var);
        }
        return rtn;
    }
    
    private void outputSamples(List<Observation<Integer>> samples, 
                               GKInstance pathway, 
                               FactorGraph fg,
                               Map<GKInstance, Variable> instToVar,
                               String fileName) throws Exception, IOException {
        List<Variable> pathwayVars = new ReactomePathwayFGRunner().getPathwayVariables(pathway, 
                                                                                       instToVar,
                                                                                       fg);
        // Data variables
        Set<Variable> observationVars = getObservationVariables(fg,
                                                                DataType.CNV,
                                                                DataType.mRNA_EXP);
        Set<Variable> geneVars = getGeneVariables(fg, observationVars);
        Set<Variable> neededVars = new HashSet<Variable>(pathwayVars);
        neededVars.addAll(observationVars);
        neededVars.addAll(geneVars);
        
        // Do sorting
        List<Variable> varList = new ArrayList<Variable>(neededVars);
        Collections.sort(varList, new Comparator<Variable>() {
            public int compare(Variable var1, Variable var2) {
                return var1.getName().compareTo(var2.getName());
            }
        });
        
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample\tLogLikelihood");
        for (Variable var : varList)
            builder.append("\t").append(var.getName());
        fu.printLine(builder.toString());
        for (Observation sample : samples) {
            builder.setLength(0);
            builder.append(sample.getName());
            Map<Variable, Integer> varToState = sample.getVariableToAssignment();
            double loglikelihood = fg.getLogLikelihood(varToState);
            builder.append("\t").append(loglikelihood);
            for (Variable var : varList) {
                Integer state = varToState.get(var);
                builder.append("\t").append(state);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    private GKInstance getPathway(Long dbId) throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "gk_current_ver50", 
                                            "root", 
                                            "macmysql01");
//        XMLFileAdaptor dba = new XMLFileAdaptor("results/paradigm/ReactomePathways/SimpleYAP1Pathway.rtpj");
        GKInstance pathway = dba.fetchInstance(dbId);
        return pathway;
    }
    
    /**
     * Call this method to generate a list of Observations (aka samples) for the passed
     * FactorGraph object.
     * @return
     */
    public List<Observation<Integer>> generateSamples(int totalSamples) {
        if (factorGraph == null)
            throw new IllegalStateException("FactorGraph has not been specified.");

        // totalSamples should be the same as maxIteration.
        setMaxIteration(totalSamples);
        
        // To be returned.
        List<Observation<Integer>> samples = new ArrayList<Observation<Integer>>();
        resetCache();
        Observation<Integer> sample = initializeAssignment();
        // This sample is not the same one passed to the burn method
        sample = burn(sample);
        sample(sample, samples);
        
        return samples;
    }
    
    /**
     * Initialize the first sample.
     */
    private Observation<Integer> initializeAssignment() {
        // The first assignment
        Observation<Integer> observation = new Observation<Integer>();
        Map<Variable, Integer> varToAssign = new HashMap<Variable, Integer>();
        observation.setVariableToAssignment(varToAssign);
        FactorValueAssignment helper = new FactorValueAssignment();
        // We want to have a random start factor
        //        List<Factor> factors = new ArrayList<Factor>(factorGraph.getFactors());
        //        factors = MathUtilities.permutate(factors, randomizer);
        //        for (Factor factor : factors) {
        for (Factor factor : factorGraph.getFactors()) {
            Factor slicedFactor = slice(factor, 
                                        varToAssign,
                                        helper);
            if (slicedFactor == null)
                continue; // Just in case nothing new in this factor.
            initializeAssignment(slicedFactor,
                                 varToAssign,
                                 helper);
        }
        System.out.println("Initialize state: " + varToAssign);
        return observation;
    }
    
    /**
     * Create a new Factor from the specified factor object by assigning
     * values from the passed varToAssign.
     * @param factor
     * @param varToAssign
     * @return
     */
    private Factor slice(Factor factor,
                         Map<Variable, Integer> varToAssign,
                         FactorValueAssignment helper) {
        if (varToAssign.size() == 0)
            return factor; // There is no need to create a new factor object.
        // Get variables that have not been covered
        List<Variable> variables = factor.getVariables();
        List<Variable> slicedVars = new ArrayList<Variable>(variables);
        slicedVars.removeAll(varToAssign.keySet());
        if (slicedVars.size() == 0)
            return factor; // If nothing is shared with known variables, just return the current one.
        // To be returned
        Factor slicedFactor = new Factor();
        slicedFactor.setVariables(slicedVars);
        // No need to read out the values from remaining
        List<Double> slicedValues = new ArrayList<Double>();
        helper.setFactor(slicedFactor);;
        List<Map<Variable, Integer>> assignments = helper.iterate();
        Map<Variable, Integer> assign = new HashMap<Variable, Integer>();
        for (int i = 0; i < assignments.size(); i++) {
            Map<Variable, Integer> slicedAssign = assignments.get(i);
            assign.putAll(slicedAssign);
            // Get the rest
            for (Variable var : variables) {
                Integer value = varToAssign.get(var);
                if (value != null)
                    assign.put(var, value);
            }
            double factorValue = factor.getValue(assign);
            slicedValues.add(factorValue);
        }
        slicedFactor.setValues(slicedValues);
        return slicedFactor;
    }
    
    private void initializeAssignment(Factor factor,
                                      Map<Variable, Integer> varToAssignment,
                                      FactorValueAssignment helper) {
        helper.setFactor(factor);
        List<Map<Variable, Integer>> list = helper.iterate();
        // Want to get a random value
        int index = 0;
        while (true) {
//            index = randomizer.nextInt(0, list.size() - 1); // Note: the end points is included.
            double value = helper.getFactorValue(index);
            if (value > MINIMUM_PROB) {
                Map<Variable, Integer> assign = list.get(index);
                varToAssignment.putAll(assign);
                break;
            }
            index ++;
        }
    }
    
    private Set<String> getOutputEntityNames() {
        Set<String> entityNames = new HashSet<String>();
        // The following entity names are for YAP1 and WWTR1-stimulated gene expression
//        entityNames.add("CTGF [extracellular region]");
//        entityNames.add("NPPA(1-153) [extracellular region]");
        // The following entity names are for Transcriptional regulation of pluripotent stem cells
        String[] names = new String[] {
                "NANOG [nucleoplasm]",
                "POU5F1 [nucleoplasm]",
                "SOX2 [nucleoplasm]"
        };
        for (String name : names)
            entityNames.add(name);
        return entityNames;
    }
    
    @Test
    public void testGetTypeToSamples() throws IOException {
        Map<String, Set<String>> typeToSamples = getTypeToSamples(false);
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            System.out.println(type + ": " + samples.size());
        }
    }
    
    public Map<String, Set<String>> getTypeToSamples(boolean randomization) throws IOException {
        Map<String, Set<String>> typeToSamples = new HashMap<String, Set<String>>();
        FileUtility fu = new FileUtility();
        fu.setInput(SAMPLE_FILE_NAME);
        String line = fu.readLine();
        String[] names = line.split("\t");
        Set<String> values = new HashSet<String>();
        Set<String> entityNames = getOutputEntityNames();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            values.clear();
            for (int i = 2; i < tokens.length; i++) {
                if (entityNames.contains(names[i])) {
                    values.add(tokens[i]);
                }
            }
            if (values.size() > 1)
                continue;
            String value = values.iterator().next();
            if (value.equals("0"))
                InteractionUtilities.addElementToSet(typeToSamples, "0", tokens[0]);
            else if (value.equals("2")) // We don't need type == 1
                InteractionUtilities.addElementToSet(typeToSamples, "2", tokens[0]);
        }
        fu.close();
        if (!randomization)
            return typeToSamples;
        // Here is for randomization
        Map<String, Set<String>> rtn = new HashMap<String, Set<String>>();
        Set<String> allSamples = new HashSet<String>();
        for (Set<String> samples : typeToSamples.values())
            allSamples.addAll(samples);
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            if (allSamples.size() <= samples.size()) {
                rtn.put(type, allSamples);
            }
            else {
                Set<String> randomSamples = MathUtilities.randomSampling(allSamples, samples.size());
                rtn.put(type, randomSamples);
                allSamples.removeAll(randomSamples);
            }
        }
        return rtn;
    }
    
}
