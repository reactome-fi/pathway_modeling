/*
 * Created on May 28, 2014
 *
 */
package org.reactome.fi.pgm;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.junit.Test;
import org.reactome.factorgraph.ContinuousVariable;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.factorgraph.InferenceCannotConvergeException;
import org.reactome.factorgraph.InferenceType;
import org.reactome.factorgraph.LoopyBeliefPropagation;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.MutationEmpiricalFactorHandler;
import org.reactome.factorgraph.common.ObservationFactorHandler;
import org.reactome.factorgraph.common.ObservationFileLoader;
import org.reactome.factorgraph.common.ObservationHelper;
import org.reactome.factorgraph.common.ObservationRandomizer;
import org.reactome.fi.pgm.FIPGMConstructor.PGMType;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.UniProtAnalyzer;

/**
 * The main class for running the FI PGM model.
 * @author gwu
 *
 */
public class FIPGMRunner {
    private static final Logger logger = Logger.getLogger(FIPGMRunner.class);
    
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Please provide a choice (learn|inference) for the program.");
            System.exit(1);
        }
        FIPGMRunner runner = new FIPGMRunner();
        if (args[0].equals("learn"))
            runner.learnParameters();
        else
            runner.runInference();
    }
    
    @Test
    public void checkFIScores() throws IOException, MathException {
        String fileName = "results/FI_PGM/icgc_pancancer/PGM_FI_Inference_Results_072115.txt";
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double score = new Double(tokens[1]);
            if (score > 19.0d)
                geneToScore.put(tokens[0], new Double(tokens[1]));
        }
        fu.close();
        Map<String, Integer> geneToLength = new UniProtAnalyzer().loadGeneToProteinLength();
        Set<String> sharedGenes = InteractionUtilities.getShared(geneToScore.keySet(), geneToLength.keySet());
        List<Double> scores = new ArrayList<Double>();
        List<Double> lengths = new ArrayList<Double>();
        for (String gene : sharedGenes) {
            scores.add(geneToScore.get(gene));
            lengths.add(new Double(geneToLength.get(gene)));
//            System.out.println(gene + "\t" + geneToScore.get(gene) + "\t" + geneToLength.get(gene));
        }
        PearsonsCorrelation corr = MathUtilities.constructPearsonCorrelation(scores, lengths);
        System.out.println("Correlations: " + corr.getCorrelationMatrix().getEntry(0, 1));
        System.out.println("p-value: " + corr.getCorrelationPValues().getEntry(0, 1));
        
    }
    
    /**
     * Default constructor.
     */
    public FIPGMRunner() {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    @Test
    public void testEmpiricalDistribution() {
        // Sample from a normal distribution
        NormalDistribution normal = new NormalDistribution();
        double[] sample = new double[1000];
        for (int i = 0; i < sample.length; i++)
            sample[i] = normal.sample();
//        for (double value : sample)
//            System.out.println(value);
        // Check with EmpiricalDistribution
        EmpiricalDistribution empirical = new EmpiricalDistribution();
        empirical.load(sample);
        System.out.println("Mean: " + empirical.getNumericalMean());
        System.out.println("Variance: " + empirical.getNumericalVariance());
        // Check cumulative distribution
        System.out.println("\nValue\tNormalCP\tEmpiricalCP");
        for (double value : sample) {
            System.out.println(value + "\t" + 
                               normal.cumulativeProbability(value) + "\t" +
                               empirical.cumulativeProbability(value));
        }
    }
    
    @Test
    public void testFactorGraphLoad() throws Exception {
        // This is for test
        String fileName = "tmp/PGM_FI_072214.xml";
        FactorGraph fg = new FactorGraph();
        fg.importFG(new FileInputStream(fileName));
        System.out.println("Loaded factor graph: " + fg.getFactors().size() + " factors and " + fg.getVariables().size());
    }
    
    private Map<DataType, String> getEvidenceFiles() {
        return FIPGMConfiguration.getConfig().getTypeToEvidenceFile();
    }
    
    @Test
    public void learnParameters() throws Exception {
        
//        FIPGMConstructor constructor = new FIPGMConstructor();
//        constructor.setEvidenceFiles(getEvidenceFiles());
//        
//        FactorGraph fg = constructor.constructFactorGraph(PGMType.PairwiseMRF);
//        logger.info("Converted factor graph: " + fg.getFactors().size() + " factors and " + 
//                fg.getVariables().size() + " variables.");
//        Map<DataType, SharedEMFactors> typeToSharedFactors = constructor.getTypeToSharedFactors();
//        logger.info("Learning shared factors: " + typeToSharedFactors.size());
//        for (DataType type : typeToSharedFactors.keySet()) {
//            System.out.println(type + ": " + typeToSharedFactors.get(type));
//        }
//        // This is for test
////        String fileName = "tmp/PGM_FI_072214.xml";
////        fg.exportFG(new FileOutputStream(fileName));
////        if (true)
////            return;
//        
//        ExpectationMaximization em = new ExpectationMaximization();
//        LoopyBeliefPropagation lbp = (LoopyBeliefPropagation) em.getInferencer();
//        lbp.setDebug(true);
//        em.setDebug(true);
//        lbp.setTolerance(1.0e-5);
//        lbp.setMaxIteration(75);
//        lbp.setUseLogSpace(false);
//        em.setTolerance(1.0e-5);
//        em.setMaxIteration(25);
//        List<Observation<Float>> observations = constructor.getObservationLoader().getObservations();
//        em.setEvidences(Observation.convertListFromFloatToInteger(observations));
//        logger.info("Starting learning...");
//        long time1 = System.currentTimeMillis();
//        em.learn(fg, new ArrayList<EMFactor>(typeToSharedFactors.values()));
//        long time2 = System.currentTimeMillis();
//        logger.info("Learning is done: " + (time2 - time1) / 1000.0d + " seconds.");
//        StringBuilder builder = new StringBuilder();
//        // Save the learned parameters into a file
//        String fileName = FIPGMConfiguration.getConfig().getProperties().get("learnedResult");
//        if (fileName == null)
//            fileName = FIPGMConfiguration.RESULT_DIR + "LearnedParameters.txt";
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        for (DataType type : typeToSharedFactors.keySet()) {
//            SharedEMFactors factor = typeToSharedFactors.get(type);
//            builder.append(type);
//            for (double value : factor.getValues())
//                builder.append("\t").append(value);
//            logger.info("Learned parameters: " + builder.toString());
//            System.out.println(builder.toString());
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
    }
    
//    private Map<Long, String> getVarLabelToName(FIPGMConstructor constructor) {
//        Map<Long, String> varLabelToName = new HashMap<Long, String>();
//        Map<String, Var> nameToVar = constructor.getNameToVar();
//        for (String name : nameToVar.keySet()) {
//            Var var = nameToVar.get(name);
//            varLabelToName.put(var.label(), name);
//        }
//        return varLabelToName;
//    }
//    
    
    @Test
    public void testRandomGeneration() throws Exception {
        FIPGMConstructor constructor = getPGMConstructor();
        
        PGMType type = PGMType.PairwiseMRF;
        //        type = PGMType.NearestNeighborGibbs;
        FactorGraph fg = constructor.constructFactorGraph(type);
        logger.info("Converted factor graph: " + fg.getFactors().size() + " factors and " + 
                fg.getVariables().size() + " variables.");
        // Want to use an Observation for testing purposes
        List<Observation<Number>> observations = constructor.getObservationLoader().getObservations();
        logger.info("Total observations: " + observations.size());
        countAssignments(observations);
        ObservationHelper helper = new ObservationHelper();
        List<Observation<Number>> copy = new ArrayList<Observation<Number>>(observations);
        helper.filterObservationsToHaveSharedDataTypes(copy);
        logger.info("Observations having shared data types: " + copy.size());
        
        // Just to test randomization
        ObservationRandomizer randomizer = new ObservationRandomizer();
        randomizer.setNumberOfPermutation(100);
        long time11 = System.currentTimeMillis();
        
        List<Observation<Number>> randomObservations = randomizer.createRandomObservations(observations);
        logger.info("Total randomized observations: " + randomObservations.size());
        long time12 = System.currentTimeMillis();
        logger.info("Time for randomization: " + (time12 - time11) + " ms");
        countAssignments(randomObservations);
        //        // Check the first randomized sample
        //        Observation<Number> randomObservation = randomObservations.get(0);
        //        for (VariableAssignment<Number> varAssgn : randomObservation.getVariableAssignments()) {
        //            System.out.println(varAssgn.getVariable().getName() + "\t" + varAssgn.getAssignment());
        //        }
    
    }

    private void countAssignments(List<Observation<Number>> observations) {
        // Check the total assignments
        int count0 = 0;
        int count1 = 0;
        for (Observation<Number> obs : observations) {
            for (VariableAssignment<Number> varAssgn : obs.getVariableAssignments()) {
                if (varAssgn.getAssignment().intValue() == 0)
                    count0 ++;
                else
                    count1 ++;
            }
        }
        System.out.println("count 0: " + count0);
        System.out.println("count 1: " + count1);
        System.out.println("ratio: " + (double)count1 / (count0 + count1));
    }
   
    @Test
    public void runToyInference() throws Exception {
        FIPGMConstructor constructor = getPGMConstructor();
//        ObservationFileLoader fileLoader = constructor.getObservationLoader();
//        ObservationFactorHandler factorHandler = new DiscreteObservationFactorhandler();
//        fileLoader.setObservationFactorHandler(DataType.mRNA_EXP, factorHandler);
//        fileLoader.setObservationFactorHandler(DataType.Mutation, factorHandler);
        
        Set<String> fis = new HashSet<String>();
        fis.add("A\tB");
        fis.add("B\tC");
        fis.add("B\tF");
        fis.add("C\tD");
        fis.add("F\tC");
        for (int i = 0; i < 100; i++) {
            fis.add("B\tB" + i);
        }
        for (int i = 0; i < 100; i++) {
            fis.add("C\tC" + i);
        }
        for (int i = 0; i < 100; i++)
            fis.add("F\tF" + i);
        PGMType type = PGMType.NearestNeighborGibbs;
        type = PGMType.PairwiseMRF;
        FactorGraph fg = constructor.constructFactorGraph(PGMType.getValuesAssigner(type),
                                                          fis);
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setFactorGraph(fg);
        lbp.runInference();
        List<Variable> sortedVariables = getSortedVariables(fg);
        Map<Variable, double[]> varToPrior = new HashMap<Variable, double[]>();
        System.out.println("Prior:");
        for (Variable var : sortedVariables) {
            if (var instanceof ContinuousVariable)
                continue;
            String name = var.getName();
            if (name.length() > 1 && !name.endsWith(DataType.Mutation.toString()) && !name.endsWith(DataType.mRNA_EXP.toString()))
                continue;
            double[] belief = var.getBelief();
            System.out.println(var.getName() + "\t" + belief[0] + "\t" + belief[1]);
            double[] prior = new double[belief.length];
            System.arraycopy(belief, 0, prior, 0, belief.length);
            varToPrior.put(var, prior);
        }
//        if (true)
//            return;
        List<Observation<Number>> observations = constructor.getObservationLoader().getObservations();
        Collections.sort(observations, new Comparator<Observation>() {
            public int compare(Observation obs1, Observation obs2) {
                int index = obs1.getName().indexOf("_");
                int name1 = new Integer(obs1.getName().substring(index + 1));
                index = obs2.getName().indexOf("_");
                int name2 = new Integer(obs2.getName().substring(index + 1));
                return name1 - name2;
            }
        });
        // First observation
        Observation obs0 = new ObservationHelper().createBaseObservation(observations, 
                                                                         null,
                                                                         0);
        lbp.setObservation(obs0);
        lbp.runInference();
        for (Variable var : sortedVariables) {
            if (var instanceof ContinuousVariable)
                continue;
            String name = var.getName();
            if (name.length() > 1 && !name.endsWith(DataType.Mutation.toString()) && !name.endsWith(DataType.mRNA_EXP.toString()))
                continue;
            double[] belief = var.getBelief();
            double[] prior = new double[belief.length];
            System.arraycopy(belief, 0, prior, 0, belief.length);
            varToPrior.put(var, prior);
        }
        
        for (Observation observation : observations) {
            System.out.println("\nSample: " + observation.getName());
            lbp.setObservation(observation);
            lbp.runInference();
            Map<Variable, Number> varToAssgn = observation.getVariableToAssignment();
            for (Variable var : sortedVariables) {
                Number value = varToAssgn.get(var);
                if (value == null)
                    continue;
                System.out.println(var.getName() + "\t" + varToAssgn.get(var));
            }
            
            for (Variable var : sortedVariables) {
                String name = var.getName();
                if (name.length() > 1)
                    continue;
                double[] belief = var.getBelief();
                double[] prior = varToPrior.get(var);
                double logRatio = Math.log10(belief[1] * prior[0] / (belief[0] * prior[1]));
                System.out.println(var.getName() + "\t" + belief[0] + "\t" + belief[1] + "\t" + logRatio);
            }
        }
    }

    private FIPGMConstructor getPGMConstructor() {
        FIPGMConstructor constructor = new FIPGMConstructor();
        constructor.setEvidenceFiles(getEvidenceFiles());
        ObservationFileLoader fileLoader = constructor.getObservationLoader();
        Map<DataType, ObservationFactorHandler> typeToHandler = FIPGMConfiguration.getConfig().getTypeToHandler();
        if (typeToHandler != null && typeToHandler.containsKey(DataType.Mutation))
            fileLoader.setObservationFactorHandler(DataType.Mutation,
                                                   typeToHandler.get(DataType.Mutation));
        else {
            //        fileLoader.setObservationFactorHandler(DataType.Mutation,
            //                                               new DiscreteObservationFactorhandler());
            MutationEmpiricalFactorHandler mutationHandler = new MutationEmpiricalFactorHandler();
            mutationHandler.setConfiguration(FIPGMConfiguration.getConfig());
            fileLoader.setObservationFactorHandler(DataType.Mutation,
                                                   mutationHandler);
        }
//        fileLoader.setObservationFactorHandler(DataType.Mutation,
//                                               mutationHandler);
        
//        fileLoader.setObservationFactorHandler(DataType.mRNA_EXP,
//                                               new DiscreteObservationFactorhandler());
        // Try CLG for mRNA Expression
//        fileLoader.setObservationFactorHandler(DataType.mRNA_EXP,
//                                               new CLGObservationFactorHandler());
        // Try empirical distribution for mRNA expression
//        fileLoader.setObservationFactorHandler(DataType.CNV,
//                                               new DiscreteObservationFactorhandler());
        return constructor;
    }
    
    private List<Variable> getSortedVariables(FactorGraph fg) {
        List<Variable> sortedVariables = new ArrayList<Variable>((fg.getVariables()));
        Collections.sort(sortedVariables, new Comparator<Variable>() {
           public int compare(Variable var1, Variable var2) { 
               return var1.getName().compareTo(var2.getName());
           }
        });
        return sortedVariables;
    }
    
    /**
     * Filter observations to have all types information.
     * @param observations
     */
    private void filterObservations(List<Observation<Number>> observations) {
        new ObservationHelper().filterObservationsToHaveSharedDataTypes(observations);
    }

    /**
     * To avoid interference among different data types (e.g. a large mRNA expression may be cancelled out 
     * by a lacking of mutation in the central dogma model), we run inference based on individual data types,
     * and then add scores together.
     * @throws Exception
     */
    @Test
    public void runInferenceBasedOnDataTypes() throws Exception {
//        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
//        
//        Set<String> fis = config.getFIs();
//        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
//        Map<String, Set<String>> geneToParnters = InteractionUtilities.generateProteinToPartners(fis);
//        String target = "TP53";
//        Set<String> gnbPartners = geneToParnters.get(target);
////        gnbPartners = MathUtilities.randomSampling(genes, gnbPartners.size());
//        // Randomized target
////        target = MathUtilities.randomSampling(genes, 1).iterator().next();
//        
//        Set<String> gnbFIs = new HashSet<String>();
//        for (String partner : gnbPartners)
//            gnbFIs.add(target + "\t" + partner);
//        config.setFIs(gnbFIs);
        
        //        Set<String> randomFIs = InteractionUtilities.generateRandomFIs(config.getFIs());
        //        config.setFIs(randomFIs);
        
        FIPGMConstructor constructor = getPGMConstructor();
        // Initialize an inference algorithm
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
//        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.setInferenceType(InferenceType.SUM_PRODUCT);
        lbp.setUseLogSpace(true);
        lbp.setTolerance(1.0e-5);
//        lbp.setMaxIteration(100);
        lbp.setDebug(true);
        
        PGMType type = PGMType.PairwiseMRF;
        //        type = PGMType.NearestNeighborGibbs;
        FactorGraph fg = constructor.constructFactorGraph(type);
        logger.info("Converted factor graph: " + fg.getFactors().size() + " factors and " + 
                fg.getVariables().size() + " variables.");
        lbp.setFactorGraph(fg);
        // Get a list of variables for output
        List<Variable> sortedVariables = getSortedVariables(fg);
        List<Observation<Number>> observations = constructor.getObservationLoader().getObservations();
        logger.info("Total observations: " + observations.size());
        filterObservations(observations);
        logger.info("After filtering: " + observations.size());
        
        // Save the output into a file
        String fileName = FIPGMConfiguration.getConfig().getProperties().get("inferenceResult");
        if (fileName == null)
            throw new IllegalStateException("inferenceResult file has not been specified!");
        logger.info("inferenceResult: " + fileName);
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        int count = 0;
        
        // We want to perform impact inference based on data types one by one
        // and them sum them together.
        Map<DataType, Map<Variable, double[]>> dataTypeToVarToPrior = new HashMap<DataType, Map<Variable,double[]>>();
        // Total data types in the observations
        ObservationHelper observationHelper = new ObservationHelper();
        Set<DataType> dataTypes = observationHelper.getDataTypesFromObservations(observations);
        logger.info("Total data types: " + dataTypes.size());
        for (DataType dataType : dataTypes) {
            Map<Variable, double[]> varToPrior = new HashMap<Variable, double[]>();
            dataTypeToVarToPrior.put(dataType, varToPrior);
            Observation<Number> baseObs = observationHelper.createBaseObservation(observations,
                                                                                  dataType,
                                                                                  0);
            lbp.setObservation(baseObs);
            lbp.runInference();
            for (Variable var : fg.getVariables()) {
                if (var instanceof ContinuousVariable)
                    continue; // Don't support this query.
                double[] belief = var.getBelief();
                double[] prior = new double[belief.length];
                System.arraycopy(belief, 0, prior, 0, belief.length);
                varToPrior.put(var, prior);
            }
        }
        
        //        // This is for testing
        //        Observation<Number> baseObs = createBaseObservation(observations, DataType.mRNA_EXP);
        //        List<VariableAssignment<Number>> varAssgns = baseObs.getVariableAssignments();
        //        for (VariableAssignment<Number> varAssgn : varAssgns) {
        //            Variable var = varAssgn.getVariable();
        //            if (var.getName().equals("GNB1_Mutation"))
        //                varAssgn.setAssignment(-3.505000114440918);
        //            else if(var.getName().equals("GNB1_mRNA_EXP"))
        //                varAssgn.setAssignment(2.3319245968204547);
        //        }
        //        observations.clear();
        //        baseObs.setName("TCGA-20-0991");
        //        observations.add(baseObs);
        
        List<DataType> sortedDataTypes = new ArrayList<DataType>(dataTypes);
        Collections.sort(sortedDataTypes);
        for (Observation<Number> observation : observations) {
//                        if (!observation.getName().equals("TCGA-20-0991"))
//                            continue;
            //            if (!observation.getName().equals("TCGA-24-0980"))
            //                continue; // TTN has high value
            //                        if (!observation.getName().equals("TCGA-10-0927"))
            //                            continue; // Running extremely long time!!!
//            if (!observation.getName().equals("d3d65db3-36f9-41c7-8e5e-1683ce94dfcb"))
//                continue; // TTN has four mutations here. We want to see the effects of TTN with meaning FI scores
            long time1 = System.currentTimeMillis();
            fu.printLine("Sample:" + observation.getName());
            logger.info("Sample: " + observation.getName());
            
//            for (VariableAssignment<Number> varAssgn : observation.getVariableAssignments()) {
////                if (!varAssgn.getVariable().getName().startsWith("GNB1") &&
////                        !varAssgn.getVariable().getName().startsWith("MC5R"))
////                    continue;
//                if (!varAssgn.getVariable().getName().endsWith("_" + DataType.Mutation))
//                    continue;
//                EmpiricalDistribution dist = varAssgn.getDistribution();
//                if (dist == null)
//                    System.out.println(varAssgn.getVariable().getName() + "\t" + varAssgn.getAssignment());
//                else {
//                    double pvalue = dist.cumulativeProbability(varAssgn.getAssignment().doubleValue());
//                    System.out.println(varAssgn.getVariable().getName() + "\t" + varAssgn.getAssignment() + "\t" + pvalue);
//                }
//            }
            // Split the observation into multiple ones according to data types
            Map<DataType, Observation<Number>> dataTypeToObs = observationHelper.splitObservationBasedOnDataTypes(observation);
            logger.info("Split samples: " + dataTypeToObs.size());
            String sampleLine = "";
            String header = "";
            Map<Variable, String> varToLine = new HashMap<Variable, String>();
            boolean isFirst = true;
            Map<Variable, Double> varToScore = new HashMap<Variable, Double>();
            for (DataType dataType : sortedDataTypes) {
                Observation<Number> dataObs = dataTypeToObs.get(dataType);
                lbp.setObservation(dataObs);
                //        lbp.setUpdateViaFactors(true);
                lbp.runInference();
                Map<String, Number> geneVarToValue = observationHelper.getGeneToAssignment(dataObs, dataType);
                long time2 = System.currentTimeMillis();
                logger.info("Sample: " + dataObs.getName());
                //            logger.info("LogZ: " + lbp.calculateLogZ());
                logger.info("Total time: " + (time2 - time1) / 1000.0d + " seconds.");
                //        Map<Variable, Integer> maxima = lbp.findMaximum();
                sampleLine += dataObs.getName() + "\t";
                if (isFirst) 
                    header +=("Genes\tBelief(0)\tBelief(1)\tPrior(0)\tPrior(1)\tLogRatio\tLinkedVar\tObservation");
                else
                    header +=("\tBelief(0)\tBelief(1)\tLogRatio\tObservation");
                Map<Variable, double[]> varToPrior = dataTypeToVarToPrior.get(dataType);
                for (Variable var : sortedVariables) {
                    String name = var.getName();
                    if (name.contains("_"))
                        continue;
                    Set<Variable> linkedVars = getLinkedVariables(var);
                    double[] belief = var.getBelief();
                    double[] prior = varToPrior.get(var);
                    double logRatio = Math.log10(belief[1] * prior[0] / (belief[0] * prior[1]));
                    Double oldValue = varToScore.get(var);
                    if (oldValue == null)
                        varToScore.put(var, logRatio);
                    else
                        varToScore.put(var, oldValue + logRatio);
                    String line = null;
                    Number value = geneVarToValue.get(name + "_" + dataType);
                    if (isFirst) {
                        line = (name + "\t" + belief[0] + "\t" + belief[1] + "\t" +  + prior[0] +
                                "\t" + prior[1] + "\t" + logRatio + "\t" + linkedVars.size() + "\t" + value);
                    }
                    else
                        line = (belief[0] + "\t" + belief[1] + "\t" + logRatio + "\t" + value);
                    String oldLine = varToLine.get(var);
                    if (oldLine == null)
                        varToLine.put(var, line);
                    else
                        varToLine.put(var, oldLine + "\t" + line);
                }
                if (isFirst)
                    isFirst = false;
            }
            // Output
            fu.printLine(sampleLine);
            fu.printLine(header + "\tSum");
            for (Variable var : sortedVariables) {
                String line = varToLine.get(var);
                if (line == null)
                    continue;
                Double score = varToScore.get(var);
                fu.printLine(line + "\t" + score);
            }
            count ++;
            logger.info("Counts: " + count);
//            if (count == 1)
//                break;
        }
        fu.close();
//        convertToMatrix(fileName);
    }
    
    @Test
    public void runInference() throws IOException, InferenceCannotConvergeException {
        
        FIPGMConstructor constructor = getPGMConstructor();
        
        PGMType type = PGMType.PairwiseMRF;
//        type = PGMType.NearestNeighborGibbs;
        FactorGraph fg = constructor.constructFactorGraph(type);
        logger.info("Converted factor graph: " + fg.getFactors().size() + " factors and " + 
                            fg.getVariables().size() + " variables.");
        // Get a list of variables for output
        List<Variable> sortedVariables = getSortedVariables(fg);
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (Variable var : fg.getVariables())
            nameToVar.put(var.getName(), var);
        
//        // This is just a test
//        Var tp53 = nameToVar.get("TP53");
//        infAlg.clamp(tp53.label(), 1L);
//        Var egfr = nameToVar.get("EGFR");
//        infAlg.clamp(egfr.label(), 1L);
        // Want to use an Observation for testing purposes
        List<Observation<Number>> observations = constructor.getObservationLoader().getObservations();
        logger.info("Total observations: " + observations.size());
        
        // Just to test randomization
        ObservationRandomizer randomizer = new ObservationRandomizer();
        randomizer.setNumberOfPermutation(100);
        long time11 = System.currentTimeMillis();
        List<Observation<Number>> randomObservations = randomizer.createRandomObservations(observations);
        logger.info("Total randomized observations: " + randomObservations.size());
        long time12 = System.currentTimeMillis();
        logger.info("Time for randomization: " + (time12 - time11) + " ms");
        // Check the first randomized sample
        Observation<Number> randomObservation = randomObservations.get(0);
        for (VariableAssignment<Number> varAssgn : randomObservation.getVariableAssignments()) {
            System.out.println(varAssgn.getVariable().getName() + "\t" + varAssgn.getAssignment());
        }
        
//        if (true)
//            return;
        
        filterObservations(observations);
        logger.info("After filtering: " + observations.size());
        
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.setDebug(true);
        lbp.setFactorGraph(fg);
        // Save the output into a file
        String fileName = FIPGMConfiguration.getConfig().getProperties().get("inferenceResult");
        if (fileName == null)
            throw new IllegalStateException("inferenceResult file has not been specified!");
        logger.info("inferenceResult: " + fileName);
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        int count = 0;
        // Calculate prior
        lbp.setUseLogSpace(true);
        lbp.setTolerance(1.0e-5);
        Map<Variable, double[]> varToPrior = new HashMap<Variable, double[]>();
      
//      lbp.setInferenceType(InferenceType.MAX_PRODUCT);
//        lbp.setInferenceType(InferenceType.SUM_PRODUCT);
//        lbp.runInference();
//        for (Variable var : fg.getVariables()) {
//            if (var instanceof ContinuousVariable)
//                continue; // Don't support this query.
//            double[] belief = var.getBelief();
//            double[] prior = new double[belief.length];
//            System.arraycopy(belief, 0, prior, 0, belief.length);
//            varToPrior.put(var, prior);
////            Set<Variable> variables = getLinkedVariables(var);
////            fu.printLine(var.getName() + "\t" + belief[0] + "\t" + belief[1] + "\t" + variables.size());
//        }
      
        Observation<Number> baseObs = new ObservationHelper().createBaseObservation(observations,
                                                                                    null,
                                                                                    0);
        lbp.setObservation(baseObs);
        lbp.runInference();
        for (Variable var : fg.getVariables()) {
            if (var instanceof ContinuousVariable)
                continue; // Don't support this query.
            double[] belief = var.getBelief();
            double[] prior = new double[belief.length];
            System.arraycopy(belief, 0, prior, 0, belief.length);
            varToPrior.put(var, prior);
        }
//        if (true)
//            return;
        
        // This is for testing
//        List<VariableAssignment<Number>> varAssgns = baseObs.getVariableAssignments();
//        for (VariableAssignment<Number> varAssgn : varAssgns) {
//            Variable var = varAssgn.getVariable();
//            if (var.getName().equals("GNB1_Mutation"))
//                varAssgn.setAssignment(-3.505000114440918);
//            else if(var.getName().equals("GNB1_mRNA_EXP"))
//                 varAssgn.setAssignment(2.3319245968204547);
//        }
//        observations.clear();
//        observations.add(baseObs);
        Map<Variable, List<Double>> varToScores = new HashMap<Variable, List<Double>>();
        for (Observation<Number> observation : randomObservations) {
//            if (!observation.getName().equals("TCGA-20-0991"))
//                continue;
//            if (!observation.getName().equals("TCGA-24-0980"))
//                continue; // TTN has high value
//            if (!observation.getName().equals("TCGA-10-0927"))
//                continue; // Running extremely long time!!!
            long time1 = System.currentTimeMillis();
            fu.printLine("Sample:" + observation.getName());
            
            lbp.setObservation(observation);
            //        lbp.setUpdateViaFactors(true);
            lbp.runInference();
            
            long time2 = System.currentTimeMillis();
            logger.info("Sample: " + observation.getName());
//            logger.info("LogZ: " + lbp.calculateLogZ());
            logger.info("Total time: " + (time2 - time1) / 1000.0d + " seconds.");
            //        Map<Variable, Integer> maxima = lbp.findMaximum();
            fu.printLine("Genes\tBelief(0)\tBelief(1)\tPrior(0)\tPrior(1)\tLogRatio\tLinkedVar\tMutation");
            for (Variable var : sortedVariables) {
                String name = var.getName();
                if (name.contains("_"))
                    continue;
                Set<Variable> linkedVars = getLinkedVariables(var);
                double[] belief = var.getBelief();
                double[] prior = varToPrior.get(var);
                double logRatio = Math.log10(belief[1] * prior[0] / (belief[0] * prior[1]));
                //Variable mutVar = nameToVar.get(name + "_" + DataType.Mutation);
                Variable mutVar = nameToVar.get(name + "_" + DataType.CNV);
                VariableAssignment<Number> mutAssgn = observation.getVariableAssignment(mutVar);
                fu.printLine(name + "\t" + belief[0] + "\t" + belief[1] + "\t" +  + prior[0] +
                             "\t" + prior[1] + "\t" + logRatio + "\t" + linkedVars.size() + "\t" + 
                             (mutAssgn == null ? "null" : mutAssgn.getAssignment()));
                List<Double> scores = varToScores.get(var);
                if (scores == null) {
                    scores = new ArrayList<Double>();
                    varToScores.put(var, scores);
                }
                scores.add(logRatio);
            }
            count ++;
//            break;
            if (count == 5)
                break;
        }
        fu.close();
        fu.setOutput("tmp/output.txt");
        fu.printLine("Gene\tScores\tMean");
        for (Variable var : varToScores.keySet()) {
            List<Double> scores = varToScores.get(var);
            Double mean = MathUtilities.calculateMean(scores);
            fu.printLine(var.getName()  + "\t" + scores + "\t" + mean);
        }
        fu.close();
//        convertToMatrix(fileName);
    }
    
    private Set<Variable> getLinkedVariables(Variable var) {
        Set<Variable> variables = new HashSet<Variable>();
        for (Factor factor : var.getFactors()) {
            variables.addAll(factor.getVariables());
        }
        variables.remove(var);
        return variables;
    }
    
    public void convertToMatrix(String srcFileName,
                                String targetFileName,
                                int scoreCol) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        // Get all samples
        String line = null;
        List<String> samples = new ArrayList<String>();
        List<String> genes = new ArrayList<String>();
        List<List<Double>> values = new ArrayList<List<Double>>();
        // Used to hold values
        List<Double> valueList = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Sample:")) {
                int index = line.indexOf(":");
                String sample = line.substring(index + 1);
                samples.add(sample);
                if (valueList != null)
                    values.add(valueList);
                valueList = new ArrayList<Double>();
                // Escape two lines
                for (int i = 0; i < 2; i++)
                    line = fu.readLine();
            }
            else { // value line
                String[] tokens = line.split("\t");
                if (samples.size() == 1) {
                    // First sample: get the genes
                    genes.add(tokens[0]);
                }
                valueList.add(new Double(tokens[scoreCol])); // Take the last value
            }
        }
        fu.close();
        values.add(valueList); // Add the last one
        // Keep the original file 
        exportToMatrix(targetFileName, samples, genes, values);
    }

    private void exportToMatrix(String fileName, 
                                List<String> samples,
                                List<String> genes,
                                List<List<Double>> values) throws IOException {
        // Save the original file for debugging
        File file = new File(fileName);
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Gene\tSUM");
        for (String sample : samples) {
            builder.append("\t").append(sample);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (int i = 0; i < genes.size(); i++) {
            String gene = genes.get(i);
            builder.append(gene);
            double sum = 0.0d;
            for (List<Double> valueList : values) {
                sum += valueList.get(i);
            }
            builder.append("\t").append(sum);
            for (List<Double> valueList : values) {
                builder.append("\t").append(valueList.get(i));
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private Map<String, Factor> getGeneToFactor(FactorGraph fg) {
        Map<String, Factor> geneToFactor = new HashMap<String, Factor>();
        for (Factor factor : fg.getFactors()) {
            if (factor.getVariables().size() == 1)
                geneToFactor.put(factor.getName(), factor);
        }
        return geneToFactor;
    }
    
    @Test
    public void runInferenceWithTemplate() throws IOException, InferenceCannotConvergeException {
        FIPGMConstructor constructor = new FIPGMConstructor();
        FactorGraph fg = constructor.constructFactorGraph(PGMType.PairwiseMRF);
        Map<String, Factor> geneToFactor = getGeneToFactor(fg);
        System.out.println("Converted factor graph: " + fg.getFactors().size() + " factors and " + 
                fg.getVariables().size() + " variables.");
        
        ObservationFileLoader fileLoader = new FIObservationFileLoader();
//        String cnvFileName = "tmp/ov.CNV.txt";
//        String mRNAExpFileName = "tmp/ov.mRNA.txt";
//        String cnvFileName = FIPGMConstants.DATA_DIR + "ov.CNV.100.txt";
//        String mRNAExpFileName = FIPGMConstants.DATA_DIR + "ov.mRNA.100.txt";
        String cnvFileName = FIPGMConfiguration.DATA_DIR + "ov.CNV.10.txt";
        String mRNAExpFileName = FIPGMConfiguration.DATA_DIR + "ov.mRNA.10.txt";
        
        Map<String, Map<String, Float>> cnvSampleToGeneToAssgn = fileLoader.loadObservationData(cnvFileName,
                                                                                                DataType.CNV);
        Map<String, Map<String, Float>> mRNAExpSampleToGeneToAssgn = fileLoader.loadObservationData(mRNAExpFileName,
                                                                                                    DataType.mRNA_EXP);
        // Get a union of samples
        Set<String> samples = new HashSet<String>(cnvSampleToGeneToAssgn.keySet());
        samples.addAll(mRNAExpSampleToGeneToAssgn.keySet());
        
        long time1 = System.currentTimeMillis();
        CentralDogmaSubGraphHandler templateHandler = new CentralDogmaSubGraphHandler();
        templateHandler.enableDataType(DataType.CNV);
        templateHandler.enableDataType(DataType.mRNA_EXP);
        templateHandler.setInferenceType(InferenceType.MAX_PRODUCT);
        Map<DataType, Integer> typeToAssgn = new HashMap<DataType, Integer>();
        for (String sample : samples) {
//            if (!sample.equals("TCGA-30-1880"))
//                continue;
            Map<String, Float> cnvGeneToAssgn = cnvSampleToGeneToAssgn.get(sample);
            Map<String, Float> mRNAExpGeneToAssgn = mRNAExpSampleToGeneToAssgn.get(sample);
            for (String gene : geneToFactor.keySet()) {
                typeToAssgn.clear();
                if (cnvGeneToAssgn != null)
                    typeToAssgn.put(DataType.CNV, cnvGeneToAssgn.get(gene).intValue());
                if (mRNAExpGeneToAssgn != null)
                    typeToAssgn.put(DataType.mRNA_EXP, mRNAExpGeneToAssgn.get(gene).intValue());
                double[] belief = templateHandler.calculateBelief(typeToAssgn);
                Factor factor = geneToFactor.get(gene);
                double[] values = factor.getValues();
                templateHandler.mergeBeliefToFactorValues(values, belief);
            }
        }
        checkFactorValues(geneToFactor);
        if (true)
            return;
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setDebug(true);
        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.setFactorGraph(fg);
        lbp.setUseLogSpace(false);
        lbp.runInference();

        long time2 = System.currentTimeMillis();
//        Use all data files, NaN was thrown when calculating logZ!
        System.out.println("LogZ: " + lbp.calculateLogZ());
        System.out.println("Total time: " + (time2 - time1) / 1000.0d + " seconds.");
        
        Map<Variable, Integer> maxAssign = lbp.findMaximum();
        int totalStates = 0;
        for (Variable var : fg.getVariables()) {
            // Remove some weird entries
            if (var.getName().contains("_"))
                continue;
            Integer assgn = maxAssign.get(var);
            double[] belief = var.getBelief();
            System.out.println(var.getName() + "\t" + assgn + "\t" + 
                               belief[0] + "\t" + belief[1]);
            totalStates += assgn;
        }
        System.out.println("Total states: " + totalStates);
    }
    
    private void checkFactorValues(Map<String, Factor> geneToFactor) {
        for (Factor factor : geneToFactor.values())
            System.out.println(factor);
    }
    
}
