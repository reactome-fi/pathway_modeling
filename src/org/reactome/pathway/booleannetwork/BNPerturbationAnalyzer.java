/*
 * Created on Jul 27, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.booleannetwork.BooleanNetwork;
import org.reactome.booleannetwork.BooleanNetworkUtilities;
import org.reactome.booleannetwork.BooleanVariable;
import org.reactome.booleannetwork.FuzzyLogicSimulator;
import org.reactome.booleannetwork.IdentityFunction;
import org.reactome.booleannetwork.SimulationComparator;
import org.reactome.booleannetwork.SimulationConfiguration;
import org.reactome.booleannetwork.SimulationResults;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * Use Boolean network based simulation to analyze perturbation caused by drugs or somatic muations.
 * @author wug
 *
 */
public class BNPerturbationAnalyzer {

	public BNPerturbationAnalyzer() {
	}
	
	public MySQLAdaptor getDBA() throws Exception {
		MySQLAdaptor dba = new MySQLAdaptor("localhost",
										   "reactome_63_plus_i", 
										   "root",
										   "macmysql01");
		return dba;
	}
	
	@Test
	public void checkInitials() throws Exception {
	    FileUtility.initializeLogging();
	    Long dbId = 5693567L;

	    BooleanNetwork network = getNetwork(dbId);
	    Map<BooleanVariable, Number> varToInit = createInitials(network, 1.0d);
	    Integer count = 0;
	    for (BooleanVariable var : varToInit.keySet()) {
	        Number init = varToInit.get(var);
	        if (init.intValue() == 1)
	            count ++;
	        System.out.println(var.getId() + "\t" + var.getName() + "\t" + init);
	    }
	    System.out.println("Total variabls having 1: " + count);
	}
	
	/**
	 * Check how stable the differences between two simulation results by expanding.
	 * @throws Exception
	 */
	public void checkPerturbationScores(SimulationResults perturbed,
	                                    SimulationResults base) throws Exception {
	    int initStep = Math.max(perturbed.getTimeSteps(), base.getTimeSteps());
	    int increase = 25;
	    int maxStep = 500;
	    List<Integer> steps = new ArrayList<>();
	    
	    // For later use
	    SimulationResults perturbedCopy = perturbed.copy();
	    SimulationResults baseCopy = base.copy();
	    
	    SimulationComparator comparator = new SimulationComparator();
	    Map<BooleanVariable, List<Double>> varToDiffs = new HashMap<>();
	    for (int step = initStep; step < maxStep; step += increase) {
	        steps.add(step);
	        Map<BooleanVariable, Double> varToDiff = comparator.calculateDiff(perturbed,
	                                                                          base,
	                                                                          step);
	        // Record these values
	        varToDiff.forEach((var, diff) -> {
	            varToDiffs.compute(var, (key, diffs) -> {
	                if (diffs == null)
	                    diffs = new ArrayList<>();
	                diffs.add(diff);
	                return diffs;
	            });
	        });
	    }
	    
	    // Check the implemented method in SimulationComparator
	    Map<BooleanVariable, Double> varToDiff = comparator.calculateDiff(perturbedCopy,
	                                                                      baseCopy,
	                                                                      increase,
	                                                                      0.001d,
	                                                                      maxStep);
	    
	    // Want to output all results
	    System.out.println("\nOutput from check perturbation:");
	    StringBuilder builder = new StringBuilder();
	    builder.append("id\tVariable");
	    steps.forEach(step -> builder.append("\t").append(step));
	    builder.append("\tConverged");
	    System.out.println(builder.toString());
	    builder.setLength(0);
	    
	    varToDiffs.forEach((var, diffs) -> {
	        builder.append(var.getId()).append("\t");
	        builder.append(var.getName());
	        diffs.forEach(diff -> builder.append("\t").append(diff));
	        double converged = varToDiff.get(var);
	        builder.append("\t").append(converged);
	        System.out.println(builder.toString());
	        builder.setLength(0);
	    });
	}
	
//	@Test
//	public void testPerturbationAnalysis() throws Exception {
//		FileUtility.initializeLogging();
//		Long dbId = 5693567L; // A sub-pathway in DNA repair
//		dbId = 400253L; // Circadian Clock
//	      // Circadian Clock
//        Map<String, Double> varToInhibition = new HashMap<>();
//		String varName = "NCOR1";
//        varToInhibition.put(varName,  0.99d);
//        varName = "HDAC3";
//        varToInhibition.put(varName, 0.99d);
//        
//        GKInstance pathway = getDBA().fetchInstance(dbId);
//		PathwayToBooleanNetworkConverter converter = new PathwayToBooleanNetworkConverter();
//		
//		PathwayImpactAnalysisResult result = performDrugImpactAnalysis(pathway, converter, drug, mapper)(pathway,
//		        converter,
//		        varToInhibition,
//		        null);
//		        
//		System.out.println("Perturbation results: " + result);
//	}

    private BooleanNetwork getNetwork(Long dbId) throws Exception {
        MySQLAdaptor dba = getDBA();
		
		GKInstance pathway = dba.fetchInstance(dbId);
		System.out.println("Working on " + pathway);
		
		PathwayToBooleanNetworkConverter converter = new PathwayToBooleanNetworkConverter();
		BooleanNetwork network = converter.convert(pathway);
        return network;
    }
	
    /**
     * Generate initial values for all BooleanVariables in the passed network.
     * @param network
     * @param defaultValue
     * @return
     */
	public Map<BooleanVariable, Number> createInitials(BooleanNetwork network,
	                                                   Number defaultValue) {
	    Map<BooleanVariable, Number> varToInial = new HashMap<>();
	    Set<BooleanVariable> cycleVars = BooleanNetworkUtilities.getVariablesInCycles(network);
	    network.getVariables().forEach(var -> varToInial.put(var, getInitial(var, 
	                                                                         cycleVars,
	                                                                         defaultValue)));
	    return varToInial;
	}
	
	private Number getInitial(BooleanVariable var,
	                          Set<BooleanVariable> cycleVars,
	                          Number defaultValue) {
	    if (cycleVars.contains(var)) {
	        // Check if it is a native variable, which is converted from original Reactome PE
	        String reactomeId = var.getProperty("reactomeId");
	        if (reactomeId != null)
	            return defaultValue; // Want to enable native variables in cycles first to avoid any problem. 
	                         // This may be overkill. However, currently there is no better way to figure
	                         // out what variables in cycles should be initialized first.
	                         // TODO: create a better algorithm to figure this out!
	    }
	    if (var.getInRelations() == null || var.getInRelations().size() == 0)
            return defaultValue;
        return 0.0d; // There is no need to activate them initially
    }
	
	/**
	 * Get a list of BooleanVariables that are converted from entities used as outputs.
	 * @param pathway
	 * @param varToDiff
	 * @return
	 * @throws Exception
	 */
	public List<BooleanVariable> getOutputVariables(Map<BooleanVariable, Double> varToDiff) {
        List<BooleanVariable> outputVars = varToDiff.keySet()
                                                    .stream()
                                                    .filter(var -> {
                                                        String output = var.getProperty("output");
                                                        if (output != null && output.equals(Boolean.TRUE + ""))
                                                            return true;
                                                        return false;
                                                    })
                                                    .collect(Collectors.toList());
        return outputVars;
	}
	
	/**
	 * Summarize the perturbation results for the whole pathway.
	 * @param varToDiff
	 * @param pathway
	 * @return
	 * @throws Exception
	 */
	public double summarizePerturbations(Map<BooleanVariable, Double> varToDiff) {
	    // Check the targets: We want to check outputs only
	    List<BooleanVariable> outputVars = getOutputVariables(varToDiff);
	    List<Double> values = new ArrayList<>();
	    outputVars.forEach(var -> {
	        Double diff = varToDiff.get(var);
	        if (diff != null)
	            values.add(diff);
	    });
	    if (values.size() == 0)
	        throw new IllegalArgumentException("The passed varToDiff doesn't contain output values!");
	    double sum = values.stream().mapToDouble(value -> Math.abs(value)).sum();
	    return sum / values.size();
	}
	
    private Map<BooleanVariable, String> getVarToGene(BooleanNetwork network) {
        Map<BooleanVariable, String> varToGene = new HashMap<>();
        network.getVariables().forEach(var -> {
            String gene = (String) var.getProperty("gene");
            if (gene != null)
                varToGene.put(var, gene);
        });
        return varToGene;
    }
	
	public PathwayImpactAnalysisResult performDrugImpactAnalysis(GKInstance pathway,
	                                                             PathwayToBooleanNetworkConverter converter,
	                                                             String drug,
	                                                             DrugToTargetsMapper mapper) throws Exception {
	    Set<String> targets = mapper.getDrugTargets(drug);
	    if (targets == null)
	        targets = new HashSet<>();
	    BooleanNetwork network = converter.convert(pathway, targets);
	    Map<BooleanVariable, String> varToGene = getVarToGene(network);
	    if (varToGene.size() == 0)
	        return null; // Nothing to do
	    // Most of drugs are used as inhibitors
	    Map<String, Double> geneToInhibition = mapper.getGeneToInhibition(varToGene.values(),
	                                                                      drug);
	    // However, some drugs may be used as activator.
	    Map<String, Double> geneToActivation = mapper.getGeneToActivation(varToGene.values(),
	                                                                      drug);

	    return performPerturbationAnalysis(pathway,
	            network, 
	            varToGene, 
	            geneToInhibition,
	            geneToActivation);
	}
	
	public PathwayImpactAnalysisResult performDrugHitAnalysis(GKInstance pathway,
	                                                          PathwayToBooleanNetworkConverter converter,
	                                                          String drug,
	                                                          DrugToTargetsMapper mapper) throws Exception {
	    Set<String> targets = mapper.getDrugTargets(drug);
	    if (targets == null)
	        targets = new HashSet<>();
	    BooleanNetwork network = converter.convert(pathway, targets);
	    Map<BooleanVariable, String> varToGene = getVarToGene(network);
	    if (varToGene.size() == 0)
	        return null; // Nothing to do
	    Set<String> sharedTargets = InteractionUtilities.getShared(targets, varToGene.values());
	    if (sharedTargets.size() == 0)
	        return null; // Nothing to report
	    PathwayImpactAnalysisResult result = new PathwayImpactAnalysisResult();
	    result.setPathwayName(pathway.getDisplayName());
	    result.setDbId(pathway.getDBID());
	    result.setTargetGenes(new ArrayList<>(sharedTargets));
	    return result;
	}
	
    private PathwayImpactAnalysisResult performPerturbationAnalysis(GKInstance pathway, 
                                                                    BooleanNetwork network,
                                                                    Map<BooleanVariable, String> varToGene, 
                                                                    Map<String, Double> geneToInhibition,
                                                                    Map<String, Double> geneToActivation) {
        // Convert gene to BooleanVariable
	    Map<BooleanVariable, Double> varToInhibition = new HashMap<>();
	    Map<BooleanVariable, Double> varToActivation = new HashMap<>();
	    BNPerturbationInjector injectors = new BNPerturbationInjector();
	    injectors.inject(network, geneToInhibition, geneToActivation, varToInhibition, varToActivation);
	    if (varToInhibition.size() == 0 && varToActivation.size() == 0)
	        return null; // There is no need to perform perturbation analysis

	    FuzzyLogicSimulator simulator = new FuzzyLogicSimulator();
	    simulator.setIteration(500); // Use 500 as the default to make sure an attractor can be reached actually.
	    //      simulator.enableDebug(true);
	    // Use the Identity function as the default here
	    simulator.setTransferFunction(new IdentityFunction());
	    SimulationConfiguration configuration = new SimulationConfiguration();
	    configuration.setDefaultValue(1.0d);
	    Map<BooleanVariable, Number> varToInitial = createInitials(network, 
	                                                               configuration.getDefaultValue());
	    configuration.setInitial(varToInitial);
	    // Background simulation
	    simulator.simulate(network, configuration);
	    SimulationResults reference = new SimulationResults();
	    reference.recordResults(network, simulator);

	    // Perform a perturbation analysis
	    configuration.setInhibition(varToInhibition);
	    configuration.setActivation(varToActivation);
	    
	    simulator.simulate(network, configuration);
	    SimulationResults drugSimulation = new SimulationResults();
	    drugSimulation.recordResults(network, simulator);

	    // Compare these two results
	    SimulationComparator comparator = new SimulationComparator();
	    String propKey = "reactomeId";
	    comparator.setVarPropertyKey(propKey);
	    // varToDiff is keyed based on reference variables
	    Map<BooleanVariable, Double> varToDiff = comparator.calculateDiff(drugSimulation, 
	            reference,
	            20, // These three parameters are arbitrary
	            0.001d,
	            1000);

	    List<BooleanVariable> outputVars = getOutputVariables(varToDiff);
	    double sum = 0.0d;
	    int counter = 0;
	    for (BooleanVariable var : outputVars) {
	        Double diff = varToDiff.get(var);
	        sum += Math.abs(diff);
	        counter ++;
	    }
	    if (counter == 0)
	        return null; 
	    // Want to list target genes too
	    Set<String> targets = new HashSet<>();
	    targets.addAll(geneToInhibition.keySet());
	    targets.addAll(geneToActivation.keySet());
	    List<String> geneList = new ArrayList<>(targets);
	    geneList.sort(Comparator.naturalOrder());
	    PathwayImpactAnalysisResult result = new PathwayImpactAnalysisResult();
	    result.setDbId(pathway.getDBID());
	    result.setPathwayName(pathway.getDisplayName());
	    result.setSum(sum);
	    result.setAverage(sum / counter);
	    result.setTargetGenes(geneList);
	    return result;
    }

}
