/*
 * Created on Jul 27, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

/**
 * This class is used to compare two simulation results.
 * @author wug
 *
 */
public class SimulationComparator {
    // Apparently using Double.MIN_VALUE can still generate Infinity
    // Use this hard-coded minimum value for comparison
    private final static double MIN_VALUE = 1.0E-8; // This is an arbitray value
    // Is this Object is used to compare values for two different BooleanNetwork,
    // The following propertyKey should be set to compare variables in two networks.
    private String varPropertyKey;
	
	public SimulationComparator() {
	}
	
	public String getVarPropertyKey() {
        return varPropertyKey;
    }

    public void setVarPropertyKey(String varPropertyKey) {
        this.varPropertyKey = varPropertyKey;
    }

    /**
	 * Use this method to calculate converged differences between two SimulationResults. 
	 * @param perturbed
	 * @param base Used as the reference
	 * @param increaseStep 
	 * @param threshold stop condition when the maximum change is less than this value
	 * @param maxStep maxStep
	 * @return
	 */
	public Map<BooleanVariable, Double> calculateDiff(SimulationResults perturbed,
	                                                  SimulationResults base,
	                                                  int increaseStep,
	                                                  double threshold,
	                                                  int maxStep) {
	    int start = Math.max(perturbed.getTimeSteps(), base.getTimeSteps());
	    Map<BooleanVariable, Double> preVarToDiff = null;
	    double maxChange = 1.0d; // This temp value will not be used actually
	    int currentStep = start;
	    while (true) {
	        Map<BooleanVariable, Double> varToDiff = calculateDiff(perturbed, base, currentStep);
	        if (preVarToDiff == null) 
	            preVarToDiff = varToDiff;
	        else {
	            maxChange = calculateMaxChange(varToDiff, preVarToDiff);
	            preVarToDiff = varToDiff;
	        }
	        currentStep += increaseStep;
	        // Check stop conditions based on total step and
	        // converging threshold
	        if (currentStep > maxStep || maxChange < threshold)
	            break;
	    }
	    return preVarToDiff;
	}
	
	private double calculateMaxChange(Map<BooleanVariable, Double> varToDiff,
	                                Map<BooleanVariable, Double> preVarToDiff) {
	    double maxChange = Double.NEGATIVE_INFINITY;
	    for (BooleanVariable var : varToDiff.keySet()) {
	        Double diff = varToDiff.get(var);
	        Double preDiff = preVarToDiff.get(var);
	        double change = Math.abs(diff - preDiff);
	        maxChange = Math.max(maxChange, change);
	    }
	    return maxChange;
	}
	
	/**
	 * Compare two SimulationResults.
	 * @param perturbed
	 * @param base
	 * @return
	 */
	public Map<BooleanVariable, Double> calculateDiff(SimulationResults perturbed,
											         SimulationResults base) {
		int timeSteps1 = perturbed.getTimeSteps();
		int timeSteps2 = base.getTimeSteps();
		int maxStep = Math.max(timeSteps1, timeSteps2);
		return calculateDiff(perturbed, base, maxStep);
	}
	
	/**
	 * Compare two results after expanding to the same max time steps specified by the parameter.
	 * @param perturbed
	 * @param base
	 * @param maxTimeStep
	 * @return
	 */
	public Map<BooleanVariable, Double> calculateDiff(SimulationResults perturbed,
											   		 SimulationResults base,
											   		 int maxTimeStep) {
	    // BooleanNetworks can be compared (e.g. minor modification in a BooleanNetwork).
	    // Make sure two Results are generated from the same network
	    if (perturbed.getNetwork() != base.getNetwork() && varPropertyKey == null)
	        throw new IllegalArgumentException("The passed two parameters don't have the same BooleanNetwork! "
	                + "To compare two different BooleanNetworks, set varPropertyKey.");
		perturbed.expandTracks(maxTimeStep);
		base.expandTracks(maxTimeStep);
		// Calculate AUC
		Map<BooleanVariable, Double> pertVarToAuc = perturbed.calculateAUC();
		Map<BooleanVariable, Double> baseVarToAuc = base.calculateAUC();
		Map<String, Double> pertPropToAuc = null;
		if (varPropertyKey != null) {
		    pertPropToAuc = convertVarToAucMap(pertVarToAuc, varPropertyKey);
		}
		final Map<String, Double> finalPertPropToAuc = pertPropToAuc;
		Map<BooleanVariable, Double> varToDiff = new HashMap<>();
		// Get the percent change
		baseVarToAuc.forEach((var, baseAuc) -> {
		    Double pertAuc = null;
		    if (varPropertyKey != null) {
		        String property = var.getProperty(varPropertyKey);
		        pertAuc = finalPertPropToAuc.get(property);
		    }
		    else
		        pertAuc = pertVarToAuc.get(var);
			if (pertAuc == null)
			    return; 
			// Assign the minimum value to avoid NaN
			if (baseAuc < MIN_VALUE)
			    baseAuc = MIN_VALUE;
			if (pertAuc < MIN_VALUE)
			    pertAuc = MIN_VALUE;
//			double change = Math.abs(pertAuc - baseAuc) / baseAuc;
			double change = (pertAuc - baseAuc) / (pertAuc + baseAuc); // Use this formula to control small values
			varToDiff.put(var,  change);
		});
		
		return varToDiff;
	}
	
	private Map<String, Double> convertVarToAucMap(Map<BooleanVariable, Double> varToAuc,
	                                               String propertyKey) {
	    Map<String, Double> propertyToAuc = new HashMap<>();
	    varToAuc.forEach((var, auc) -> {
	        String prop = var.getProperty(propertyKey);
	        if (prop == null)
	            return;
	        propertyToAuc.put(prop, auc);
	    });
	    return propertyToAuc;
	}
	
	@Test
	public void testCalculateDiff() {
		FuzzyLogicSimulator simulator = new FuzzyLogicSimulator();
		BooleanNetwork network = BooleanNetworkUtilities.generateFeedbackLoopBN();
		Map<String, Number> varToStimulation = new HashMap<>();
		
		// There is a cycle track when A = 0.8
		double stimulation = 0.8d;
		varToStimulation.put("A", stimulation);
		System.out.println("Stimulation A = " + stimulation);
		
		simulator.simulate(network, varToStimulation);
		SimulationResults base = new SimulationResults();
		base.recordResults(network, simulator);
		System.out.println("Total timesteps: " + base.getTimeSteps());
		
		stimulation = 0.6d;
		varToStimulation.put("A", stimulation);
		System.out.println("Stimulation A = " + stimulation);
		
		simulator.simulate(network, varToStimulation);
		SimulationResults perturbed = new SimulationResults();
		perturbed.recordResults(network, simulator);
		System.out.println("Total timesteps: " + perturbed.getTimeSteps());
		
		System.out.println("Calculate diff with default timesteps:");
		SimulationComparator comparator = new SimulationComparator();
		Map<BooleanVariable, Double> varToDiff = comparator.calculateDiff(perturbed, base);
		varToDiff.forEach((var, diff) -> {
			System.out.println(var + "\t" + diff);
		});
		
		// Check different time steps
		int[] timeSteps = {20, 40, 50, 100, 200, 300};
		for (int step : timeSteps) {
			System.out.println("Timestep: " + step);
			varToDiff = comparator.calculateDiff(perturbed,
										        base,
										        step);
			varToDiff.forEach((var, diff) -> {
				System.out.println(var + "\t" + diff);
			});
		}
	}

}
