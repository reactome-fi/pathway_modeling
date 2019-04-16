/*
 * Created on Jul 27, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

/**
 * This class is used to model simulation results generated from a Simulator.
 * @author wug
 *
 */
public class SimulationResults {
	
	private BooleanNetwork network;
	private Map<BooleanVariable, List<Number>> varToTrack;
	private Attractor attractor;
	
	public SimulationResults() {
	}

	public BooleanNetwork getNetwork() {
		return network;
	}

	public Map<BooleanVariable, List<Number>> getVarToTrack() {
		return this.varToTrack;
	}
	
	public Attractor getAttractor() {
		return attractor;
	}
	
	/**
	 * Copy this SimulationResults into another object. Attractor and BooleanNetwork
	 * in this SimulationResults are not copied. Only varToTrack is copied.
	 * @return
	 */
	public SimulationResults copy() {
	    SimulationResults copy = new SimulationResults();
	    copy.network = network;
	    copy.attractor = attractor;
	    if (varToTrack != null) {
	        Map<BooleanVariable, List<Number>> copyVarToTrack = new HashMap<>();
	        varToTrack.forEach((var, track) -> {
	            copyVarToTrack.computeIfAbsent(var, key -> new ArrayList<>(track));
	        });
	        copy.varToTrack = copyVarToTrack;
	    }
	    return copy;
	}
	
	/**
	 * Call this method to copy simulation results to this Object once the simulation is done
	 * to avoid simulation results overwritten.
	 * @param simulator
	 */
	public void recordResults(BooleanNetwork network,
						     Simulator simulator) {
		this.network = network;
		Attractor attractor = simulator.getAttractor();
		this.attractor = attractor;
		// Keep tracks
		varToTrack = new HashMap<>();
		network.getVariables().forEach(var -> {
			Number[] tracks = var.getTrack();
			// Make sure the list can be expanded
			List<Number> copy = new ArrayList<>();
			for (int i = 0; i < tracks.length; i++)
				copy.add(tracks[i]);
			varToTrack.put(var, copy);
		});
	}
	
	/**
	 * Calculate area under cover for each variable.
	 * @return
	 */
	public Map<BooleanVariable, Double> calculateAUC() {
		Map<BooleanVariable, Double> varToAUC = new HashMap<>();
		varToTrack.forEach((var, track) -> {
			double total = 0.0d;
			for (int i = 0; i < track.size(); i++) {
				if (i == 0 || i == track.size() - 1)
					total += track.get(i).doubleValue() / 2.0d;
				else
					total += track.get(i).doubleValue();
			}
			varToAUC.put(var,  total);
		});
		return varToAUC;
	}
	
	public int getTimeSteps() {
		List<Number> track = varToTrack.values().stream().findFirst().get();
		return track.size();
	}
	
	/**
	 * Expand the recorded tracks to more time steps as defined.
	 * @param maxTimeStep
	 */
	public void expandTracks(int maxTimeStep) {
		// Check if expand is needed
		// Get the first track
		int originalSize = getTimeSteps();
		if (originalSize >= maxTimeStep)
			return; // Do nothing is the passed maxTimeStep is less than the current size
		List<BooleanVariable> variables = attractor.getVariables();
		List<Number[]> values = attractor.getValues();
		int attractorSize = values.size();
		for (int i = 0; i < maxTimeStep - originalSize; i++) {
			int index = i % attractorSize;
			Number[] currentValues = values.get(index);
			for (int j = 0; j < variables.size(); j++) {
				BooleanVariable var = variables.get(j);
				List<Number> varTrack = varToTrack.get(var);
				varTrack.add(currentValues[j]);
			}
		}
	}
	
	@Test
	public void testExpandTracks() {
		BooleanNetwork network = BooleanNetworkUtilities.generateFeedbackLoopBN();
		Map<String, Number> varToStimulation = new HashMap<>();
		// There is a cycle track when A = 0.8
		double stimulation = 0.8d;
		varToStimulation.put("A", stimulation);
		System.out.println("Stimulation A = " + stimulation);
		FuzzyLogicSimulator simulator = new FuzzyLogicSimulator();
		simulator.simulate(network, varToStimulation);
		recordResults(network, simulator);
		
		System.out.println("\nAttactor: " + attractor.outputAsText());
		
		System.out.println("\nBefore expanding:");
		BooleanVariable var = network.getNameToVar().get("D");
		Map<BooleanVariable, Double> varToAUC = calculateAUC();
		checkVariable(var, varToAUC);
		
		System.out.println("\nExpand to 15 time steps:");
		expandTracks(15);
		varToAUC = calculateAUC();
		checkVariable(var, varToAUC);
	}
	
	private void checkVariable(BooleanVariable var,
						      Map<BooleanVariable, Double> varToAUC) {
		List<Number> track = varToTrack.get(var);
		int size = track.size();
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < size; i++)
			builder.append(i).append("\t");
		builder.append("Area");
		builder.append("\n");
		for (int i = 0; i < size; i++)
			builder.append(track.get(i)).append("\t");
		builder.append(varToAUC.get(var));
		System.out.println(builder.toString());
	}

}
