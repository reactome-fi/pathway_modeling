/*
 * Created on Mar 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;
import java.util.Map;

/**
 * An abstract class for performing Boolean network simulation.
 * @author gwu
 *
 */
public abstract class Simulator {
    
    public Simulator() {
    }
    
    /**
     * Perform simulation for the passed BooleanNetwork object with the values specified
     * in the passed name to value map.
     * @param network
     */
    public void simulate(BooleanNetwork network,
                         Map<String, Number> stimulation) {
        SimulationConfiguration configutation = new SimulationConfiguration();
        Map<String, BooleanVariable> nameToVar = network.getNameToVar();
        Map<BooleanVariable, Number> varStimulation = new HashMap<>();
        for (String name : stimulation.keySet()) {
            BooleanVariable var = nameToVar.get(name);
            if (var == null)
                continue;
            varStimulation.put(var, stimulation.get(name));
        }
        configutation.setStimulation(varStimulation);
        simulate(network, configutation);
    }
    
    /**
     * Perform simulation for the passed BooleanNetwork based on a complex configuration
     * specified by a SimulationConfiguration object.
     * @param network
     * @param configuration
     */
    public abstract void simulate(BooleanNetwork network,
                                  SimulationConfiguration configuration);
    
    /**
     * Get the attractor after the simulation is done. This method should be called
     * after a simulation has been performed. Otherwise, a null should be returned.
     * @return
     */
    public abstract Attractor getAttractor();
    
}
