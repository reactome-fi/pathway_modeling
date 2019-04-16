/*
 * Created on May 6, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;
import java.util.Map;

/**
 * This class is used to model a configuraton for performing Boolean network simulation. The configuration
 * may contain information about simulation, activation, inhibition, and initial values.
 * Note: The implementation of this class is based on Morris et al, Biotechnol. J. 2012, Section 2.3.
 * @author gwu
 *
 */
public class SimulationConfiguration {
    // Stimulation: Their values should be kept or larger values should be used based on flag (useLargerValuesForStimulation)
    private Map<BooleanVariable, Number> stimulation;
    private boolean useLargerValueForStimulation;
    // Activation: Simulated values should be increased by the specified values
    private Map<BooleanVariable, Double> activation;
    // Inhibition: Simulated values should be decreased by the specified values
    private Map<BooleanVariable, Double> inhibition;
    // Initial values: Variables listed in activation or inhibition may have intial values still
    private Map<BooleanVariable, Number> initial;
    // Default value in case no value is assigned for a Boolean variable
    private Number defaultValue;
    
    /**
     * Default constructor.
     */
    public SimulationConfiguration() {
    }

    public Number getDefaultValue() {
        return defaultValue;
    }

    public void setDefaultValue(Number defaultValue) {
        this.defaultValue = defaultValue;
    }

    public Map<BooleanVariable, Number> getStimulation() {
        if (stimulation == null)
            return new HashMap<>();
        return stimulation;
    }

    public void setStimulation(Map<BooleanVariable, Number> stimulation) {
        this.stimulation = stimulation;
    }

    public boolean isUseLargerValueForStimulation() {
        return useLargerValueForStimulation;
    }

    public void setUseLargerValueForStimulation(boolean useLargerValueForStimulation) {
        this.useLargerValueForStimulation = useLargerValueForStimulation;
    }

    public Map<BooleanVariable, Double> getActivation() {
        if (activation == null)
            return new HashMap<>();
        return activation;
    }

    public void setActivation(Map<BooleanVariable, Double> activation) {
        this.activation = activation;
    }

    public Map<BooleanVariable, Double> getInhibition() {
        if (inhibition == null)
            return new HashMap<>();
        return inhibition;
    }

    public void setInhibition(Map<BooleanVariable, Double> inhibition) {
        this.inhibition = inhibition;
    }

    public Map<BooleanVariable, Number> getInitial() {
        if (initial == null)
            return new HashMap<>();
        return initial;
    }

    public void setInitial(Map<BooleanVariable, Number> initial) {
        this.initial = initial;
    }
    
}
