/*
 * Created on Mar 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

/**
 * Perform constrained Fuzzy logic simulation based on the implementation in CNOFuzzy.
 * @author gwu
 *
 */
public class FuzzyLogicSimulator extends Simulator {
    private boolean debug = false;
    private TransferFunction transferFunction;
    
    private double stopDiffValue = 1.0e-6;
    private int iteration = 100;
    // Cache values to check attractors
    private List<Number[]> values;
    private List<BooleanVariable> variables;
    // A flag to indicate if an attractor is reached
    private boolean isAttractorReached;
    // How to calculate AND. Default is minimum. However,
    // we may need to use product for AND converting from inputs
    // in reactions.
    private ANDGateMode andGateMode = ANDGateMode.MIN;
    
    /**
     * Default constructor.
     */
    public FuzzyLogicSimulator() {
    }
    
    public void setTransferFunction(TransferFunction function) {
        this.transferFunction = function;
    }
    
    public TransferFunction getTransferFunction() {
        return this.transferFunction;
    }
    
    public double getStopDiffValue() {
        return stopDiffValue;
    }

    public void setStopDiffValue(double stopDiffValue) {
        this.stopDiffValue = stopDiffValue;
    }

    public int getIteration() {
        return iteration;
    }

    public void setIteration(int iteration) {
        this.iteration = iteration;
    }

    public void enableDebug(boolean enabled) {
        this.debug = enabled;
    }
    
    public ANDGateMode getAndGateMode() {
        return andGateMode;
    }

    public void setAndGateMode(ANDGateMode andGateMode) {
        this.andGateMode = andGateMode;
    }

    /**
     * This method implements a synchronous update scheme using constrained fuzzy logic
     * framework as implemented in CNOFuzzy. The value of the stimulation should not be
     * changed during the update.
     */
    @Override
    public void simulate(BooleanNetwork network,
                         SimulationConfiguration configuration) {
        if (transferFunction == null)
            throw new IllegalStateException("TransferFunction has not specified!");
        Map<BooleanVariable, Number> stimulation = configuration.getStimulation();
        Map<BooleanVariable, Number> initial = configuration.getInitial();
        Map<BooleanVariable, Double> activation = configuration.getActivation();
        Map<BooleanVariable, Double> inhibition = configuration.getInhibition();
        for (BooleanVariable var : network.getVariables()) {
            Number value = stimulation.get(var);
            if (value != null)
                var.setValue(value);
            else if (initial.containsKey(var)) {
                var.setValue(initial.get(var));
            }
            else if (var.getValue() == null) {
                if (configuration.getDefaultValue() != null)
                    var.setValue(configuration.getDefaultValue());
                else
                    var.setValue(0.0d); // 0.0 will be used as the default value: all nodes are not activated!
            }
        }
        Map<BooleanVariable, Double> varToPostValue = new HashMap<>();
        // Cache calculated transfer values to avoid calculating them
        // repeated when nodes are involved in multiple relations.
        Map<BooleanVariable, Double> varToTransfer = new HashMap<>();
        // Start update
        int runCounter = 0;
        
        // For debugging
        List<BooleanVariable> vars = BooleanNetworkUtilities.getSortedVariables(network);
        this.variables = vars; // Keep it for returning
        values = new ArrayList<>();
        
        if (debug) {
            System.out.print("Iteration\t");
            for (BooleanVariable var : vars)
                System.out.print(var.getName() + "\t");
            System.out.println();
            System.out.print(runCounter + "\t");
            for (BooleanVariable var : vars)
                System.out.print(var.getValue() + "\t");
            System.out.println();
        }
        
        // Need to change input values first in case they are inhibited or activated
        modifyInputValues(network, inhibition, activation);
        // Keep the initial values as time = 0
        recordValues(vars);
        
        while (true) {
            runCounter ++;
//            System.out.println("Iteration: " + runCounter);
            varToTransfer.clear();
            for (BooleanVariable var : vars) {
                Set<BooleanRelation> inRelations = var.getInRelations();
                if (inRelations == null || inRelations.size() == 0)
                    continue;
                // Get the maximum value since OR date should be used
                double value = Double.NEGATIVE_INFINITY;
                boolean isTransfered = false;
                double negativeOnlyValue = 1.0d; // For negative only relation, we will dump 
                                                 // the calculated value or the intial value
                for (BooleanRelation rel : inRelations) {
                    Double transferValue = getTransferValue(rel, varToTransfer);
                    Map<BooleanVariable, Boolean> var2isNegated = rel.getInputVarToIsNegated();
                    if (var2isNegated.size() == 1 && var2isNegated.values().stream().findAny().get()) {
                    	// This should be negative value
                    	if (transferValue < negativeOnlyValue)
                    		negativeOnlyValue = transferValue;
                    }
                    else if (transferValue > value) {
                        value = transferValue;
                        isTransfered = true;
                    }
                }
                if (!isTransfered) {
                	value = var.getValue().doubleValue();
                }
                value *= negativeOnlyValue; // Apply dumping is any from negation
                // Apply inhibition and activation if any
                value = modifyValue(var, value, inhibition, activation);
                varToPostValue.put(var, value);
            }
            // Need to assign post value to variables
            double maxDiff = Double.NEGATIVE_INFINITY;
            for (BooleanVariable var : network.getVariables()) {
                // Don't change the stimulation values
                Number postValue = stimulation.get(var);
                if (postValue == null)
                    postValue = varToPostValue.get(var);
                else if (configuration.isUseLargerValueForStimulation() && varToPostValue.get(var) != null) {
                    if (varToPostValue.get(var) > postValue.doubleValue())
                        postValue = varToPostValue.get(var);
                }
                // Some variables are input that don't get updated so they will not have postValues
                if (postValue == null)
                    continue;
                double diff = Math.abs(var.getValue().doubleValue() - postValue.doubleValue());
                if (diff > maxDiff)
                    maxDiff = diff;
                var.setValue(postValue);
            }
            
            if (debug) {
                //System.out.println("MaxDiff: " + maxDiff);
                System.out.print(runCounter + "\t");
                for (BooleanVariable var : vars)
                    System.out.print(var.getValue() + "\t");
                System.out.println();
            }
            
            recordValues(vars);
            
            // The following three cases should stop this while loop
            if (_isAttractorReached()) {
                isAttractorReached = true;
                break;
            }
            // There is no need to do this
            if (maxDiff < stopDiffValue)
                break;
            if (runCounter > Math.max(iteration, network.getVariables().size() * 1.2d))
                break;
        }
        
        keepTracks();
    }
    
    private void modifyInputValues(BooleanNetwork network,
                                   Map<BooleanVariable, Double> inhibition,
                                   Map<BooleanVariable, Double> activation) {
        Map<BooleanVariable, Number> varToValue = new HashMap<>();
        network.getVariables().stream().forEach(var -> {
            if (var.isInputVariable()) {
                Number value = var.getValue();
                if (value == null)
                    return;
                value = modifyValue(var, 
                                    value.doubleValue(), 
                                    inhibition, 
                                    activation);
                var.setValue(value);
            }
        });
    }

    private double modifyValue(BooleanVariable var, 
                               double value, 
                               Map<BooleanVariable, Double> inhibition,
                               Map<BooleanVariable, Double> activation) {
        if (inhibition.containsKey(var)) // Inhibition is preferred
            value *= (1.0 - inhibition.get(var));
        else if (activation.containsKey(var)) {
            value *= (1.0 + activation.get(var));
        }
        if (value > 1.0d)
            value = 1.0d; // The maximum value
        return value;
    }

    private void keepTracks() {
        for (int i = 0; i < variables.size(); i++) {
            BooleanVariable var = variables.get(i);
            Number[] varValues = new Number[values.size()];
            for (int j = 0; j < values.size(); j++) {
                Number value = values.get(j)[i];
                varValues[j] = value;
            }
            var.setTrack(varValues);
        }
    }
    
    /**
     * This check is based on algorithm Floyd "tortoise and hare" cycle detection algorithm 
     * (https://en.wikipedia.org/wiki/Cycle_detection), and should be called for each step 
     * in stimulation.
     * @return
     */
    private boolean _isAttractorReached() {
        if (values.size() < 2)
            return false; // We need at least two data points
        // A simple case having one single value attractor
        int size = values.size();
        if (checkIsValuesEqual(values.get(size - 1), values.get(size - 2)))
            return true; // Single value attractor
        // After two values are added and at least have 3 values
        if (values.size() % 2 != 1) // Detect a loop based on Floyd algorithm
            return false;
        // Check last item and the middle one and see if they are the same
        // If it is the same, the attractor has been reached
        int middle = values.size() / 2; // e.g. 1 in 3 elements
        // Compare 1 and 2 in 3 elements or 2 and 4 in 5 elements
        return checkIsValuesEqual(values.get(middle), values.get(values.size() - 1));
    }
    
    /**
     * Check if there is any significant difference between two number arrays.
     * @param values1
     * @param values2
     * @return
     */
    private boolean checkIsValuesEqual(Number[] values1, Number[] values2) {
        for (int i = 0; i < values1.length; i++) {
            double diff = Math.abs(values1[i].doubleValue() - values2[i].doubleValue());
            if (diff > stopDiffValue)
                return false;
        }
        return true;
    }
    
    /**
     * Check if an attractor is reached after the simulation is finished.
     * This method should be called after {@link simulate(BooleanNetwork, Map<String, Number>) simulate}.
     * @return
     */
    public boolean isAttractorReached() {
        return this.isAttractorReached;
    }
    
    /**
     * Get the attractor after simulation. If an attractor has not been reached,
     * an exception will be thrown. This method should be called after
     * {@link simulate(BooleanNetwork, Map<String, Number>) simulate}.
     * @return
     */
    public Attractor getAttractor() {
        if (!isAttractorReached) {
            throw new IllegalStateException("An attractor is not reached!");
        }
        Attractor attractor = new Attractor();
        attractor.setVariables(variables);
        List<Number[]> attractorValues = new ArrayList<>();
        attractor.setValues(attractorValues);
        // Get the values from the end of the cached values
        for (int i = values.size() - 1; i >= 0; i --) {
            Number[] value = values.get(i);
            if (isContained(attractorValues, value))
                break;
            attractorValues.add(0, value);
        }
        return attractor;
    }
    
    private boolean isContained(List<Number[]> attractorValues,
                                Number[] values) {
        for (Number[] tmp : attractorValues) {
            if (checkIsValuesEqual(tmp, values))
                return true;
        }
        return false;
    }

    private void recordValues(List<BooleanVariable> variables) {
        Number[] currentValues = new Number[variables.size()];
        for (int i = 0; i < variables.size(); i++) 
            currentValues[i] = variables.get(i).getValue();
        values.add(currentValues);
    }
    
    private double getTransferValue(BooleanRelation rel,
                                    Map<BooleanVariable, Double> varToValue) {
        Map<BooleanVariable, Boolean> inputVarToIsNegated = rel.getInputVarToIsNegated();
        double rtn;
        if (andGateMode == ANDGateMode.PROD)
            rtn = 1.0d;
        else
            rtn = Double.MAX_VALUE;
        for (BooleanVariable var : inputVarToIsNegated.keySet()) {
            Boolean isNegated = inputVarToIsNegated.get(var);
            Double transferValue = getTransferValue(var, varToValue, rel.getTransferFunction());
            if (transferValue == null)
                continue;
            if (isNegated)
                transferValue = 1.0 - transferValue; 
            if (andGateMode == ANDGateMode.PROD)
                rtn *= transferValue;
            else { // Get the minimum value
                if (transferValue < rtn)
                    rtn = transferValue;
            }
        }
        return rtn;
    }
    
    private double getMaximumValue(BooleanRelation rel,
                                   Map<BooleanVariable, Double> varToValue) {
        Map<BooleanVariable, Boolean> inputVarToIsNegated = rel.getInputVarToIsNegated();
        // Get the minimum value
        double rtn = Double.MAX_VALUE;
        for (BooleanVariable var : inputVarToIsNegated.keySet()) {
            Boolean isNegated = inputVarToIsNegated.get(var);
            Double transferValue = getTransferValue(var, varToValue, rel.getTransferFunction());
            if (transferValue == null)
                continue;
            if (isNegated)
                transferValue = 1.0 - transferValue; 
            if (transferValue < rtn)
                rtn = transferValue;
        }
        return rtn;
    }
    
    private double getProductValue(BooleanRelation rel,
                                   Map<BooleanVariable, Double> varToValue) {
        Map<BooleanVariable, Boolean> inputVarToIsNegated = rel.getInputVarToIsNegated();
        // Get the minimum value
        double rtn = 1.0d;
        for (BooleanVariable var : inputVarToIsNegated.keySet()) {
            Boolean isNegated = inputVarToIsNegated.get(var);
            Double transferValue = getTransferValue(var, varToValue, rel.getTransferFunction());
            if (transferValue == null)
                continue;
            if (isNegated)
                transferValue = 1.0 - transferValue; 
            if (transferValue < rtn)
                rtn = transferValue;
        }
        return rtn;
    }
    
    private Double getTransferValue(BooleanVariable var,
                                    Map<BooleanVariable, Double> varToValue,
                                    TransferFunction transferFunction) {
    	if (transferFunction == null)
    		transferFunction = this.transferFunction; // Use the default TransferFunction is nothing here.
        Double value = varToValue.get(var);
        if (value == null) {
            if (var.isUseIdenticalTransfer())
                value = var.getValue() != null ? var.getValue().doubleValue() : null; 
            else
                value = transferFunction.transfer(var.getValue());
            if (value == null)
                return null; // Don't put null there
            varToValue.put(var, value);
        }
        return value;
    }

    @Test
    public void testSimulate() throws IOException {
        TransferFunction function = new HillFunction();
        function = new IdentityFunction();
        setTransferFunction(function);
        debug = true;
        BooleanNetworkReader reader = new BooleanNetworkReader();
        String fileName = "results/BooleanNetwork/ToyModel.txt";
        BooleanNetwork network = reader.readToyModel(fileName);
        Map<String, Number> varToStimulation = new HashMap<>();
        varToStimulation.put("EGF", 0.8d);
//        varToStimulation.put("TNFa", 1.0d);
//        varToStimulation.put("Ras", 1.0d);
        testSimulation(network, varToStimulation);
    }
    
    private void testSimulation(BooleanNetwork network,
                                Map<String, Number> varToStimulation) {
        simulate(network, varToStimulation);
//        for (BooleanVariable var : network.getVariables())
//            System.out.println(var.getName() + "\t" + var.getValue());
        if (isAttractorReached()) {
            System.out.println("Reached the following attractor: ");
            Attractor attractor = getAttractor();
            System.out.println(attractor.outputAsText());
        }
        else {
            System.out.println("Cannot reach an attractor!");
        }
    }
    
    @Test
    public void testSimpleFeedbackLoop() {
        BooleanNetwork network = BooleanNetworkUtilities.generateFeedbackLoopBN();
        TransferFunction function = new HillFunction();
//        function = new IdentityFunction();
        setTransferFunction(function);
        debug = true;
        Map<String, Number> varToStimulation = new HashMap<>();
        for (int i = 0; i < 10; i++) {
            // Reset
            network.getVariables().forEach(var -> var.setValue(0.0d));
            double stimulation = (i + 1) * 0.1;
            varToStimulation.put("A", stimulation);
            System.out.println("Stimulation A = " + stimulation);
            testSimulation(network, varToStimulation);
            System.out.println();
        }
    }
    
    /**
     * There are two ways to aggregate information from an AND gate: product or minimum of
     * all input values.
     * @author gwu
     *
     */
    public static enum ANDGateMode {
        MIN,
        PROD
    }
    
}
