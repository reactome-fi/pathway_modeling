/*
 * Created on Oct 27, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class implements a Gibbs sampling based inference. The actual implementation is ported
 * from the mallet java library (http://mallet.cs.umass.edu), and borrow some implementations 
 * from libdai (https://staff.fnwi.uva.nl/j.m.mooij/libDAI/).
 * The current implementation supports discrete variables only so all Observation objects used should
 * be parameterized by Integer.
 * @author gwu
 *
 */
public class GibbsSampling extends AbstractInferencer {
    private final int DEFAULT_MAX_ITERATION = 500;
    private static final Logger logger = Logger.getLogger(GibbsSampling.class);
    // Some Gibbs sampling specific attributes
    private int burnin;
    // The total iteration has been decided by the value of maxIteration.
    // Every restart will have its own burn-in.
    private int restart = 1; // Default is 1 means one restart only (without restart actually). 
    // A random object for randomization. Using the
    // Apache math package to get better performance.
    private RandomDataGenerator randomizer;
    // Cache calculate probabilitities
    private Map<Variable, Factor> varToMBFactor;
    // Cache variables and their neighbors including itself for fast access
    private Map<Variable, Set<Variable>> varToMbVars;
    
    /**
     * Default constructor.
     */
    public GibbsSampling() {
        randomizer = new RandomDataGenerator();
        maxIteration = DEFAULT_MAX_ITERATION;
    }
    
    public void setRandomGenerator(RandomDataGenerator randomizer) {
        this.randomizer = randomizer;
    }
    
    public RandomDataGenerator getRandomGenerator() {
        return this.randomizer;
    }
    
    public void setBurnin(int burnin) {
        this.burnin = burnin;
    }

    public int getBurnin() {
        return this.burnin;
    }
    
    public void setRestart(int restart) {
        this.restart = restart;
    }
    
    public int getRestart() {
        return this.restart;
    }

    @Override
    public synchronized void runInference() throws InferenceCannotConvergeException {
        super.runInference();
//        if (maxIteration < DEFAULT_MAX_ITERATION)
//            throw new IllegalStateException("MaxIteration should be greater than 500");
        resetCache();
        truncateContinuousFactors();
        attachObservation();
        // Get the list of samples
        List<List<Observation<Integer>>> allSamples = new ArrayList<List<Observation<Integer>>>();
        iteration = 0;
        for (int i = 0; i < restart; i++) {
            Observation<Integer> sample = initializeAssignment();
            // This sample is not the same one passed to the burn method
            sample = burn(sample);
            // Hold samples
            List<Observation<Integer>> samples = new ArrayList<Observation<Integer>>();
            sample(sample, samples);
            allSamples.add(samples);
        }
        double measure = calculateConvergence(allSamples);
        if (debug)
            logger.info("Measurement of Gibbs Sampling: " + measure);
        calculateMarginals(allSamples);
        detachObservation();
        addBackContinuosFactors();
    }

    protected void resetCache() {
        if (varToMBFactor == null)
            varToMBFactor = new HashMap<Variable, Factor>();
        else
            varToMBFactor.clear();
        if (varToMbVars == null)
            varToMbVars = new HashMap<Variable, Set<Variable>>();
        else
            varToMbVars.clear();
    }
    
    /**
     * This method implements a measure to check the quality of MCMC described in the
     * PGM book: page 523.
     * @param allSamples
     * @return
     */
    private double calculateConvergence(List<List<Observation<Integer>>> allSamples) {
        // Get the required average beliefs
        List<List<Double>> valuesList = new ArrayList<List<Double>>();
        int count = 0;
        // Pick up a Variable randomly
        Variable var = null;
        for (List<Observation<Integer>> samples : allSamples) {
            List<Double> values = new ArrayList<Double>();
            count = 0;
            for (int i = 0; i < samples.size(); i++) {
                Observation<Integer> sample = samples.get(i);
                if (var == null)
                    var = sample.getVariableToAssignment().keySet().iterator().next();
                Integer state = sample.getVariableToAssignment().get(var);
                if (state == 0)
                    count ++;
                if ((i + 1) % 100 == 0) {
                    values.add((double)count/(i + 1));
                }
            }
            valuesList.add(values);
        }
        // We use the counts of the first state in each sample
        List<Double> fkList = new ArrayList<Double>();
        double fkTotal = 0.0d;
        for (List<Double> values : valuesList) {
            double value = 0.0d;
            for (Double tmp : values) 
                value += tmp;
            fkList.add(value / values.size());
            fkTotal += fkList.get(fkList.size() - 1);
        }
        double f = fkTotal / valuesList.size();
        double B = 0.0d;
        for (int i = 0; i < fkList.size(); i++) {
            B += (fkList.get(i) - f) * (fkList.get(i) - f);
        }
        int M = valuesList.get(0).size();
        int k = valuesList.size();
        B *= M / (k == 1 ? 2 : k  - 1); // To avoid dividing by 0.
        double W = 0.0d;
        for (int i = 0; i < valuesList.size(); i++) {
            List<Double> values = valuesList.get(i);
            double fk = fkList.get(i);
            for (Double value : values) {
                W += (value - fk) * (value - fk);
            }
        }
        W /= (valuesList.size() * (M - 1));
        double V = (M - 1) / M * W + B / M;
        double R = Math.sqrt(V / W);
        return R;
    }
    
    /**
     * Initialize the first sample.
     */
    private Observation<Integer> initializeAssignment() {
        // The first assignment
        Observation<Integer> observation = new Observation<Integer>();
        Map<Variable, Integer> varToAssign = new HashMap<Variable, Integer>();
        for (Variable var : factorGraph.getVariables()) {
            int state = randomizer.nextInt(0, var.getStates() - 1);
            varToAssign.put(var, state);
        }
        observation.setVariableToAssignment(varToAssign);
        return observation;
    }
    
    /**
     * Perform the burnin phase.
     */
    protected Observation<Integer> burn(Observation<Integer> observation) {
        // Repeat sampling
        for (int i = 0; i < burnin; i++) {
            observation = sampleOnce(observation);
        }
        return observation;
    }
    
    /**
     * Generate samples after burn-in.
     */
    protected void sample(Observation<Integer> sample, List<Observation<Integer>> samples) {
        for (int i = 0; i < maxIteration; i++) {
            sample = sampleOnce(sample);
            sample.setName("Sample" + iteration);
            samples.add(sample);
            iteration++;
        }
    }
    
    /**
     * If the user set the restart, this maxIteration is used for each restart. The total iteration
     * should be got from getIteration(). At least 500 iterations is needed.
     */
    @Override
    public void setMaxIteration(int iteration) {
//        if (iteration < DEFAULT_MAX_ITERATION)
//            throw new IllegalArgumentException("Iteration should be greater than " + DEFAULT_MAX_ITERATION + ".");
        this.maxIteration = iteration;
    }
    
    /**
     * Calculate marginals for variables inside the FactorGraph object.
     */
    private void calculateMarginals(List<List<Observation<Integer>>> allSamples) {
        List<Observation<Integer>> samples = new ArrayList<Observation<Integer>>();
        for (List<Observation<Integer>> samples1 : allSamples)
            samples.addAll(samples1);
        // Marginals for variables
        // counts
        Map<Variable, int[]> varToCounts = new HashMap<Variable, int[]>();
        Map<Factor, int[]> factorToCounts = new HashMap<Factor, int[]>();
        for (Observation<Integer> sample : samples) {
            Map<Variable, Integer> varToAssgn = sample.getVariableToAssignment();
            countVariables(varToCounts, varToAssgn);
            countFactors(factorToCounts, varToAssgn);
        }
        boolean needNormalize = false;
        // Calculate beliefs for variables
        for (Variable var : varToCounts.keySet()) {
            needNormalize = false;
            double[] belief = new double[var.getStates()];
            int[] counts = varToCounts.get(var);
            for (int i = 0; i < counts.length; i++) {
                // To avoid belief is 0 that causes downstream analysis error,
                // 1 is used as a minimum.
                if (counts[i] == 0) {
                    belief[i] = 1.0d / samples.size();
                    needNormalize = true;
                }
                else
                    belief[i] = (double) counts[i] / samples.size();
            }
            if (needNormalize)
                var.normalize(belief, false, false);
            var.setBelief(belief);
        }
        // Calculate beliefs for factors
        for (Factor factor : factorToCounts.keySet()) {
            needNormalize = false;
            int[] counts = factorToCounts.get(factor);
            double[] belief = new double[counts.length];
            for (int i = 0; i < counts.length; i++) {
                if (counts[i] == 0) {
                    belief[i] = 1.0d / samples.size();
                    needNormalize = true;
                }
                else
                    belief[i] = (double) counts[i] / samples.size();
            }
            if (needNormalize)
                factor.normalize(belief, false, false);
            factor.setBelief(belief);
        }
    }
    
    
    
    private void countFactors(Map<Factor, int[]> factorToCounts,
                              Map<Variable, Integer> varToAssign) {
        Map<Variable, Integer> factorVarToAssign = new HashMap<Variable, Integer>();
        for (Factor factor : factorGraph.getFactors()) {
            int[] counts = factorToCounts.get(factor);
            if (counts == null) {
                counts = new int[factor.getValues().length]; // Get the number of state from the values
                factorToCounts.put(factor, counts);
            }
            factorVarToAssign.clear();
            for (Variable var : factor.getVariables()) {
                factorVarToAssign.put(var, varToAssign.get(var));
            }
            int index = factor.getIndexForAssignment(factorVarToAssign);
            counts[index] ++;
        }
    }

    private void countVariables(Map<Variable, int[]> varToCounts,
                               Map<Variable, Integer> varToAssgn) {
        for (Variable var : varToAssgn.keySet()) {
            int[] states = varToCounts.get(var);
            if (states == null) {
                states = new int[var.getStates()];
                varToCounts.put(var, states);
            }
            int assgn = varToAssgn.get(var);
            states[assgn] ++;
        }
    }

    /**
     * Perform one sampling.
     * @param obs
     */
    private Observation<Integer> sampleOnce(Observation<Integer> obs) {
        Observation<Integer> newSample = new Observation<Integer>();
        Map<Variable, Integer> obsVarToAssign = obs.getVariableToAssignment();
        Map<Variable, Integer> newVarToAssign = new HashMap<Variable, Integer>(obsVarToAssign);
        newSample.setVariableToAssignment(newVarToAssign);
        // Start sampling
        for (Variable var : factorGraph.getVariables()) {
            double[] probs = calculateProbabilities(var,
                                                    newVarToAssign);
            int assign = sampleOnce(probs);
            newVarToAssign.put(var, assign);
        }
        return newSample;
    }
    
    private int sampleOnce(double[] probs) {
        double sampled = randomizer.nextUniform(0.0d, 1.0d, true); // Include the lower point
        double cum = 0.0d;
        for (int i = 0; i < probs.length; i++) {
            cum += probs[i];
            if (sampled <= cum)
                return i;
        }
        // Otherwise, we need to throw an exception
        throw new IllegalStateException("Cannot find a legal state: " + Arrays.toString(probs));
    }
    
    private double[] calculateProbabilities(Variable var,
                                            Map<Variable, Integer> obsVarToAssign) {
        double[] rtn = getProbabilities(var, obsVarToAssign);
        if (rtn != null)
            return rtn;
        rtn = new double[var.getStates()];
        Map<Variable, Integer> assignment = new HashMap<Variable, Integer>();
        // Use log space for small values
        for (int i = 0; i < var.getStates(); i++) {
            double prob = 0.0d;
            for (Factor factor : var.getFactors()) {
                assignment.clear();
                // Other variables
                for (Variable var1 : factor.getVariables()) {
                    assignment.put(var1, obsVarToAssign.get(var1));
                }
                // This variable
                assignment.put(var, i);
                // Use log-space to avoid possible underflow.
                prob += Math.log(factor.getValue(assignment));
            }
            rtn[i] = Math.exp(prob);
        }
        // Normalize
        normalize(rtn);
        cacheProbabilities(rtn, var, obsVarToAssign);
        return rtn;
    }
    
    /**
     * Cache a calculated probability to save time.
     * @param var
     * @param varState
     * @param obsVarToAssgn
     */
    private void cacheProbabilities(double[] probs,
                                    Variable var,
                                    Map<Variable, Integer> obsVarToAssgn) {
        Set<Variable> mbVars = getMBVars(var);
        // If there are too many variables, don't cache it to comsume memory
        if (mbVars.size() > 6)
            return;
        // Check if the MB factor has been created
        Factor mbFactor = varToMBFactor.get(var);
        if (mbFactor == null) {
            mbFactor = new Factor();
            mbFactor.setVariables(new ArrayList<Variable>(mbVars));
            // Need to mark values as null
            double[] values = mbFactor.getValues();
            for (int i = 0; i < values.length; i++)
                values[i] = -1.0d; // Negative as not assigned
            varToMBFactor.put(var, mbFactor);
        }
        for (int i = 0; i< probs.length; i++) {
            obsVarToAssgn.put(var, i);
            mbFactor.setValue(probs[i], obsVarToAssgn);
        }
    }
    
    /**
     * Get a cache probability value to avoid re-calculation.
     * @param var
     * @param varState
     * @param obsVarToAssgn
     * @return
     */
    private double[] getProbabilities(Variable var,
                                      Map<Variable, Integer> obsVarToAssgn) {
        Factor mbFactor = varToMBFactor.get(var);
        if (mbFactor == null)
            return null;
//        Set<Variable> mbVars = getMBVars(var);
//        Map<Variable, Integer> mbAssgn = createMbAssgn(obsVarToAssgn, mbVars);
        double[] probs = new double[var.getStates()];
        for (int i = 0; i < var.getStates(); i++) {
//            mbAssgn.put(var, i);
            obsVarToAssgn.put(var, i);
            double prob = mbFactor.getValue(obsVarToAssgn);
            if (prob < 0.0d) // Not assigned yet
                return null;
            probs[i] = prob;
        }
        return probs;
    }

    private Set<Variable> getMBVars(Variable var) {
        Set<Variable> mbVars = varToMbVars.get(var);
        if (mbVars != null)
            return mbVars;
        // Query if there is a value for the combined state
        mbVars = new HashSet<Variable>();
        for (Factor factor : var.getFactors())
            mbVars.addAll(factor.getVariables());
        varToMbVars.put(var, mbVars);
        return mbVars;
    }
    
    private void normalize(double[] probs) {
        double sum = 0.0d;
        for (double d : probs)
            sum += d;
        for (int i = 0; i < probs.length; i++)
            probs[i] /= sum;
    }
    
    @Test
    public void testInitializeAssignment() {
        FactorGraph fg = TestUtilities.createSimpleFG();
        setFactorGraph(fg);
        Observation<Integer> obs = initializeAssignment();
        Map<Variable, Integer> varToAssign = obs.getVariableToAssignment();
        for (Variable var : varToAssign.keySet()) {
            System.out.println(var + ": " + varToAssign.get(var));
        }
    }
    
    @Test
    public void testInference() throws Exception {
        FileUtility.initializeLogging();
        FactorGraph fg = TestUtilities.createSimpleFG();
        System.out.println("FG is a tree: " + fg.isTree());
        
        Variable variable = TestUtilities.getVariable(fg, "protein");
        Map<Variable, Integer> observation = new HashMap<Variable, Integer>();
        observation.put(variable, 2);
        setObservation(observation);
        
        System.out.println("Inference using Gibbs sampling:");
        setDebug(true);
        setFactorGraph(fg);
        setBurnin(10000);
        setMaxIteration(10000);
        setRestart(5);
        long time1 = System.currentTimeMillis();
        runInference();
        long time2 = System.currentTimeMillis();
        System.out.println("Time: " + (time2 - time1));
        for (Variable var : fg.getVariables()) {
            System.out.println(var + ": " + Arrays.toString(var.getBelief()));
        }
        double logZ = calculateLogZ();
        System.out.println("LogZ: " + logZ);
        
        System.out.println("\nInfernece using LBP:");
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setFactorGraph(fg);
        lbp.setObservation(observation);
        lbp.setDebug(true);
        time1 = System.currentTimeMillis();
        lbp.runInference();
        time2 = System.currentTimeMillis();
        System.out.println("Time: " + (time2 - time1));
        for (Variable var : fg.getVariables()) {
            System.out.println(var + ": " + Arrays.toString(var.getBelief()));
        }
        logZ = calculateLogZ();
        System.out.println("LogZ: " + logZ);
    }
    
}
