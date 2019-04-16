/*
 * Created on Oct 27, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * An abstract Inferencer to provide some common attributes and methods.
 * @author gwu
 *
 */
public abstract class AbstractInferencer implements Inferencer {
    // FactorGraph upon which this LBP is performed
    protected FactorGraph factorGraph;
    // Observation
    private Observation<? extends Number> observation;
    private List<Factor> observationFactors;
    protected double tolerance = 1.0e-6;
    protected int maxIteration = 10000;
    protected int iteration;
    protected double maxDiff;
    protected boolean debug;
    // These two properties are related to CLGFactors and CLGVariables
    private List<Factor> factorsFromContinuous; // Marginalized from CLGFactors
    private List<ContinuousFactor> continuousFactors; // CLGFactors in the original FactorGraph
    
    protected AbstractInferencer() {
    }
    
    /* (non-Javadoc)
     * @see org.reactome.factorgraph.Inferencer#setFactorGraph(org.reactome.factorgraph.FactorGraph)
     */
    @Override
    public void setFactorGraph(FactorGraph fg) {
        this.factorGraph = fg;
    }
    
    /* (non-Javadoc)
     * @see org.reactome.factorgraph.Inferencer#getFactorGraph()
     */
    @Override
    public FactorGraph getFactorGraph() {
        return this.factorGraph;
    }
    
    /* (non-Javadoc)
     * @see org.reactome.factorgraph.Inferencer#setObservation(java.util.Map)
     */
    @Override
    public <T extends Number> void setObservation(Map<Variable, T> varToAssign) {
        if (varToAssign == null) {
            clearObservation();
            return;
        }
        // Create an observation
        Observation<T> obs = new Observation<T>();
        obs.setVariableToAssignment(varToAssign);
        this.observation = obs; // Have to go through this intermediate step. Otherwise, the above assignment cannot work.
    }
    
    @Override
    public <T extends Number> void setObservation(Observation<T> observation) {
        if (observation == null) {
            clearObservation();
            return;
        }
        this.observation = observation;
    }
    
    @Override
    public void clearObservation() {
        this.observation = null;
    }

    @Override
    public Observation<? extends Number> getObservation() {
        return this.observation;
    }
    
    public double getTolerance() {
        return tolerance;
    }

    public void setTolerance(double tolerance) {
        this.tolerance = tolerance;
    }

    public int getMaxIteration() {
        return maxIteration;
    }

    public void setMaxIteration(int maxIteration) {
        this.maxIteration = maxIteration;
    }

    /**
     * A read-only property.
     * @return
     */
    public int getIteration() {
        return iteration;
    }

    public boolean getDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }
    
    public double getMaxDiff() {
        return this.maxDiff;
    }
    
//  private double calculateMaxDiff(FactorGraph fg) {
//  double maxDiff = 0.0d;
//  Set<FGNode> nodes = new HashSet<FGNode>(fg.getVariables());
//  nodes.addAll(fg.getFactors());
//  for (FGNode node : nodes) {
//      double[] belief = node.getBelief();
//      double[] copy = new double[belief.length];
//      System.arraycopy(belief, 0, copy, 0, belief.length);
//      node.updateBelief(logSpace);
//      belief = node.getBelief();
//      for (int i = 0; i < copy.length; i++) {
//          double diff = Math.abs(copy[i] - belief[i]);
//          if (diff > maxDiff)
//              maxDiff = diff;
//      }
//  }
//  return maxDiff;
//}

    /* (non-Javadoc)
     * @see org.reactome.factorgraph.Inferencer#runInference()
     */
    @Override
    public void runInference() throws InferenceCannotConvergeException {
        if (factorGraph == null)
            throw new IllegalArgumentException("The target FactorGraph has not been assigned.");
        if (!factorGraph.isInferreable())
            throw new IllegalArgumentException("This type of FactorGraph is not supported yet: probably continuous variables are not leaf nodes.");
    }
    
    /**
     * Remove CLGNodes so that we can use inference algorithms implemented for discrete variables only.
     */
    protected void truncateContinuousFactors() {
        if (continuousFactors == null)
            continuousFactors = new ArrayList<ContinuousFactor>();
        else
            continuousFactors.clear();
        if (factorsFromContinuous == null)
            factorsFromContinuous = new ArrayList<Factor>();
        else
            factorsFromContinuous.clear();
        // Go through all factors
        for (Factor factor : factorGraph.getFactors()) {
            if (!(factor instanceof ContinuousFactor))
                continue;
            ContinuousFactor contFactor = (ContinuousFactor) factor;
            continuousFactors.add(contFactor);
            Factor convertedFactor = convertContinuousFactor(contFactor);
            factorsFromContinuous.add(convertedFactor);
        }
        factorGraph.getFactors().removeAll(continuousFactors);
        for (Factor factor : factorsFromContinuous)
            factorGraph.addFactor(factor);
        factorGraph.validatVariables(); // Some variables are not included
    }
    
    protected void addBackContinuosFactors() {
        if (continuousFactors.size() == 0 && factorsFromContinuous.size() == 0)
            return; // Nothing to do
        // Remove converted factors
        detachFactors(factorsFromContinuous);
        // Add back CLGFactors
        for (ContinuousFactor factor : continuousFactors) {
            factorGraph.addFactor(factor);
        }
        factorGraph.validatVariables();
    }
    
    private Factor convertContinuousFactor(ContinuousFactor continuousFactor) {
        ContinuousVariable continuousVar = continuousFactor.getContinuousVariable();
        Variable disVar = continuousFactor.getDiscreteVariable();
        Factor factor = new Factor();
        List<Variable> varList = new ArrayList<Variable>();
        varList.add(disVar);
        factor.setVariables(varList);
        disVar.addFactor(factor);
        // Factor value
        VariableAssignment<? extends Number> assignment = null;
        if (observation != null)
            assignment = observation.getVariableAssignment(continuousVar);
        double[] factorValues = continuousFactor.marginalizeForDiscrete(assignment);
        if (Double.isNaN(factorValues[0]))
            System.out.println("Wrong factor values!");
        factor.setValues(factorValues);
        // Modify it a little bit
        disVar.removeFactor(continuousFactor);
        return factor;
    }

    /**
     * This method is used to calculate the Bethe approximation of the partition function
     * for the passed FactorGraph. The client should call runInference() first for the
     * passed FactorGraph object. Otherwise, the result will not be correct. The actual
     * calculation is based on this paper: Graph polynomials and approximation of partition 
     * functions with Loopy Belief Propagation by Watanabe & Fukumizu (http://arxiv.org/pdf/0903.4527.pdf).
     * @return
     */
    public double calculateLogZ() {
        double logZ = 0.0d;
        double tmp = 0.0d; // Just a tmp variable
        for (Variable variable : factorGraph.getVariables()) {
            tmp = 0.0d;
            for (double value : variable.getBelief()) {
                if (value == 0.0d)
                    continue;
                tmp += value * Math.log(value);
            }
            logZ += (variable.getFactors().size() - 1) * tmp;
            if (Double.isNaN(logZ)) {
                throw new IllegalStateException("NaN encountered in variable: " + variable + " with belief " + variable.getBelief());
            }
        }
        for (Factor factor : factorGraph.getFactors()) {
            double[] belief = factor.getBelief();
            double[] values = factor.getValues();
            for (int i = 0; i < belief.length; i++) {
                if (belief[i] == 0.0d || values[i] == 0.0d)
                    continue; // Escape these cases
                // In case belief[i] is extremely small, using Math.log(values[i]/belief[i])
                // may generate an infinity double. So have to split into two logarithm calculation.
                logZ += belief[i] * (Math.log(values[i]) - Math.log(belief[i]));
                if (Double.isNaN(logZ)) {
                    throw new IllegalStateException("NaN encountered in factor: " + factor + " with belief " + belief[i]);
                }
            }
        }
        return logZ;
    }
    
    protected void attachObservation() {
        if (observation == null)
            return;
        if (observationFactors == null)
            observationFactors = new ArrayList<Factor>();
        else
            observationFactors.clear();
        Set<Variable> variables = factorGraph.getVariables();
        Map<Variable, ? extends Number> varToValue = observation.getVariableToAssignment();
        for (Variable var : varToValue.keySet()) {
            if (!variables.contains(var)) // At this stage, CLGVariable should not be here any more
                continue;
            // Add a new factor
            Factor factor = new Factor();
            List<Variable> list = new ArrayList<Variable>();
            list.add(var);
            factor.setVariables(list);
            var.addFactor(factor);
            double[] values = new double[var.getStates()];
            if (var.getStates() == 2) { // Two states for perturbation 
                // The value should be the strength. The largest should be 1.0
                double value = varToValue.get(var).doubleValue();
                if (value > 1.0d || value < 0.0d)
                    throw new IllegalStateException(var + " value should be in [0, 1]. The assigne value is: " + value);
                values[0] = 1.0d - value;
                values[1] = value;
            }
            else // More than two states
                values[varToValue.get(var).intValue()] = 1.0d;
            factor.setValues(values);
            factorGraph.addFactor(factor);
            observationFactors.add(factor);
        }
    }

    protected void detachObservation() {
        if (observation == null)
            return;
        detachFactors(observationFactors);
    }

    private void detachFactors(List<Factor> factors) {
        // Remove factors added for the observation.
        factorGraph.getFactors().removeAll(factors);
        // Remove edges in variables related to these factors
        for (Factor factor : factors) {
            detachFactor(factor);
        }
        factors.clear();
    }

    private void detachFactor(Factor factor) {
        // Need to clean-up the data structure a little bit.
        for (Variable var : factor.getVariables())
            var.removeFactor(factor);
        if (factor.getOutEdges() != null) {
            for (Edge edge : factor.getOutEdges()) {
                Variable variable = (Variable) edge.getToNode();
                variable.removeInEdge(edge);
                variable.removeOutEdge(factor);
            }
        }
    }
    
}
