/*
 * Created on Jun 13, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.log4j.PropertyConfigurator;
import org.junit.Test;
import org.reactome.factorgraph.CLGFactor.CLGFactorDistribution;
import org.reactome.factorgraph.ContinuousVariable.DistributionType;
import org.reactome.r3.util.FileUtility;

/**
 * @author gwu
 *
 */
public class LBPTester {
    private LoopyBeliefPropagation lbp;
    
    public LBPTester() {
        PropertyConfigurator.configure("resources/log4j.properties");
        lbp = new LoopyBeliefPropagation();
    }
    
    @Test
    public void testMaxProductInference() throws InferenceCannotConvergeException {
        PropertyConfigurator.configure("resources/log4j.properties");
        // Example 13.1 in the PGM book
        FactorGraph fg = new FactorGraph();
        
        Variable a = new Variable(2);
        a.setName("A");
        Variable b = new Variable(2);
        b.setName("B");
        Factor factor = new Factor();
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(a);
        double[] values = new double[2];
        values[0] = 0.4d;
        values[1] = 0.6d;
        factor.setVariables(variables);
        factor.setValues(values);
        fg.addFactor(factor);
        
        variables = new ArrayList<Variable>();
        variables.add(a);
        variables.add(b);
        factor = new Factor();
        factor.setVariables(variables);
        values = new double[4];
        values[0] = 0.1;
        values[1] = 0.55;
        values[2] = 0.9;
        values[3] = 0.45;
        factor.setValues(values);
        fg.addFactor(factor);
        
        fg.validatVariables();
        System.out.println("Example 13.1:");
        testMaxProd(fg);
        
        // Example 13.10
//        values = new double[4];
//        values[0] = 0.1;
//        values[1] = 0.4;
//        values[2] = 0.4;
//        values[3] = 0.1;
        values = new double[]{0.3, 0.4, 0.3, 0.0}; // Bishop book: page 411.
        factor.setValues(values);
        fg = new FactorGraph();
        fg.addFactor(factor);
        fg.validatVariables();
        System.out.println("\nExample 13.10:");
        testMaxProd(fg);
        
        // Example 13.11: A frustrated loops
        fg = new FactorGraph();
        factor = new Factor();
        variables = new ArrayList<Variable>();
        variables.add(a);
        variables.add(b);
        String valueText = "1, 2, 2, 1";
        factor.setVariables(variables);
        factor.setValues(parseValues(valueText));
        fg.addFactor(factor);
        
        Variable c = new Variable(2);
        c.setName("C");
        factor = new Factor();
        variables = new ArrayList<Variable>();
        variables.add(b);
        variables.add(c);
        factor.setVariables(variables);
        factor.setValues(parseValues(valueText));
        fg.addFactor(factor);
        
        factor = new Factor();
        variables = new ArrayList<Variable>();
        variables.add(a);
        variables.add(c);
        factor.setVariables(variables);
        factor.setValues(parseValues(valueText));
        fg.addFactor(factor);
        fg.validatVariables();
        System.out.println("\nExample 13.11:");
        testMaxProd(fg);
        
    }
    
    private double[] parseValues(String values) {
        String[] tokens = values.split(", ");
        double[] rtn = new double[tokens.length];
        for (int i = 0; i < rtn.length; i++)
            rtn[i] = new Double(tokens[i]);
        return rtn;
    }

    private void testMaxProd(FactorGraph fg) throws InferenceCannotConvergeException {
        lbp.setFactorGraph(fg);
        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.runInference();
        System.out.println("In the probability space:");
        Map<Variable, Integer> varToAssign = lbp.findMaximum();
        for (Variable var : varToAssign.keySet())
            System.out.println(var.getName() + ": " + varToAssign.get(var));
        lbp.setUseLogSpace(true);
        lbp.runInference();
        System.out.println("\nIn the log space:");
        varToAssign = lbp.findMaximum();
        for (Variable var : varToAssign.keySet())
            System.out.println(var.getName() + ": " + varToAssign.get(var));
        System.out.println("LogLikelihood: " + fg.getLogLikelihood(varToAssign));
        
        System.out.println("\nRun SUM_PRODUCT:");
        lbp.setFactorGraph(fg);
        lbp.setInferenceType(InferenceType.SUM_PRODUCT);
        lbp.runInference();
        outputBelief(fg);
    }
    
    @Test
    public void testEmpiricalDistribution() {
        // Use normal distribution as an example
        int count = 1000;
        NormalDistribution normal = new NormalDistribution(0.0d, 1.0d);
        double[] values = new double[count];
        for (int i = 0; i < count; i++) {
            values[i] = normal.sample();
//            System.out.println(values[i]);
        }
        EmpiricalDistribution empirical = new EmpiricalDistribution();
        empirical.load(values);
        System.out.println("\n");
        // Try these values, which have to be listed in the values array.
        // Otherwise, NaN will be returned.
        Arrays.sort(values);
        int[] indices = new int[] {
                (int) (count * 0.05d),
                (int) (count * 0.45d),
                (int) (count * 0.50d),
                (int) (count * 0.55d),
                (int) (count * 1.0d - 1)
        };
        for (int index : indices) {
            double xValue = values[index];
            System.out.println(xValue + ": " + empirical.cumulativeProbability(xValue));
        }
    }
    
    @Test
    public void testRunEmpiricalInference() throws InferenceCannotConvergeException {
        FileUtility.initializeLogging();
        // A -> B -> X: A and B are discrete and X is continuous variables
        Set<Factor> factors = new HashSet<Factor>();
        Factor factorAB = createSimpleABFactor();
        factors.add(factorAB);
        
        Variable a = getVariable(factorAB, "A");
        Variable b = getVariable(factorAB, "B");
        
        ContinuousVariable x = new ContinuousVariable();
        x.setName("X");
        x.setDistributionType(DistributionType.TWO_SIDED);
        EmpiricalFactor factorBX = new EmpiricalFactor();
        factorBX.setDiscreteVariable(b);
        factorBX.setContinuousVariable(x);
        factors.add(factorBX);
        
        FactorGraph fg = new FactorGraph();
        fg.setFactors(factors);
        fg.validatVariables();
        
        // Generate an EmpiricalDistribution
        int count = 1000;
        // Use normal distribution as an example
        NormalDistribution normal = new NormalDistribution(0.0d, 1.0d);
        double[] values = new double[count];
        for (int i = 0; i < count; i++)
            values[i] = normal.sample();
        EmpiricalDistribution empirical = new EmpiricalDistribution();
        empirical.load(values);
        
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setFactorGraph(fg);
        
        // Try these values
        // Try these values, which have to be listed in the values array.
        // Otherwise, NaN will be returned.
        Arrays.sort(values);
        int[] indices = new int[] {
                0,
                (int) (count * 0.05d),
                (int) (count * 0.25d),
                (int) (count * 0.50d),
                (int) (count * 0.75d),
                (int) (count * 1.0d - 1)
        };
        
        // Run three times to ensure the results are the same
        for (int k = 0; k < 3; k++) {
            lbp.clearObservation();
            lbp.runInference();
            double[] aValues = a.getBelief();
            double[] bValues = b.getBelief();
            System.out.println("Belief for A:");
            for (int i = 0; i < aValues.length; i++)
                System.out.println("State " + i + ": " + aValues[i]);
            System.out.println("Belief for B:");
            for (int i = 0; i < bValues.length; i++)
                System.out.println("State " + i + ": " + bValues[i]);
            VariableAssignment<Double> varAssgn = new VariableAssignment<Double>();
            varAssgn.setVariable(x);
            varAssgn.setDistribution(empirical);
            Observation<Double> observation = new Observation<Double>();
            observation.addAssignment(varAssgn);
            lbp.setObservation(observation);
            for (int index : indices) {
                double xValue = values[index];
                varAssgn.setAssignment(xValue);
                lbp.runInference();
                System.out.println("X = " + xValue);
                System.out.println("Belief for A:");
                for (int i = 0; i < aValues.length; i++)
                    System.out.println("State " + i + ": " + aValues[i]);
                System.out.println("Belief for B:");
                for (int i = 0; i < bValues.length; i++)
                    System.out.println("State " + i + ": " + bValues[i]);
                // Just the distribution for the continuous value
                System.out.println("Continuous probability: " + empirical.cumulativeProbability(xValue));
            }
            System.out.println();
        }
        
    }
    
    private Factor createSimpleABFactor() {
        Variable a = new Variable(2);
        a.setName("A");
        Variable b = new Variable(2);
        b.setName("B");
        Factor factorAB = new Factor();
        List<Variable> list = new ArrayList<Variable>();
        list.add(a);
        list.add(b);
        factorAB.setVariables(list);
        double[] values = new double[]{0.8, 0.6, 0.2, 0.4};
        factorAB.setValues(values);
        return factorAB;
    }
    
    private Variable getVariable(Factor factor,
                                 String name) {
        for (Variable var : factor.getVariables())
            if (var.getName().equals(name))
                return var;
        return null;
    }
    
    @Test
    public void testRunGaussianInference() throws InferenceCannotConvergeException {
        FileUtility.initializeLogging();
        // A -> B -> X: A and B are discrete and X is continuous variables
        Set<Factor> factors = new HashSet<Factor>();
        Factor factorAB = createSimpleABFactor();
        factors.add(factorAB);
        
        Variable a = getVariable(factorAB, "A");
        Variable b = getVariable(factorAB, "B");
        
        Map<Variable, Integer> varToState = new HashMap<Variable, Integer>();
        varToState.put(a, 1);
        varToState.put(b, 0);
        System.out.println("Index for " + varToState + ": " + factorAB.getIndexForAssignment(varToState));
        
//        // Test code
//        Factor factorB = new Factor();
//        List<Variable> list = new ArrayList<Variable>(1);
//        list.add(b);
//        factorB.setVariables(list);
//        factorB.setValues(new double[]{0.119, 0.881});
//        factors.add(factorB);
//        FactorGraph fg = new FactorGraph();
//        fg.setFactors(factors);
//        fg.validatVariables();
//        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
//        lbp.setDebug(true);
//        lbp.setFactorGraph(fg);
//        lbp.runInference();
//        double[] aValues = a.getBelief();
//        double[] bValues = b.getBelief();
//        System.out.println("Belief for A:");
//        for (int i = 0; i < aValues.length; i++)
//            System.out.println("State " + i + ": " + aValues[i]);
//        System.out.println("Belief for B:");
//        for (int i = 0; i < bValues.length; i++)
//            System.out.println("State " + i + ": " + bValues[i]);
//        
        ContinuousVariable x = new ContinuousVariable();
        x.setName("X");
        CLGFactor factorBX = new CLGFactor();
        factorBX.setDiscreteVariable(b);
        factorBX.setContinuousVariable(x);
        
        List<CLGFactorDistribution> distributions = new ArrayList<CLGFactor.CLGFactorDistribution>();
        CLGFactorDistribution factorDistribution0 = new CLGFactorDistribution(0, 
                                                                              0.5, 
                                                                              new NormalDistribution(0.0d, 1.0d));
        CLGFactorDistribution factorDistribution1 = new CLGFactorDistribution(1,
                                                                              0.5, 
                                                                              new NormalDistribution(2.0d, 1.0d));
        distributions.add(factorDistribution0);
        distributions.add(factorDistribution1);
        
        factorBX.setDistributions(distributions);
        
        factors.add(factorBX);
        
        
        FactorGraph fg = new FactorGraph();
        fg.setFactors(factors);
        fg.validatVariables();
        
        System.out.println("Total variables: " + fg.getVariables().size());
        System.out.println("Total factors: " + fg.getFactors().size());
        for (Factor factor : fg.getFactors()) {
            List<Variable> vars = factor.getVariables();
            System.out.println(factor + ": " + vars.size());
        }
        System.out.println();
        
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setDebug(true);
        lbp.setFactorGraph(fg);
        
        Map<Variable, Float> xToValue = new HashMap<Variable, Float>();
        float[] xValues = new float[]{0.0f, 2.0f};
        xValues = new float[]{2.0f};

        // Run three times to ensure the results are the same
        for (int k = 0; k < 1; k++) {
            lbp.clearObservation();
            lbp.runInference();
            double[] aValues = a.getBelief();
            double[] bValues = b.getBelief();
            System.out.println("Belief for A:");
            for (int i = 0; i < aValues.length; i++)
                System.out.println("State " + i + ": " + aValues[i]);
            System.out.println("Belief for B:");
            for (int i = 0; i < bValues.length; i++)
                System.out.println("State " + i + ": " + bValues[i]);
            for (float xValue : xValues) {
                xToValue.put(x, xValue);
                lbp.setObservation(xToValue);
                lbp.runInference();
                System.out.println("X = " + xValue);
                System.out.println("Belief for A:");
                for (int i = 0; i < aValues.length; i++)
                    System.out.println("State " + i + ": " + aValues[i]);
                System.out.println("Belief for B:");
                for (int i = 0; i < bValues.length; i++)
                    System.out.println("State " + i + ": " + bValues[i]);
                // Just the distribution for the continuous value
                System.out.println("Continuous probabilities:");
                for (int i = 0; i < distributions.size(); i++) {
                    System.out.println(distributions.get(i).getDensity(xValue));
                }
            }
            System.out.println();
        }
    }
    
    @Test
    public void testCompetitiveReactions() throws Exception {
        FileUtility.initializeLogging();
        // The BN is something like this:
        // A -> B -> C
        //   -> D -> E (A links to both B and D)
        String[] varNames = new String[] {
                "A", "B", "C", "D", "E", "F", "G"
        };
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (String varName : varNames) {
            Variable var = new Variable(3);
            var.setName(varName);
            nameToVar.put(var.getName(), var);
        }
        double[] values = new double[] {
                0.99,
                0.005,
                0.005,
                0.005,
                0.99,
                0.005,
                0.005,
                0.005,
                0.99
        };
        FactorGraph fg = new FactorGraph();
        Factor factor = new Factor(nameToVar.get("A"),
                                   nameToVar.get("B"),
                                   values);
        fg.addFactor(factor);
        factor = new Factor(nameToVar.get("B"),
                                   nameToVar.get("C"),
                                   values);
        fg.addFactor(factor);
        factor = new Factor(nameToVar.get("A"),
                                   nameToVar.get("D"),
                                   values);
        fg.addFactor(factor);
        factor = new Factor(nameToVar.get("D"),
                                   nameToVar.get("E"),
                                   values);
        fg.addFactor(factor);
        factor = new Factor(nameToVar.get("D"),
                            nameToVar.get("F"),
                            values);
        fg.addFactor(factor);
        factor = new Factor(nameToVar.get("A"),
                            nameToVar.get("G"),
                            values);
        fg.addFactor(factor);
        
        fg.validatVariables();
        fg.setIdsInFactors();
//        fg.exportFG(System.out);
        
        // Prior inference
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setFactorGraph(fg);
        lbp.runInference();
        System.out.println("Prior:");
        outputBelief(fg);
        
        System.out.println("\nIf D = 2:");
        Observation<Number> obs = new Observation<Number>();
        obs.addAssignment(nameToVar.get("D"), 2);
        lbp.setObservation(obs);
        lbp.runInference();
        outputBelief(fg);
        
        System.out.println("\nIf F = 2:");
        obs = new Observation<Number>();
        obs.addAssignment(nameToVar.get("F"), 2);
        lbp.setObservation(obs);
        lbp.runInference();
        outputBelief(fg);
        
        System.out.println("\nIf D = 2, F = 2:");
        obs = new Observation<Number>();
        obs.addAssignment(nameToVar.get("F"), 2);
        obs.addAssignment(nameToVar.get("D"), 2);
        lbp.setObservation(obs);
        lbp.runInference();
        outputBelief(fg);
        
        System.out.println("\nIf F = 2, G = 1:");
        obs = new Observation<Number>();
        obs.addAssignment(nameToVar.get("F"), 2);
        obs.addAssignment(nameToVar.get("G"), 1);
        lbp.setObservation(obs);
        lbp.runInference();
        outputBelief(fg);
    }
    
    private Factor createFactor(Variable v1, Variable v2, double[] values) {
        Factor factor = new Factor();
        List<Variable> vars = new ArrayList<Variable>();
        vars.add(v1);
        vars.add(v2);
        factor.setVariables(vars);
        factor.setValues(values);
        return factor;
    }
    
    @Test
    public void testRunCompetitiveModel() throws Exception {
        FileUtility.initializeLogging();
        FactorGraph fg = TestUtilities.createCompetitiveFB();
        fg.exportFG(System.out);
        System.out.println("\nPrior:");
        performLBP(fg, null);
        Variable a = TestUtilities.getVariable(fg, "A");
        Variable b = TestUtilities.getVariable(fg, "B");
        Variable c = TestUtilities.getVariable(fg, "C");
        Map<Variable, Integer> observation = new HashMap<Variable, Integer>();
        observation.put(b, 0);
        System.out.println("\nB = 0");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 1);
        System.out.println("\nB = 1");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 2);
        System.out.println("\nB = 2");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 0);
        observation.put(a, 1);
        System.out.println("\nA = 1, B = 0");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 1);
        observation.put(a, 1);
        System.out.println("\nA = 1, B = 1");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 2);
        observation.put(a, 1);
        System.out.println("\nA = 1, B = 2");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 2);
        observation.put(a, 2);
        System.out.println("\nA = 2, B = 2");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 1);
        observation.put(a, 1);
        observation.put(c, 1);
        System.out.println("\nA = 1, B = 1, C = 1:");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(b, 2);
        observation.put(a, 1);
        observation.put(c, 1);
        System.out.println("\nA = 1, B = 2, C = 1:");
        performLBP(fg, observation);
        
    }
    
    @Test
    public void testRunFeedbackLoop() throws Exception {
        FileUtility.initializeLogging();
        FactorGraph fg = TestUtilities.createFeedbackLoopFG();
        fg.exportFG(System.out);
        System.out.println("\nPrior:");
        performLBP(fg, null);
        Variable a = TestUtilities.getVariable(fg, "A");
        Variable b = TestUtilities.getVariable(fg, "B");
        Variable c = TestUtilities.getVariable(fg, "C");
        Map<Variable, Integer> observation = new HashMap<Variable, Integer>();
        observation.put(a, 2);
        System.out.println("\nA = 2");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(a, 0);
        System.out.println("\nA == 0");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(a, 1);
        System.out.println("\nA == 1");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(c, 0);
        System.out.println("\nC == 0");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(c, 1);
        System.out.println("\nC == 1");
        performLBP(fg, observation);
        
        observation.clear();
        observation.put(c, 2);
        System.out.println("\nC == 2");
        performLBP(fg, observation);
    }
    
    @Test
    public void testRunInference() throws InferenceCannotConvergeException {
        FileUtility.initializeLogging();
        FactorGraph fg = TestUtilities.createSimpleFG();
        
        Variable variable = TestUtilities.getVariable(fg, "mRNA");
        Map<Variable, Integer> observation = new HashMap<Variable, Integer>();
        observation.put(variable, 2);
        performLBP(fg, observation);
//        if (true)
//            return;
        
        lbp.setUseLogSpace(true);
        lbp.runInference();
        System.out.println("\n\nIn the log space:");
        System.out.println("iteration: " + lbp.getIteration());
        System.out.println("maxDiff: " + lbp.getMaxDiff());
        outputBelief(fg);
    }

    private void performLBP(FactorGraph fg, Map<Variable, Integer> observation) throws InferenceCannotConvergeException {
        lbp.setObservation(observation);
        
        lbp.setFactorGraph(fg);
        lbp.setDebug(true);
        lbp.runInference();
        
        System.out.println("In the probability space:");
        System.out.println("iteration: " + lbp.getIteration());
        System.out.println("maxDiff: " + lbp.getMaxDiff());
        outputBelief(fg);
        System.out.println("LogZ: " + lbp.calculateLogZ());
    }

    private void outputBelief(FactorGraph fg) {
        StringBuilder builder = new StringBuilder();
        List<Variable> varList = new ArrayList<Variable>(fg.getVariables());
        Collections.sort(varList, new Comparator<Variable>() {
            public int compare(Variable var1, Variable var2) {
                return var1.getName().compareTo(var2.getName());
            }
        });
        for (Variable var : varList)  {
            double[] belief = var.getBelief();
            for (double v : belief)
                builder.append(v).append("\t");
            System.out.println(var.getName() + ": " + builder.toString());
            builder.setLength(0);
        }
    }
    
}
