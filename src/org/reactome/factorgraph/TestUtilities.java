/*
 * Created on Oct 28, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * @author gwu
 *
 */
public class TestUtilities {
    
    /**
     * 
     */
    public TestUtilities() {
        // TODO Auto-generated constructor stub
    }
    
    public static Variable getVariable(FactorGraph fg, String varName) {
        Variable variable = null;
        for (Variable var : fg.getVariables()) {
            if (var.getName().equals(varName)) {
                variable = var;
                break;
            }
        }
        return variable;
    }
    
    /**
     * Create a factor graph to model the following reactions
     * A + B -> X and A + C -> Y. If only B is increased, we should
     * see Y decrease.
     * @return
     */
    public static FactorGraph createCompetitiveFB() {
        int id = 0;
        Variable a = createVariable("A", id++);
        Variable b = createVariable("B", id++);
        Variable c = createVariable("C", id++);
        Variable x = createVariable("X", id++);
        Variable y = createVariable("Y", id++);
        
        Factor abx = createReactionFactor(a, b, x);
        Factor acy = createReactionFactor(a, c, y);
        FactorGraph fg = new FactorGraph();
        fg.addFactor(abx);
        fg.addFactor(acy);
        fg.validatVariables();
        fg.setIdsInFactors();
        
        return fg;
    }
    
    private static Factor createReactionFactor(Variable a,
                                               Variable b,
                                               Variable x) {
        NormalDistribution normal = new NormalDistribution();
        List<Double> values = new ArrayList<Double>();
        for (int i = 0; i < a.getStates(); i++) {
            for (int j = 0; j < b.getStates(); j++) {
                double exp = Math.sqrt((i + 1) * (j + 1));
                for (int k = 0; k < x.getStates(); k++) {
                    double diff = Math.abs(k + 1 - exp);
                    // Use two-sides p-value
                    double value = (1.0 - normal.cumulativeProbability(diff)) * 2;
                    values.add(value);
                }
            }
        }
        Factor factor = new Factor();
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(x);
        variables.add(a);
        variables.add(b);
        factor.setVariables(variables);
        factor.setValues(values);
        return factor;
    }
    
    /**
     * Use this method to create the following feedback loop
     * A -> B -> C -| A
     * @return
     */
    public static FactorGraph createFeedbackLoopFG() {
        int id = 0;
        Variable a = createVariable("A", id++);
        Variable b = createVariable("B", id++);
        Variable c = createVariable("C", id++);
        
        Factor ab = createFactor(a, b, false);
        Factor bc = createFactor(b, c, false);
        Factor ca = createFactor(c, a, true);
        Set<Factor> factors = new HashSet<Factor>();
        factors.add(ab);
        factors.add(bc);
        factors.add(ca);
        
        FactorGraph graph = new FactorGraph();
        graph.setFactors(factors);
        graph.validatVariables();
        graph.setIdsInFactors();
        return graph;
    }
    
    private static Factor createFactor(Variable v1, 
                                       Variable v2,
                                       boolean isInhibit) {
        Factor factor = new Factor();
        if (isInhibit)
            factor.setName(v1.getName() + " -| " + v2.getName());
        else
            factor.setName(v1.getName() + " -> " + v2.getName());
        List<Variable> vars = new ArrayList<Variable>();
        vars.add(v2);
        vars.add(v1);
        factor.setVariables(vars);
        List<Double> values = new ArrayList<Double>();
        NormalDistribution normal = new NormalDistribution();
        for (int i = 0; i < v1.getStates(); i++) {
            double exp = i;
            if (isInhibit)
                exp = 2 - i;
            for (int j = 0; j < v2.getStates(); j++) {
                double diff = Math.abs(exp - j);
                double value = (1.0d - normal.cumulativeProbability(diff)) * 2.0d;
                values.add(value);
            }
        }
//        double[] values = new double[] {
//                0.90,
//                0.05,
//                0.05,
//                0.05,
//                0.90,
//                0.05,
//                0.05,
//                0.05,
//                0.90
//        };
        factor.setValues(values);
        return factor;
    }
    
//    private static Factor createInhibitFactor(Variable v1, 
//                                              Variable v2) {
//        Factor factor = new Factor();
//        factor.setName(v1.getName() + " -| " + v2.getName());
//        List<Variable> vars = new ArrayList<Variable>();
//        vars.add(v2);
//        vars.add(v1);
//        factor.setVariables(vars);
////        double[] values = new double[] {
////                0.05,
////                0.05,
////                0.90,
////                0.05,
////                0.90,
////                0.05,
////                0.90,
////                0.05,
////                0.05
////        };
//        factor.setValues(values);
//        return factor;
//    }
    
    private static Variable createVariable(String name,
                                           int id) {
        Variable var = new Variable();
        var.setStates(3);
        var.setName(name);
        var.setId(id);
        return var;
    }
    
    public static FactorGraph createSimpleFG() {
        Variable var1 = new Variable();
        var1.setStates(3);
        var1.setId(1);
        var1.setName("mRNA");
        Variable var2 = new Variable();
        var2.setStates(3);
        var2.setId(2);
        var2.setName("mRNA.tab");
        double[] values = new double[] {
                0.02,
                0.01,
                0.03,
                0.89,
                0.92,
                0.76,
                0.09,
                0.07,
                0.21
        };
        Factor factor1 = new Factor();
        factor1.setId(1);
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(var1);
        variables.add(var2);
        factor1.setVariables(variables);
        List<Double> valueList = new ArrayList<Double>();
        for (double value : values)
            valueList.add(value);
        factor1.setValues(valueList);
        
        Variable var3 = new Variable();
        var3.setName("protein");
        var3.setStates(3);
        var3.setId(3);
        values = new double[] {
                0.999,
                0.0005,
                0.0005,
                0.0005,
                0.999,
                0.0005,
                0.0005,
                0.0005,
                0.999
        };
        Factor factor2 = new Factor();
        factor2.setId(2);
        variables = new ArrayList<Variable>();
        variables.add(var1);
        variables.add(var3);
        factor2.setVariables(variables);
        valueList = new ArrayList<Double>();
        for (double value : values)
            valueList.add(value);
        factor2.setValues(valueList);
        
        FactorGraph fg = new FactorGraph();
        Set<Factor> factors = new HashSet<Factor>();
        factors.add(factor1);
        factors.add(factor2);
        fg.setFactors(factors);
        fg.validatVariables();
        return fg;
    }
    
}
