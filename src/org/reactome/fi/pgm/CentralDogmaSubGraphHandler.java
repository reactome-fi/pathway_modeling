/*
 * Created on Jun 18, 2014
 *
 */
package org.reactome.fi.pgm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.PropertyConfigurator;
import org.junit.Test;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.factorgraph.InferenceCannotConvergeException;
import org.reactome.factorgraph.InferenceType;
import org.reactome.factorgraph.LoopyBeliefPropagation;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.CentralDogmaHandler;
import org.reactome.factorgraph.common.DataType;

/**
 * This class is used to treat the central dogma part, which integrates multiple data types,
 * as a template subgraph.
 * @author gwu
 *
 */
public class CentralDogmaSubGraphHandler {
    // The subgraph related to central dogma nodes and other observation nodes
    private FactorGraph fg;
    private Map<String, Variable> nameToVariable;
    // inference engine
    private LoopyBeliefPropagation lbp;
    // Used for converting the passed evidence
    private Map<Variable, Integer> varToAssgn;
    // Cached inference results to avoid repetitive calculation for same observation types
    private Map<Integer, double[]> observationToBelief;
    
    /**
     * Default constructor.
     */
    public CentralDogmaSubGraphHandler() {
        construct();
        lbp = new LoopyBeliefPropagation();
        lbp.setFactorGraph(fg);
        // Cached values
        observationToBelief = new HashMap<Integer, double[]>();
    }
    
    public void setInferenceType(InferenceType type) {
        lbp.setInferenceType(type);
    }
    
    public InferenceType getInferenceType() {
        return lbp.getInferenceType();
    }
    
    /**
     * Calculate belief of the protein variable node that can be used to integrate into
     * the FI network FG.
     * @param dataTypeToAssignment the observed states.
     * @return
     */
    public double[] calculateBelief(Map<DataType, Integer> dataTypeToAssignment) throws InferenceCannotConvergeException {
        Integer key = generateObservationKey(dataTypeToAssignment);
        double[] belief = observationToBelief.get(key);
        if (belief != null)
            return belief;
        // Create an observation
        if (varToAssgn == null)
            varToAssgn = new HashMap<Variable, Integer>();
        else
            varToAssgn.clear();
        if (dataTypeToAssignment != null) {
            Map<Variable, Integer> variableToAssn = new HashMap<Variable, Integer>();
            for (DataType type : dataTypeToAssignment.keySet()) {
                Variable var = nameToVariable.get(type.toString());
                if (var == null)
                    throw new IllegalArgumentException(type + " is not enabled!");
                Integer assign = dataTypeToAssignment.get(type);
                if (assign != null)
                    varToAssgn.put(var, assign); // Actual assignment only
            }
        }
        lbp.setObservation(varToAssgn);
        lbp.runInference();
        belief = nameToVariable.get(FIPGMConfiguration.protein).getBelief();
        // Since belief is cahed in the variable, we have to copy this array to avoid
        // values' change
        double[] cached = new double[belief.length];
        System.arraycopy(belief, 0, cached, 0, belief.length);
        observationToBelief.put(key, cached);
        return belief;
    }
    
    private Integer generateObservationKey(Map<DataType, Integer> dataTypeToAssignment) {
        int key = 0;
        if (dataTypeToAssignment == null)
            return key;
        for (DataType type : dataTypeToAssignment.keySet()) {
            Integer value = dataTypeToAssignment.get(type);
            // In order to generate the key based on the strides,
            // null is mapped to 0, 0 to 1, and 1 to 2
            if (value == null)
                value = -1;
            key += DataType.getKeyStride(type) * (value + 1); 
        }
        return key;
    }
    
    /**
     * Calculate a prior belief without any observation.
     * @return
     */
    public double[] calculateBelief() throws InferenceCannotConvergeException {
        return calculateBelief(null);
    }
    
    /**
     * Create a subgraph.
     */
    private void construct() {
        fg = new FactorGraph();
        nameToVariable = new HashMap<String, Variable>();
        Set<Factor> factors = new HashSet<Factor>();
        // Create three central dogma variables
        Variable protein = createVariable(FIPGMConfiguration.protein);
        Variable mRNA = createVariable(FIPGMConfiguration.mRNA);
        Variable dna = createVariable(FIPGMConfiguration.DNA);
        // A helper class
        CentralDogmaHandler helper = new CentralDogmaHandler();      
        helper.setConfiguration(FIPGMConfiguration.getConfig());
        // Create factor from mRNA -> protein
        helper.createCentralDogmaFactor(mRNA, protein, factors, null);
        // Create factor from DNA -> mRNA
        helper.createCentralDogmaFactor(dna, protein, factors, null);
        fg.setFactors(new HashSet<Factor>(factors));
        fg.validatVariables();
    }
    
    public void enableDataType(DataType dataType) {
        if (nameToVariable.containsKey(dataType))
            return;
        double[] values = FIPGMConfiguration.getConfig().getDataTypeValues(dataType);
        if (values == null)
            throw new IllegalArgumentException(dataType + " has not supported yet: no values have been assigned.");
        Variable variable = createVariable(dataType.toString());
        Set<Factor> factors = fg.getFactors();
        Variable parent = null;
        switch (dataType) {
            case CNV :
                parent = nameToVariable.get(FIPGMConfiguration.DNA);
                break;
            case mRNA_EXP :
                parent = nameToVariable.get(FIPGMConfiguration.mRNA);
                break;
            case Methylation :
                parent = nameToVariable.get(FIPGMConfiguration.DNA);
                break;
            case miRNA :
                parent = nameToVariable.get(FIPGMConfiguration.mRNA);
                break;
            case Mutation :
                parent = nameToVariable.get(FIPGMConfiguration.protein);
                break;
        }
        if (parent == null)
            throw new IllegalArgumentException(dataType + " has not supported yet: no central dogma node can be assigned.");
        Factor factor = new Factor(parent, variable, values);
        factors.add(factor);
        fg.validatVariables();
    }
    
    private Variable createVariable(String name) {
        Variable var = new Variable();
        var.setName(name);;
        var.setStates(2);
        nameToVariable.put(name, var);
        return var;
    }
    
    @Test
    public void testInference() throws InferenceCannotConvergeException {
        PropertyConfigurator.configure("resources/log4j.properties");
//        enableDataType(DataType.CNV);
        enableDataType(DataType.mRNA_EXP);
//        enableDataType(DataType.miRNA);
        enableDataType(DataType.Mutation);
//        enableDataType(DataType.Methylation);
        setInferenceType(InferenceType.SUM_PRODUCT);
        StringBuilder builder = new StringBuilder();
        
        double[] belief = calculateBelief();
        for (double value : belief)
            builder.append(value).append("\t");
        System.out.println("Prior belief: " + builder.toString());
        Map<DataType, Integer> assignment = new HashMap<DataType, Integer>();
        List<DataType> types = new ArrayList<DataType>();
//        types.add(DataType.CNV);
//        types.add(DataType.Methylation);
        types.add(DataType.mRNA_EXP);
        types.add(DataType.Mutation);
        List<Integer> states = new ArrayList<Integer>();
        for (int i = 0; i < types.size(); i++)
            states.add(0);
//        System.out.println("CNV\tmRNA_Exp\tMethylation\tMutation\tp(i=0)\tp(i=1)");
        System.out.println("mRNA_Exp\tMutation\tp(i=0)\tp(i=1)");
        while (!isDone(states)) {
            for (int i = 0; i < states.size(); i++) {
                assignment.put(types.get(i), states.get(i));
            }
            for (int i = 0; i < states.size(); i++) {
                Integer state = states.get(i);
                if (state < 1) {
                    states.set(i, ++state);
                    for (int j = i - 1; j >= 0; j--) {
                        states.set(j, 0); // Just reset
                    }
                    break;
                }
            }
            testInference(builder, assignment);
        }
        for (int i = 0; i < states.size(); i++) {
            assignment.put(types.get(i), states.get(i));
        }
        testInference(builder, assignment);
//        long time1 = System.currentTimeMillis();
//        Map<DataType, Integer> assignment = new HashMap<DataType, Integer>();
//        assignment.put(DataType.mRNA_EXP, 0);
//        assignment.put(DataType.CNV, 0);
//        testInference(builder, assignment);
//        assignment.put(DataType.mRNA_EXP, 0);
//        assignment.put(DataType.CNV, 1);
//        testInference(builder, assignment);
//        assignment.put(DataType.mRNA_EXP, 1);
//        assignment.put(DataType.CNV, 0);
//        testInference(builder, assignment);
//        assignment.put(DataType.mRNA_EXP, 1);
//        assignment.put(DataType.CNV, 1);
//        testInference(builder, assignment);
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time needed for the first round of run: " + (time2 - time1));
//        // In order to check caching
//        time1 = System.currentTimeMillis();
//        for (int i = 0; i < 1; i++) {
//            assignment = new HashMap<DataType, Integer>();
//            assignment.put(DataType.mRNA_EXP, 0);
//            assignment.put(DataType.CNV, 0);
//            testInference(builder, assignment);
//            assignment.put(DataType.mRNA_EXP, 0);
//            assignment.put(DataType.CNV, 1);
//            testInference(builder, assignment);
//            assignment.put(DataType.mRNA_EXP, 1);
//            assignment.put(DataType.CNV, 0);
//            testInference(builder, assignment);
//            assignment.put(DataType.mRNA_EXP, 1);
//            assignment.put(DataType.CNV, 1);
//            testInference(builder, assignment);
//        }
//        time2 = System.currentTimeMillis();
//        System.out.println("Time needed for the second round of run: " + (time2 - time1));
    }
    
    private boolean isDone(List<Integer> states) {
        for (Integer state : states) {
            if (state < 1)
                return false;
        }
        return true;
    }
    
    private void testInference(StringBuilder builder,
                               Map<DataType, Integer> assignment) throws InferenceCannotConvergeException {
        double[] belief;
        belief = calculateBelief(assignment);
//        belief = nameToVariable.get(DataType.mRNA_EXP.toString()).getBelief();
        builder.setLength(0);
        List<DataType> types = new ArrayList<DataType>(assignment.keySet());
        Collections.sort(types);
        for (DataType type : types) {
            Integer state = assignment.get(type);
            builder.append(state + "\t");
//            builder.append(type).append("=").append(state).append(",");
        }
        builder.replace(builder.length() - 1, builder.length(), "\t");
        for (double value : belief)
            builder.append(value).append("\t");
        System.out.println(builder.toString());
    }
    
    public void mergeBeliefToFactorValues(double[] values,
                                          double[] belief) {
        List<double[]> beliefs = new ArrayList<double[]>(1);
        beliefs.add(belief);
        mergeBeliefToFactorValues(values, beliefs);
    }
    
    public void mergeBeliefToFactorValues(double[] values,
                                          List<double[]> beliefs) {
        for (double[] belief : beliefs) {
            for (int i = 0; i < values.length; i++)
                values[i] *= belief[i];
        }
        // Do a normalization
        double sum = 0.0d;
        for (double value : values)
            sum += value;
        for (int i = 0; i < values.length; i++)
            values[i] /= sum;
    }
    
}
