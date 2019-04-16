/*
 * Created on Jul 20, 2015
 *
 */
package org.reactome.factorgraph.common;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.reactome.factorgraph.ContinuousVariable;
import org.reactome.factorgraph.ContinuousVariable.DistributionType;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;

/**
 * This class groups methods related to process Observation objects.
 * @author gwu
 *
 */
public class ObservationHelper {
    
    /**
     * Default constructor.
     */
    public ObservationHelper() {
    }
    
    /**
     * Create a base Observation for the specified dataType. The dataType can be null, which
     * will generate a base Observation for all touched DataTypes. This method can be used for
     * Observation containing ContinuousVariables or two-state discrete Variables only.
     * @param observations
     * @param dataType
     * @param baseState the base (aka normal state) for discrete Variables.
     * @return
     */
    public Observation<Number> createBaseObservation(List<Observation<Number>> observations,
                                                     DataType dataType,
                                                     int baseState) {
        Map<Variable, EmpiricalDistribution> varToDist = new HashMap<Variable, EmpiricalDistribution>();
        for (Observation<Number> observation : observations) {
            List<VariableAssignment<Number>> varAssgns = observation.getVariableAssignments();
            for (VariableAssignment<Number> varAssgn : varAssgns) {
                varToDist.put(varAssgn.getVariable(), varAssgn.getDistribution());
            }
        }
        Observation<Number> base = new Observation<Number>();
        String baseName = "Base";
        if (dataType != null)
            baseName += "_" + dataType;
        base.setName(baseName);
        for (Variable var : varToDist.keySet()) {
            if (dataType != null && !var.getName().endsWith("_" + dataType))
                continue; // In this case, we want to initialize variables for the specified data type only.
            VariableAssignment<Number> varAssgn = new VariableAssignment<Number>();
            varAssgn.setVariable(var);
            EmpiricalDistribution dist = varToDist.get(var);
            if (dist != null && var instanceof ContinuousVariable) {
                ContinuousVariable cVar = (ContinuousVariable) var;
                if (cVar.getDistributionType() == DistributionType.TWO_SIDED) {
                    varAssgn.setAssignment(dist.inverseCumulativeProbability(0.50d));
//                    varAssgn.setAssignment(dist.getNumericalMean());
                }
                else
                    varAssgn.setAssignment(dist.getSupportLowerBound());
            }
            else
                varAssgn.setAssignment(baseState);
            varAssgn.setDistribution(varToDist.get(var));
            base.addAssignment(varAssgn);
        }
        return base;
    }
    
    /**
     * Use this method to split an Observation into two: one is for mutation and another for non-mutation.
     * It is possible there is only one Observation is returned if the type is either mutation or non-mutation
     * (e.g. CNV, mRNA_Expression).
     * @param observation
     * @return
     */
    public List<Observation<Number>> splitObservationIntoMutationAndNonMutation(Observation<Number> observation) {
        List<Observation<Number>> rtn = new ArrayList<Observation<Number>>();
        List<VariableAssignment<Number>> varAssgns = observation.getVariableAssignments();
        List<VariableAssignment<Number>> mutVarAssgns = new ArrayList<VariableAssignment<Number>>();
        List<VariableAssignment<Number>> nonMutVarAssgns = new ArrayList<VariableAssignment<Number>>();
        for (VariableAssignment<Number> varAssgn : varAssgns) {
            Variable variable = varAssgn.getVariable();
            if (variable.getName().endsWith(DataType.Mutation.toString())) 
                mutVarAssgns.add(varAssgn);
            else
                nonMutVarAssgns.add(varAssgn);
        }
        if (mutVarAssgns.size() == 0 || nonMutVarAssgns.size() == 0)
            rtn.add(observation); // There is no need to split
        else {
            Observation<Number> mutObs = new Observation<Number>();
            mutObs.setName(observation.getName() + "_" + DataType.Mutation);
            mutObs.setVariableAssignments(mutVarAssgns);
            rtn.add(mutObs);
            Observation<Number> nonMutObs = new Observation<Number>();
            nonMutObs.setName(observation.getName() + "_Non_" + DataType.Mutation);
            nonMutObs.setVariableAssignments(nonMutVarAssgns);
            rtn.add(nonMutObs);
        }
        return rtn;
    }
    
    /**
     * Generate observations for individual data types used in the passed observation object. For example,
     * if the passed observation object covers two data types, Mutation and mRNA_Exp, two Observation objects
     * will be generated: one for Mutation and another for mRNA_Exp, which are returned in a Map keyed by the
     * DataType objects.
     * @param observation
     * @return
     */
    public Map<DataType, Observation<Number>> splitObservationBasedOnDataTypes(Observation<Number> observation) {
        Map<DataType, Observation<Number>> typeToObs = new HashMap<DataType, Observation<Number>>();
        Set<DataType> dataTypes = getDataTypesFromObservation(observation);
        List<VariableAssignment<Number>> varAssgns = observation.getVariableAssignments();
        for (DataType dataType : dataTypes) {
            Observation<Number> dataObs = new Observation<Number>();
            dataObs.setName(observation.getName() + "_" + dataType);
            typeToObs.put(dataType, dataObs);
            // Copy related information to this new Observation object
            List<VariableAssignment<Number>> dataVarAssgns = new ArrayList<VariableAssignment<Number>>();
            for (VariableAssignment<Number> varAssgn : varAssgns) {
                Variable variable = varAssgn.getVariable();
                if (variable.getName().endsWith("_" + dataType)) {
                    dataVarAssgns.add(varAssgn); // We don't clone this varAssgn object
                }
            }
            dataObs.setVariableAssignments(dataVarAssgns);
        }
        return typeToObs;
    }
    
    /**
     * Filter observations to have all types information.
     * @param observations
     */
    public void filterObservationsToHaveSharedDataTypes(List<Observation<Number>> observations) {
        Set<DataType> dataTypes = getDataTypesFromObservations(observations);
        for (Iterator<Observation<Number>> it = observations.iterator(); it.hasNext();) {
            Observation<Number> observation = it.next();
            Set<DataType> dataTypes1 = getDataTypesFromObservation(observation);
            if (dataTypes.size() > dataTypes1.size())
                it.remove();
        }
    }
    
    public Set<DataType> getDataTypesFromObservations(List<Observation<Number>> observations) {
        // Get the total data types
        Set<DataType> dataTypes = new HashSet<DataType>();
        for (Observation<Number> observation : observations) {
            Set<DataType> dataTypes1 = getDataTypesFromObservation(observation);
            dataTypes.addAll(dataTypes1);
        }
        return dataTypes;
    }
    
    /**
     * If a Variable is an observation node, use this method to get the DataType associated with
     * it.
     * @param var
     * @return
     */
    public DataType getDataTypeForVariable(Variable var) {
        DataType[] types = DataType.values();
        String varName = var.getName();
        for (DataType type : types) {
            if (varName.endsWith("_" + type))
                return type;
        }
        return null;
    }
    
    /**
     * Get a map from genes with types to their assignments for a passed data type.
     * @param observation
     * @param dataType
     * @return
     */
    public Map<String, Number> getGeneToAssignment(Observation<Number> observation,
                                                   DataType dataType) {
        Map<Variable, Number> varToAssgns = observation.getVariableToAssignment();
        Map<String, Number> geneToAssgns = new HashMap<String, Number>();
        for (Variable var : varToAssgns.keySet()) {
            Number number = varToAssgns.get(var);
            geneToAssgns.put(var.getName(), number);
        }
        return geneToAssgns;
    }

    public Set<DataType> getDataTypesFromObservation(Observation<Number> observation) {
        Set<DataType> dataTypes = new HashSet<DataType>();
        List<VariableAssignment<Number>> varAssgns = observation.getVariableAssignments();
        for (VariableAssignment<Number> varAssgn : varAssgns) {
            Variable variable = varAssgn.getVariable();
            DataType dataType = getDataTypeForVariable(variable);
            if (dataType != null)
                dataTypes.add(dataType);
        }
        return dataTypes;
    }
    
}
