/*
 * Created on Jun 10, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlTransient;

/**
 * This class is used to describe one observation (aka a sample)
 * @author gwu
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class Observation<T extends Number> {
    private String name; // Usually sample name
    private String annoation; // Some other information except sample name (e.g. type of sample)
    // The assigned Variable values
    private List<VariableAssignment<T>> variableAssigments;
    // This map is for easy query
    @XmlTransient
    private Map<Variable, VariableAssignment<T>> varToAssign;
    
    /**
     * Default constructor.
     */
    public Observation() {
    }

    public String getAnnoation() {
        return annoation;
    }

    public void setAnnoation(String annoation) {
        this.annoation = annoation;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    public void setVariableAssignments(List<VariableAssignment<T>> variableAssignments) {
        this.variableAssigments = variableAssignments;
        varToAssign = new HashMap<Variable, VariableAssignment<T>>();
        for (VariableAssignment<T> varAssgn : variableAssignments)
            varToAssign.put(varAssgn.getVariable(), varAssgn);
    }
    
    public List<VariableAssignment<T>> getVariableAssignments() {
        return this.variableAssigments;
    }
    
    /**
     * Get the VariableAssignment object for the passed Variable object.
     * @param var
     * @return
     */
    public VariableAssignment<T> getVariableAssignment(Variable var) {
        if (variableAssigments == null)
            return null;
        return varToAssign.get(var);
    }

    public Map<Variable, T> getVariableToAssignment() {
        Map<Variable, T> varToAssign = new HashMap<Variable, T>();
        if (variableAssigments != null) {
            for (VariableAssignment<T> varAssign : variableAssigments)
                varToAssign.put(varAssign.getVariable(),
                                varAssign.getAssignment());
        }
        return varToAssign;
    }

    public void setVariableToAssignment(Map<Variable, T> variableToAssignment) {
        setUpVariableAssignments(true);
        for (Variable var : variableToAssignment.keySet()) {
            T assign = variableToAssignment.get(var);
            VariableAssignment<T> varAssgn = new VariableAssignment<T>();
            varAssgn.setVariable(var);
            varAssgn.setAssignment(assign);
            variableAssigments.add(varAssgn);
            varToAssign.put(var, varAssgn);
        }
    }

    public void setUpVariableAssignments(boolean needClear) {
        if (variableAssigments == null)
            variableAssigments = new ArrayList<VariableAssignment<T>>();
        if (varToAssign == null)
            varToAssign = new HashMap<Variable, VariableAssignment<T>>();
        if (needClear) {
            variableAssigments.clear();
            varToAssign.clear();
        }
    }
    
    public void addAssignment(Variable var, T assignment) {
        VariableAssignment<T> varAssgn = new VariableAssignment<T>();
        varAssgn.setVariable(var);
        varAssgn.setAssignment(assignment);
        addAssignment(varAssgn);
    }
    
    public void addAssignment(VariableAssignment<T> varAssgn) {
        setUpVariableAssignments(false);
        variableAssigments.add(varAssgn);
        varToAssign.put(varAssgn.getVariable(), varAssgn);
    }
    
    /**
     * Make a copy of this Observation object.
     * TODO: This copy doesn't support EmpiricalDistribution wrapped in
     * VariableAssignment object. This should be fixed soon!
     * @return
     */
    public Observation<T> copy() {
        Observation<T> copy = new Observation<T>();
        copy.name = name;
        copy.annoation = annoation;
        copy.setVariableAssignments(variableAssigments);
        return copy;
    }
    
    // Some utility methods related to Observation
    public static Observation<Integer> convertFromFloatToInteger(Observation<Number> obs) {
        Observation<Integer> rtn = new Observation<Integer>();
        rtn.annoation = obs.annoation;
        rtn.name = rtn.name;
        Map<Variable, Integer> rtnVarToState = new HashMap<Variable, Integer>();
        Map<Variable, ? extends Number> varToState = obs.getVariableToAssignment();
        for (Variable var : varToState.keySet()) {
            Number state = varToState.get(var);
            rtnVarToState.put(var, state.intValue());
        }
        rtn.setVariableToAssignment(rtnVarToState);
        return rtn;
    }
    
    public static List<Observation<Integer>> convertListFromNumberToInteger(List<Observation<Number>> list) {
        List<Observation<Integer>> rtn = new ArrayList<Observation<Integer>>();
        for (Observation<Number> obj : list) {
            Observation<Integer> obj1 = convertFromFloatToInteger(obj);
            rtn.add(obj1);
        }
        return rtn;
    }
}
