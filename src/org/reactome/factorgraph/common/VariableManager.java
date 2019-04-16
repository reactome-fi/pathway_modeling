/*
 * Created on Dec 2, 2014
 *
 */
package org.reactome.factorgraph.common;

import java.util.HashMap;
import java.util.Map;

import org.reactome.factorgraph.Variable;

/**
 * This class is used to manage variables.
 * @author gwu
 *
 */
public class VariableManager {
    private Map<String, Variable> nameToVar;
    
    /**
     * Default constructor.
     */
    public VariableManager() {
        nameToVar = new HashMap<String, Variable>();
    }
    
    public void reset() {
        nameToVar.clear();
    }
    
    /**
     * Register this variable.
     * @param variable
     */
    public void register(Variable variable) {
        nameToVar.put(variable.getName(), variable);
    }
    
    /**
     * Get the created Variable having the passed name. If such a Variable doesn't
     * exist, null will be returned.
     * @param name
     * @return
     */
    public Variable getVarForName(String name) {
        return nameToVar.get(name);
    }
    
    /**
     * Get the Variable having the passed name. If such a Variable doesn't exist,
     * it will be created and then returned.
     * @param name
     * @param variableStates
     * @return
     */
    public Variable getVarForName(String name,
                                  Integer variableStates) {
        Variable var = nameToVar.get(name);
        if (var != null)
            return var;
        var = createVar(name, variableStates);
        return var;
    }

    public Variable createVar(String name,
                              Integer variableStates) {
        Variable var = new Variable();
        if (variableStates != null)
            var.setStates(variableStates);
        var.setName(name);
        var.setId(nameToVar.size());
        nameToVar.put(name, var);
        return var;
    }
    
    public Map<String, Variable> getNameToVar() {
        return this.nameToVar;
    }
    
}
