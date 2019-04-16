/*
 * Created on Oct 15, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.VariableManager;

/**
 * This class is used to generate and manage Variable objects during pathway converting.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class PathwayVariableManager extends VariableManager {
    // A list of entity names that should be escaped during converting.
    // For example, ATP, ADP, H2O, H+, and other small molecules should
    // not be converted into the factor graph since they should be considered
    // as constant.
    private List<String> namesForEscape;
    // Cache map from GKInstance to Variable
    private Map<GKInstance, Variable> instToVar;
    // Required by other classes whose APIs are used in this package
    // The client should not use this variable to get generated variables.
    private Map<String, Variable> nameToVarInDogma;
    
    /**
     * Default constructor.
     */
    public PathwayVariableManager() {
        namesForEscape = new ArrayList<String>();
    }
    
    public void setNamesForEscape(List<String> names) {
        this.namesForEscape = names;
        if (namesForEscape == null)
            namesForEscape = new ArrayList<String>();
    }
    
    public Map<GKInstance, Variable> getInstToVarMap() {
        return this.instToVar;
    }
    
    /**
     * Clear up cache.
     */
    public void reset() {
        super.reset();
        if (instToVar == null)
            instToVar = new HashMap<GKInstance, Variable>();
        else
            instToVar.clear();
        if (nameToVarInDogma == null)
            nameToVarInDogma = new HashMap<String, Variable>();
        else
            nameToVarInDogma.clear();
    }
    
    /**
     * Get a map of from names to variables that are generated during central dogma expanding.
     * @return
     */
    public Map<String, Variable> getNameToVariableInDogma() {
        return this.nameToVarInDogma;
    }
    
    /**
     * Check if the passed GKInstance object has been converted into a Variable object.
     * @param inst
     * @return
     */
    public boolean isInstanceConverted(GKInstance inst) {
        return instToVar.containsKey(inst);
    }
    
    /**
     * Get a Variable for the passed GKInstance object.
     * @param inst
     * @return
     * @throws Exception
     */
    public Variable getVariable(GKInstance inst) throws Exception {
        if (shouldEscape(inst)) {
            return null;
        }
        Variable var = instToVar.get(inst);
        if (var != null)
            return var;
        // Attach DB_ID after name to avoid duplicated names
        // used in reactome
        var = getVarForName(inst.getDisplayName() + "_" + inst.getDBID(), 
                            PathwayPGMConfiguration.getConfig().getNumberOfStates());
        var.setProperty(ReactomeJavaConstants.DB_ID, 
                        inst.getDBID() + "");
        instToVar.put(inst, var);
        return var;
    }

    private boolean shouldEscape(GKInstance inst) throws Exception {
        if (inst.getSchemClass().isValidAttribute(ReactomeJavaConstants.name)) {
            List<String> names = inst.getAttributeValuesList(ReactomeJavaConstants.name);
            if (names != null && names.size() > 0) {
                for (String name : names) {
                    if (namesForEscape.contains(name))
                        return true;
                }
            }
        }
        return false;
    }
    
    /**
     * A helpe method to generate a factor.
     * @param varList
     * @param typeList
     * @param valueAssigner
     * @param name
     * @return
     */
    public Factor createFactor(List<Variable> varList,
                               List<FactorEdgeType> typeList,
                               FactorValueAssigner valueAssigner,
                               String name) {
        List<Double> values = valueAssigner.generateFactorValues(varList, typeList);
        Factor factor = new Factor();
        // To avoid a danger that varList may be changed by other unaccidently, clone
        // this list.
        factor.setVariables(new ArrayList<Variable>(varList));
        factor.setValues(values);
        factor.setName(name);
        return factor;
    }
}
