/*
 * Created on Apr 18, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.HashMap;
import java.util.Map;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.booleannetwork.BooleanVariable;

/**
 * This class is used to manage converting from physical entities to BooleanVariable.
 * @author gwu
 *
 */
public class BNVariableManager {
    private Map<GKInstance, BooleanVariable> entityToVar;
    
    /**
     * Default constructor.
     */
    public BNVariableManager() {
        entityToVar = new HashMap<>();
    }
    
    public void reset() {
        entityToVar.clear();
    }
    
    public BooleanVariable getVariable(GKInstance inst) throws Exception {
        BooleanVariable var = entityToVar.get(inst);
        if (var != null)
            return var;
        var = new BooleanVariable();
        var.setName(inst.getDisplayName());
        var.addProperty("reactomeId", inst.getDBID() + "");
        
        attachGeneName(var, inst);
        
        entityToVar.put(inst, var);
        return var;
    }
    
    private void attachGeneName(BooleanVariable var,
                                GKInstance inst) throws Exception {
        if (!inst.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) 
            return;
        GKInstance refGene = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.referenceEntity);
        if (!refGene.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)) // for protein only
            return;
        String gene = (String) refGene.getAttributeValue(ReactomeJavaConstants.geneName);
        if (gene == null)
            return;
        var.addProperty("gene", gene);
    }
    
    public BooleanVariable createAccessoryNode(String name) {
        BooleanVariable var = new BooleanVariable();
        var.setName(name);
        var.setUseIdenticalTransfer(true);
        return var;
    }
    
}
