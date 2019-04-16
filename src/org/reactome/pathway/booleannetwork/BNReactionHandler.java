/*
 * Created on Apr 18, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.booleannetwork.BooleanNetwork;
import org.reactome.booleannetwork.BooleanRelation;
import org.reactome.booleannetwork.BooleanVariable;

/**
 * This class is used to handle converting a reaction into a set of BooleanRelations.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class BNReactionHandler {
    
    /**
     * Default constructor.
     */
    public BNReactionHandler() {
    }
    
    /**
     * Inputs, catalyst, activator, and inhibitor (negated) will be linked together
     * as AND. If multiple activators or inhibitors exist, a virtual node will be created
     * and activators will be linked to this virtual node using OR.
     * @param reaction
     * @param varManager
     * @throws Exception
     */
    public void handleReaction(GKInstance reaction,
                               BNVariableManager varManager,
                               BooleanNetwork network) throws Exception {
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        BooleanRelation relation = new BooleanRelation();
        network.addRelation(relation);
        for (GKInstance input : inputs) {
            BooleanVariable var = varManager.getVariable(input);
            relation.addInputVariable(var, false);
        }
        List<GKInstance> cas = reaction.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        for (GKInstance ca : cas) {
            GKInstance catalyst = (GKInstance) ca.getAttributeValue(ReactomeJavaConstants.physicalEntity);
            if (catalyst == null)
                continue;
            BooleanVariable var = varManager.getVariable(catalyst);
            relation.addInputVariable(var, false);
        }
        handleRegulations(reaction, varManager, relation, network);
        handleOutputs(reaction, varManager, relation, network);
    }
    
    private void handleOutputs(GKInstance reaction,
                               BNVariableManager varManager,
                               BooleanRelation relation,
                               BooleanNetwork network) throws Exception {
        List<GKInstance> outputs = reaction.getAttributeValuesList(ReactomeJavaConstants.output);
        if (outputs == null || outputs.size() == 0)
            return;
        if (outputs.size() == 1) {
            GKInstance output = outputs.get(0);
            BooleanVariable var = varManager.getVariable(output);
            relation.setOutputVariable(var);
        }
        else {
            BooleanVariable accessory = varManager.createAccessoryNode(reaction.getDBID() + "_output");
            relation.setOutputVariable(accessory);
            for (GKInstance output : outputs) {
                BooleanVariable outputVar = varManager.getVariable(output);
                BooleanRelation outputRelation = new BooleanRelation();
                outputRelation.addInputVariable(accessory, false);
                outputRelation.setOutputVariable(outputVar);
                network.addRelation(outputRelation);
            }
        }
    }
    
    private void handleRegulations(GKInstance reaction,
                                   BNVariableManager varManager,
                                   BooleanRelation relation,
                                   BooleanNetwork network) throws Exception {
        // Sort Regulation first
        Collection<GKInstance> regulations = InstanceUtilities.getRegulations(reaction);
        if (regulations == null || regulations.size() == 0)
            return;
        List<GKInstance> positiveRegulations = new ArrayList<>();
        List<GKInstance> negativeRegulations = new ArrayList<>();
        for (GKInstance regulation : regulations) {
            // There are only two types of regulations
            if (regulation.getSchemClass().isa(ReactomeJavaConstants.PositiveRegulation))
                positiveRegulations.add(regulation);
            else if (regulation.getSchemClass().isa(ReactomeJavaConstants.NegativeRegulation))
                negativeRegulations.add(regulation);
        }
        if (positiveRegulations.size() == 1) {
            GKInstance regulation = positiveRegulations.get(0);
            handleRegulation(regulation, varManager, false, relation);
        }
        else if (positiveRegulations.size() > 1) {
            handleRegulations(reaction, varManager, relation, positiveRegulations, false, network);
        }
        if (negativeRegulations.size() == 1) {
            GKInstance regulation = negativeRegulations.get(0);
            handleRegulation(regulation, varManager, true, relation);
        }
        else if (negativeRegulations.size() > 1) {
            handleRegulations(reaction, varManager, relation, negativeRegulations, true, network);
        }
    }

    private void handleRegulations(GKInstance reaction, 
                                   BNVariableManager varManager, 
                                   BooleanRelation relation,
                                   List<GKInstance> positiveRegulations,
                                   boolean isNegative,
                                   BooleanNetwork network) throws Exception {
        // Need to create a virtual node 
        String name = reaction.getDBID() + "";
        if (isNegative)
            name += "_inhibitor";
        else
            name += "_activator";
        BooleanVariable accessoryVar = varManager.createAccessoryNode(name);
        relation.addInputVariable(accessoryVar, isNegative);
        for (GKInstance regulation : positiveRegulations) {
            GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
            if (regulator != null) {
                BooleanVariable var = varManager.getVariable(regulator);
                // New relation is needed
                BooleanRelation relation1 = new BooleanRelation();
                relation1.addInputVariable(var, false);
                relation1.setOutputVariable(accessoryVar);
                network.addRelation(relation1);
            }
        }
    }

    private void handleRegulation(GKInstance regulation, 
                                  BNVariableManager varManager,
                                  boolean isNegative,
                                  BooleanRelation relation) throws Exception {
        GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
        if (regulator != null) {
            BooleanVariable var = varManager.getVariable(regulator);
            relation.addInputVariable(var, isNegative);
        }
    }
    
}
