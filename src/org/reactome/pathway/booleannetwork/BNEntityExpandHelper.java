/*
 * Created on Apr 18, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.reactome.booleannetwork.BooleanNetwork;
import org.reactome.booleannetwork.BooleanRelation;
import org.reactome.booleannetwork.BooleanVariable;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.ReactomeDataUtilities;

/**
 * This class is used to expand entities in a pathway that have not been annotated explicitly,
 * e.g., Complexes and EntitySets.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class BNEntityExpandHelper {
    private static final Logger logger = Logger.getLogger(BNEntityExpandHelper.class);
    // Choose only members having entities having the following names
    // A pathway using a set can be expanded into multiple pathways materialized by each member
    // of a set. Using this control to select a sub-set of expanded pathways for mimicking 
    // tissue specific pathway, assuming that tissue uses the following set of entities
    private Set<String> focusedEntities;
    
    /**
     * Default constructor.
     */
    public BNEntityExpandHelper() {
    }
    
    public Set<String> getFocusedEntities() {
        return focusedEntities;
    }

    public void setFocusedEntities(Set<String> focusedEntities) {
        this.focusedEntities = focusedEntities;
    }

    public void augumentEntities(List<HyperEdge> edges,
                                 BNVariableManager varManager,
                                 PersistenceAdaptor dba,
                                 BooleanNetwork network) throws Exception {
        Set<GKInstance> outputInstances = new HashSet<GKInstance>();
        Set<Node> inputs = ReactomeDataUtilities.getNodesForAuguemnt(edges, dba, outputInstances);
        logger.info("Total inputs (inputs, catalysts, and regulators): " + inputs.size());
        
        // Record all instances that have been processed to avoid duplication. For example,
        // an input EntitySet may contain an input Complex.
        Set<GKInstance> processed = new HashSet<GKInstance>();
        for (Node node : inputs) {
            if (node.getReactomeId() == null)
                continue; // This should not be possible. But just in case.
            GKInstance input = dba.fetchInstance(node.getReactomeId());
            if (input == null) {
                logger.error(node.getDisplayName() + " with DB_ID = " + node.getReactomeId() + " is not in the database!");
                continue;
            }
            augumentEntity(input,
                           processed,
                           outputInstances,
                           varManager,
                           network);
        }
    }
    
    private void augumentEntity(GKInstance entity,
                                Set<GKInstance> processed,
                                Set<GKInstance> outputInstances,
                                BNVariableManager varManager,
                                BooleanNetwork network) throws Exception {
        if (processed.contains(entity) || // This has been processed
            outputInstances.contains(entity))  // It is used as an output, and should not be processed
            return;
        processed.add(entity);
        if (entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            augmentEntitySet(entity, 
                             processed,
                             outputInstances,
                             varManager,
                             network);
        else if (entity.getSchemClass().isa(ReactomeJavaConstants.Complex))
            augmentComplex(entity, 
                           processed,
                           outputInstances,
                           varManager,
                           network);
        else { // Mark the variable for later use
            BooleanVariable var = varManager.getVariable(entity);
            var.addProperty("role", "CycleInput");
        }
    }
    
    private void augmentComplex(GKInstance complex,
                                Set<GKInstance> processed,
                                Set<GKInstance> outputInstances,
                                BNVariableManager variableManager,
                                BooleanNetwork network) throws Exception {
        logger.info("Augment Complex: " + complex);
        BooleanVariable complexVar = variableManager.getVariable(complex);
        List<GKInstance> list = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        if (list == null || list.size() == 0)
            return;
        BooleanRelation relation = new BooleanRelation();
        network.addRelation(relation);
        relation.setOutputVariable(complexVar);
        for (GKInstance inst : list) {
            BooleanVariable var = variableManager.getVariable(inst);
            relation.addInputVariable(var, false);
        }
        // Need to call recursive
        for (GKInstance inst : list) {
            augumentEntity(inst,
                           processed,
                           outputInstances,
                           variableManager,
                           network);
        }
    }
    
    private void augmentEntitySet(GKInstance set,
                                  Set<GKInstance> processed,
                                  Set<GKInstance> outputInstances,
                                  BNVariableManager variableManager,
                                  BooleanNetwork network) throws Exception {
        logger.info("Augment EntitySet: " + set);
        BooleanVariable setVar = variableManager.getVariable(set);
        List<GKInstance> list = set.getAttributeValuesList(ReactomeJavaConstants.hasMember);
        Set<GKInstance> members = new HashSet<GKInstance>();
        if (list != null && list.size() > 0) {
            members.addAll(list);
        }
        if (set.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasCandidate)) {
            list = set.getAttributeValuesList(ReactomeJavaConstants.hasCandidate);
            if (list != null && list.size() > 0) {
                members.addAll(list);
            }
        }
        int preSize = members.size();
        filterMembers(members);
        if (members.size() == 0)
            throw new IllegalStateException("Cannot find any member after filtering: " + set);
        if (preSize > members.size())
            logger.info("Some members have been filtered out: " + set);
        for (GKInstance member : members) {
            BooleanVariable memberVar = variableManager.getVariable(member);
            BooleanRelation relation = new BooleanRelation();
            network.addRelation(relation);
            relation.addInputVariable(memberVar, false);
            relation.setOutputVariable(setVar);
        }
        // Need to call recursive
        for (GKInstance member : members) {
            augumentEntity(member,
                           processed,
                           outputInstances,
                           variableManager,
                           network);
        }
    }
    
    private void filterMembers(Set<GKInstance> members) throws Exception {
        if (focusedEntities == null || focusedEntities.size() == 0)
            return; // Do nothing
        Map<GKInstance, Set<String>> memberToNames = getMemberToNames(members);
        // If there are no names mapped to focusedEntities, use all members
        boolean useAll = true;
        for (GKInstance member : memberToNames.keySet()) {
            Set<String> names = memberToNames.get(member);
            Set<String> shared = InteractionUtilities.getShared(names, focusedEntities);
            if (shared.size() > 0) {
                useAll = false;
                break;
            }
        }
        if (useAll)
            return; // Nothing matched
        // Want to keep members having matched instances only
        for (Iterator<GKInstance> it = members.iterator(); it.hasNext();) {
            GKInstance member = it.next();
            Set<String> names = memberToNames.get(member);
            Set<String> shared = InteractionUtilities.getShared(names, focusedEntities);
            if (shared.size() == 0)
                it.remove();
        }
    }
    
    private Map<GKInstance, Set<String>> getMemberToNames(Set<GKInstance> members) throws Exception {
        Map<GKInstance, Set<String>> memberToNames = new HashMap<>();
        for (GKInstance member : members) {
            Set<GKInstance> refEntities = InstanceUtilities.grepReferenceEntitiesForPE(member);
            Set<String> names = new HashSet<>();
            for (GKInstance refEntity : refEntities) {
                String name = null;
                if (refEntity.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName))
                    name = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
                else
                    name = (String) refEntity.getAttributeValue(ReactomeJavaConstants.name);
                if (name != null)
                    names.add(name);
            }
            memberToNames.put(member, names);
        }
        return memberToNames;
    }
    
}
