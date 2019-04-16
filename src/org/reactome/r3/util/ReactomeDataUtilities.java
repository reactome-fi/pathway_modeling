/*
 * Created on Apr 18, 2017
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.gk.render.Renderable;
import org.gk.render.RenderableInteraction;
import org.gk.render.RenderablePathway;
import org.gk.render.RenderableReaction;

/**
 * Group some utilities related to Reactome data handling.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactomeDataUtilities {
    
    public ReactomeDataUtilities() {
    }
    
    /**
     * We may have to consider alias nodes.
     */
    private static void removeOutputsInInputs(Set<Node> inputs,
                                              Set<Node> outputs) {
        inputs.removeAll(outputs);
        // Remove based on Reactome ids for aliases
        Set<Long> outputIds = outputs.stream()
                .filter(node -> node.getReactomeId() != null)
                .map(node -> node.getReactomeId())
                .collect(Collectors.toSet());
        for (Iterator<Node> it = inputs.iterator(); it.hasNext();) {
            Node node = it.next();
            if (outputIds.contains(node.getReactomeId()))
                it.remove();
        }
    }
    
    /**
     * Get a set of displayed nodes for augumenting during mathematical modeling.
     * @param edges
     * @return
     * TODO: There is a problem with this simple implementation: If an Entity is involved
     * in a cycle, it will not be picked up. A better way is needed to find these entities.
     */
    public static Set<Node> getNodesForAuguemnt(List<HyperEdge> edges,
                                                PersistenceAdaptor dba,
                                                Set<GKInstance> outputInstances) throws Exception {
        // Search for inputs (inputs, catalysts and regulators) that are used as inputs only in this pathway so that we
        // can add something more to these entities
        Set<Node> inputs = new HashSet<Node>();
        Set<Node> outputs = new HashSet<Node>();
        for (HyperEdge edge : edges) {
            inputs.addAll(edge.getConnectedNodes());
            // Special cases: gene regulatory interactions are represented as BlackboxEvent
            // without input. Cases like those, we want to augment outputs as though they
            // are inputs.
//            if (isGeneRegulatoryInteraction(edge, dba)) 
//                continue;
            // As of Nov 11, 2014, there is no need to check the above. A gene regulatory reaction will be handled
            // during the stage of handling reactions.
            outputs.addAll(edge.getOutputNodes());
        }
        // The following line is not correct since aliases may be used
        //inputs.removeAll(outputs);
        removeOutputsInInputs(inputs, outputs);
        
        // Need to pick up some nodes from cycles. Otherwise, some of Entities involved in cycles cannot be 
        // augumented (e.g. UBE2N:UBE2V1 in http://www.reactome.org/PathwayBrowser/#/R-HSA-5607756).
        ReactomePathwayCycleHelper cycleHelper = new ReactomePathwayCycleHelper();
        Set<Node> cycleInputs = cycleHelper.grepInputNodesFromCycles(edges);
        // Now need to remove them from outputs 
        removeOutputsInInputs(outputs, cycleInputs); // Call this method in a reverse way
        
        for (Node cycleInput : cycleInputs) {
            inputs.add(cycleInput);
        }
        
        if (outputInstances != null) {
            // Record all instances that have been processed to avoid duplication. For example,
            // An input EntitySet may contain an input Complex.
            for (Node output : outputs) {
                if (output.getReactomeId() == null)
                    continue;
                outputInstances.add(dba.fetchInstance(output.getReactomeId()));
            }
        }
        return inputs;
    }
    
    /**
     * Get edges displayed in a PathwayDiagram.
     * TODO: This method is copied from PathwayToFactorGraphConverter and should be refactored
     * to a common place in the future.
     * @param pathway
     * @return
     * @throws Exception
     */
    public static List<HyperEdge> getDisplayedEdges(GKInstance pathway) throws Exception {
        Collection<GKInstance> referrers = pathway.getReferers(ReactomeJavaConstants.representedPathway);
        if (referrers == null || referrers.size() == 0)
            return null; 
        GKInstance diagramInst = referrers.iterator().next();
        RenderablePathway diagram = new DiagramGKBReader().openDiagram(diagramInst);
        // Check the displayed object
        List<HyperEdge> edges = new ArrayList<HyperEdge>();
        for (Object o : diagram.getComponents()) {
            Renderable r = (Renderable) o;
            // In this implementation, RenderableInteraction will not be handled since
            // they are actually not annotated in Reactome pathways
            if (r.getReactomeId() == null)
                continue;
            GKInstance inst = pathway.getDbAdaptor().fetchInstance(r.getReactomeId());
            if (inst == null)
                continue; // Nothing to be processed
            if (r instanceof RenderableReaction ||
                r instanceof RenderableInteraction) {
                // Want to handle these two types HyperEdges only
                edges.add((HyperEdge)r);
            }
        }
        return edges;
    }
    
}
