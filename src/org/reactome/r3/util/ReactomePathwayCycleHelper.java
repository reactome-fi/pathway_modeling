/*
 * Created on May 16, 2017
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.gk.render.RenderableChemical;
import org.jgrapht.DirectedGraph;
import org.jgrapht.alg.cycle.SzwarcfiterLauerSimpleCycles;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedPseudograph;

/**
 * This class is used to handle cycles in Reactome pathways.
 * @author gwu
 *
 */
public class ReactomePathwayCycleHelper {
    public static boolean DEBUG = false;
    
    /**
     * Default constructor.
     */
    public ReactomePathwayCycleHelper() {
    }
    
    /**
     * Get edges invovled in cycles.
     * @param edges
     * @param dba
     * @throws Exception
     */
    private List<List<HyperEdge>> findCycles(DirectedGraph<HyperEdge, DefaultEdge> graph) {
//        graph.edgeSet().forEach(edge -> {
//            HyperEdge src = graph.getEdgeSource(edge);
//            HyperEdge target = graph.getEdgeTarget(edge);
//            System.out.println(src.getReactomeId() + "\t" + 
//                               target.getReactomeId());
//        });
        // Get cycles
        if (DEBUG)
            analyzeEdges(graph);
        SzwarcfiterLauerSimpleCycles<HyperEdge, DefaultEdge> cycleAlg = new SzwarcfiterLauerSimpleCycles<>(graph);
        return cycleAlg.findSimpleCycles();
    }
    
    private void analyzeEdges(DirectedGraph<HyperEdge, DefaultEdge> graph) {
        // Get a list of edge loops
        Map<String, DefaultEdge> keyToLoopEdge = new HashMap<>();
        Set<String> checkedKeys = new HashSet<>();
        graph.edgeSet().forEach(edge -> {
            HyperEdge src = graph.getEdgeSource(edge);
            HyperEdge target = graph.getEdgeTarget(edge);
            String key = InteractionUtilities.generateFIFromGene(src.getReactomeId() + "",
                                                                 target.getReactomeId() + "");
            if (checkedKeys.contains(key)) 
                keyToLoopEdge.put(key, edge);
            else
                checkedKeys.add(key);
        });
        System.out.println("Total looped edges: " + keyToLoopEdge.size());
    }
    
    private DirectedPseudograph<HyperEdge, DefaultEdge> createGraph(List<HyperEdge> edges) {
        // Create a jGraphT first to use the algorithms implemented there
        // The graph may have both loops and multiple edges, therefore, DirectedPseudograph is used here.
    		DirectedPseudograph<HyperEdge, DefaultEdge> graph = new DirectedPseudograph<>(DefaultEdge.class);
        // Add all edges first
        for (HyperEdge edge : edges) {
            if (edge.getReactomeId() == null)
                continue;
            graph.addVertex(edge);
        }
        // It is possible there is self edge
        for (HyperEdge edge1 : edges) {
            Long id1 = edge1.getReactomeId();
            if (id1 == null)
                continue;
            for (HyperEdge edge2 : edges) {
                Long id2 = edge2.getReactomeId();
                if (id2 == null)
                    continue;
                if (isPrecedingTo(edge1, edge2)) {
                    graph.addEdge(edge1, edge2);
                }
            }
        }
        
        return graph;
    }
    
    /**
     * Check if edge1 is preceding to edge2.
     * @param edge1
     * @param edge2
     * @return
     */
    private boolean isPrecedingTo(HyperEdge edge1, HyperEdge edge2) {
        List<Node> outputs1 = edge1.getOutputNodes();
        Set<Node> inputs2 = grepInputNodes(edge2);
        // Some of outputs in edges1 are used by edge2 as inputs
        return isSharing(outputs1, inputs2);
    }
    
    /**
     * Have to check based on Reactome Ids since aliaes may be used.
     * @param c1
     * @param c2
     * @return
     */
    private boolean isSharing(Collection<Node> c1, Collection<Node> c2) {
        for (Node node1 : c1) {
            if (escape(node1))
                continue;
            Long id1 = node1.getReactomeId();
            if (id1 == null)
                continue; // Just in case
            for (Node node2 : c2) {
                if (id1.equals(node2.getReactomeId()))
                    return true;
            }
        }
        return false;
    }
    
    /**
     * Since the method here is to get inputs for expanding, simple entity is
     * escaped. Cycles involved will not be useful for expanding.
     * @param node
     * @return
     */
    private boolean escape(Node node) {
        return node instanceof RenderableChemical;
    }
    
    private void checkCycles(List<List<HyperEdge>> cycles) {
        for (int i = 0; i < cycles.size(); i++) {
            List<HyperEdge> edges = cycles.get(i);
            StringBuilder builder = new StringBuilder();
            edges.stream().forEach(edge -> builder.append(edge).append(", "));
//            System.out.println("Cycle " + i + ": " + builder.toString());
        }
    }
    
    /**
     * Get the appropriate nodes that can be regarded as inputs in cycles.
     * @param edges
     * @return
     */
    Set<Node> grepInputNodesFromCycles(List<HyperEdge> edges) {
        Set<Node> rtn = new HashSet<>();
        DirectedGraph<HyperEdge, DefaultEdge> graph = createGraph(edges);
        List<List<HyperEdge>> cycles = findCycles(graph);
//        checkCycles(cycles);
        for (int i = 0; i < cycles.size(); i++) {
            Set<Node> collected = new HashSet<>();
            List<HyperEdge> cycle = cycles.get(i);
            if (cycle.size() == 1) {
                grepInputNodesFromCycle(cycle.get(0), collected, cycle);
            }
            else { // More than one edge in the cycle
                // Try to find the first edge that is connected to outside of the cycle
                boolean hasBeenHandled = false;
                // We want to sort all HyperEdges in a cycle based on Reactome's internal ids
                // assuming curators annotate them based on their orders in the cycle. However,
                // This may not be reliable, but should result in a consistent result from different
                // runs.
                // Assume all edges should have Reactome ids, which should be reliable
                cycle.sort((edge1, edge2) -> edge1.getReactomeId().compareTo(edge2.getReactomeId()));
                for (HyperEdge cycleEdge : cycle) {
                    Set<DefaultEdge> inEdges = new HashSet<>(graph.incomingEdgesOf(cycleEdge));
                    inEdges.removeAll(cycle);
                    if (inEdges.size() > 0) {
                        // This edge has inputs from outside
                        grepInputNodesFromCycle(cycleEdge, collected, cycle);
                        hasBeenHandled = true;
                        break;
                    }
                }
                if (!hasBeenHandled) { // This is a closed cycle
                    // Randomly pick up some nodes
                    HyperEdge cycleEdge = cycle.get(0);
                    grepInputNodesFromCycle(cycleEdge, collected, cycle);
                }
            }
            rtn.addAll(collected);
//            System.out.println("Cycle " + i + " " + cycle + "; Selected nodes: " + collected);
        }
        return rtn;
    }

    protected void grepInputNodesFromCycle(HyperEdge targetEdge,
                                           Set<Node> rtn,
                                           List<HyperEdge> cycle) {
        // Self edge
        Set<Node> inputs = grepInputNodes(targetEdge);
        for (Node input : inputs) {
            List<HyperEdge> inEdges = input.getConnectedReactions();
            inEdges.removeAll(cycle);
            if (inEdges.size() == 0) { // This node is used for this cycle only
                rtn.add(input);
            }
        }
    }
    
    private Set<Node> grepInputNodes(HyperEdge edge) {
        Set<Node> inputs = new HashSet<>();
        inputs.addAll(edge.getInputNodes());
        inputs.addAll(edge.getActivatorNodes());
        inputs.addAll(edge.getHelperNodes());
        return inputs;
    }
    
}
