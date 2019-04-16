/*
 * Created on Apr 11, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.DirectedGraph;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedPseudograph;
import org.reactome.r3.util.InteractionUtilities;

/**
 * A list of utility methods.
 * @author gwu
 *
 */
public class BooleanNetworkUtilities {
    
    /**
     * Default constructor.
     */
    public BooleanNetworkUtilities() {
    }
    
    /**
     * Get a list of sorted BooleanVariables contained by the passed BooleanNetwork object.
     * @param network
     * @return
     */
    public static List<BooleanVariable> getSortedVariables(BooleanNetwork network) {
        List<BooleanVariable> vars = new ArrayList<>(network.getVariables());
        Collections.sort(vars, new Comparator<BooleanVariable>() {
            public int compare(BooleanVariable var1, BooleanVariable var2) {
                return var1.getName().compareTo(var2.getName());
            }
        });
        return vars;
    }
    
    /**
     * This method is used to get BooleanVariables inside cycles.
     * @param network
     * @return
     */
    public static Set<BooleanVariable> getVariablesInCycles(BooleanNetwork network) {
        DirectedGraph<BooleanVariable, DefaultEdge> graph = createGraph(network);
//        analyzeEdges(graph);
//        if (true)
//            System.exit(0);
        // Get cycles
        // The following method has performance issue if there are many self loops
//        SzwarcfiterLauerSimpleCycles<BooleanVariable, DefaultEdge> cycleAlg = new SzwarcfiterLauerSimpleCycles<>(graph);
//        List<List<BooleanVariable>> cycles = cycleAlg.findSimpleCycles();
//        Set<BooleanVariable> cycleVars = cycles.stream().flatMap(list -> list.stream()).collect(Collectors.toSet());
//        return cycleVars;
        // Try to use this method to avoid the performance problem
        CycleDetector<BooleanVariable, DefaultEdge> cycleDetector = new CycleDetector<>(graph);
        Set<BooleanVariable> cycleVars = graph.vertexSet()
                                              .stream()
                                              .filter(v -> cycleDetector.detectCyclesContainingVertex(v))
                                              .collect(Collectors.toSet());
        return cycleVars;
    }
    
    private static void analyzeEdges(DirectedGraph<BooleanVariable, DefaultEdge> graph) {
        // Get a list of edge loops
        Map<String, DefaultEdge> keyToLoopEdge = new HashMap<>();
        Set<String> checkedKeys = new HashSet<>();
        graph.edgeSet().forEach(edge -> {
            BooleanVariable src = graph.getEdgeSource(edge);
            BooleanVariable target = graph.getEdgeTarget(edge);
            String key = InteractionUtilities.generateFIFromGene(src.id + "",
                                                                 target.id + "");
            if (checkedKeys.contains(key)) 
                keyToLoopEdge.put(key, edge);
            else
                checkedKeys.add(key);
        });
        System.out.println("Total looped edges: " + keyToLoopEdge.size());
        keyToLoopEdge.values().forEach(System.out::println);
    }
    
    private static DirectedGraph<BooleanVariable, DefaultEdge> createGraph(BooleanNetwork network) {
    		// Create a jGraphT first to use the algorithms implemented there
    		// The graph may have both loops and multiple edges, therefore, DirectedPseudograph is used here.
    		DirectedPseudograph<BooleanVariable, DefaultEdge> graph = new DirectedPseudograph<>(DefaultEdge.class);
    		// Add all variables first
    		network.getVariables().stream().forEach(var -> graph.addVertex(var));
    		// Convert relations into edges. One Relation may need to be converted into several edges
    		// if multiple inputs are created for the Relation.
    		network.getRelations().stream().forEach(relation -> {
    			Set<BooleanVariable> inputs = relation.getInputVariables();
    			BooleanVariable output = relation.getOutputVariable();
    			inputs.stream().forEach(input -> {
    			    // To avoid duplicated edges
    			    if(!graph.containsEdge(input, output)) 
    			        graph.addEdge(input, output);
    			});
    		});
    		return graph;
    }
    
    /**
     * A simple feedback loop network: A -> B -> C -> D -| C
     * @return
     */
    public static BooleanNetwork generateFeedbackLoopBN() {
        String[] names = new String[] {"A", "B", "C", "D"};
        Map<String, BooleanVariable> nameToVar = new HashMap<>();
        for (String name : names) {
            BooleanVariable var = new BooleanVariable();
            var.setName(name);
            nameToVar.put(name, var);
        }
        BooleanNetwork network = new BooleanNetwork();
        // Add A -> B
        addRelation(nameToVar.get("A"),
                    false,
                    nameToVar.get("B"),
                    network);
        BooleanRelation bToC = addRelation(nameToVar.get("B"),
                                           false,
                                           nameToVar.get("C"),
                                           network);
        addRelation(nameToVar.get("C"),
                    false,
                    nameToVar.get("D"),
                    network);
        // The inhibition should be always combined to the activation to 
        // avoid an empty activation because of inactivity of the inhibitor.
        // The inactivity of an inhibitor doesn't really mean its inhibited
        // target is in an activate state, which results from upstream stimulation.
        bToC.addInputVariable(nameToVar.get("D"), true);
//        addRelation(nameToVar.get("D"),
//                    true,
//                    nameToVar.get("B"),
//                    network);
        network.validateVariables();
        
        return network;
    }
    
    private static BooleanRelation addRelation(BooleanVariable input,
                                               boolean isNegated,
                                               BooleanVariable output,
                                               BooleanNetwork network) {
        BooleanRelation relation = new BooleanRelation();
        relation.addInputVariable(input, isNegated);
        relation.setOutputVariable(output);
        network.addRelation(relation);
        return relation;
    }
    
}
