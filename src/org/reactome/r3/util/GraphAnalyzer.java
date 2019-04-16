/*
 * Created on Mar 21, 2012
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

/**
 * This class is used to do some graph related file generation, e.g. generating
 * a largest graph component.
 * @author gwu
 *
 */
public class GraphAnalyzer {
    
    public GraphAnalyzer() {
    }
    
    /**
     * Calculate linked graph components from a set of FIs.
     * @param interactions
     * @return
     */
    public List<Set<String>> calculateGraphComponents(Set<String> interactions) {
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        SimpleGraph<String, DefaultEdge> graph = (SimpleGraph<String, DefaultEdge>) createGraph(ids, interactions);
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<String, DefaultEdge>(graph);
        List<Set<String>> components = inspector.connectedSets();
        List<Set<String>> componentList = new ArrayList<Set<String>>(components);
        Collections.sort(componentList, new Comparator<Set<String>>() {
            public int compare(Set<String> set1, Set<String> set2) {
                return set2.size() - set1.size();
            }
        });
        return componentList;
    }
    
    private Graph<String, DefaultEdge> createGraph(Set<String> reactomeIds, 
                                                   Set<String> interactions) {
        int index = 0;
        // Just want to use String as edge class
        Graph<String, DefaultEdge> graph = new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
        for (String id : reactomeIds)
            graph.addVertex(id);
        // Need to find the delimit index
        String tmp = interactions.iterator().next();
        String delimit = " ";
        index = tmp.indexOf("\t");
        if (index > 0)
            delimit = "\t"; 
        for (String pair : interactions) {
            index = pair.indexOf(delimit);
            //System.out.println(pair);
            graph.addEdge(pair.substring(0, index), pair.substring(index + 1));
        }
//        System.out.printf("Graph: vertics %d edges %d%n", 
//                          graph.vertexSet().size(), 
//                          graph.edgeSet().size());        
        return graph;
    }
}
