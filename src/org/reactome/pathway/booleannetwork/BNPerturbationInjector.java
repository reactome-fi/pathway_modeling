/*
 * Created on Oct 2, 2018
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.booleannetwork.BooleanNetwork;
import org.reactome.booleannetwork.BooleanVariable;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to inject peturbation caused by drugs or mutations by looking for appropriate nodes
 * for injecting perturbation based on the following:
 * 1). Multiple nodes may be marked for a single gene (e.g. EWASes have been PTM). These nodes may be used
 * as inputs in multiple places, where they should be injected. 
 * 2). However, nodes refering to the same gene may be involved in a cycle. In this case, find one node with
 * the smallest DB_ID for injection.
 * @author wug
 *
 */
public class BNPerturbationInjector {
    
    public BNPerturbationInjector() {
        
    }
    
    @Test
    public void testInject() throws Exception {
        FileUtility.initializeLogging();
        Long dbId = 5693567L; // A sub-pathway in DNA repair
        dbId = 400253L; // Circadian Clock
        dbId = 5673001L; // Signaling by FGFR1
//        Map<String, Double> varToInhibition = new HashMap<>();
//        String varName = "NCOR1";
//        varToInhibition.put(varName,  0.99d);
//        varName = "HDAC3";
//        varToInhibition.put(varName, 0.99d);

        MySQLAdaptor dba = new MySQLAdaptor("localhost", "reactome_63_plus_i", "root", "macmysql01");
        GKInstance pathway = dba.fetchInstance(dbId);
        PathwayToBooleanNetworkConverter converter = new PathwayToBooleanNetworkConverter();
        BooleanNetwork network = converter.convert(pathway);
        inject(network, null, null, null, null);
    }
    
    public void inject(BooleanNetwork network,
                       Map<String, Double> geneToInhibition,
                       Map<String, Double> geneToActivation,
                       Map<BooleanVariable, Double> varToInhibition,
                       Map<BooleanVariable, Double> varToActivation) {
        if (geneToInhibition == null)
            geneToInhibition = new HashMap<>();
        if (geneToActivation == null)
            geneToActivation = new HashMap<>();
        if (varToInhibition == null)
            varToInhibition = new HashMap<>();
        if (varToActivation == null)
            varToActivation = new HashMap<>();
        Map<String, Set<BooleanVariable>> geneToVars = getGeneToVars(network);
        for (String gene : geneToVars.keySet()) {
            if (!geneToInhibition.containsKey(gene) && !geneToActivation.containsKey(gene))
                continue;
            Double inhibition = geneToInhibition.get(gene);
            Double activation = geneToActivation.get(gene);
            Set<BooleanVariable> vars = geneToVars.get(gene);
            Set<List<BooleanVariable>> sortedLists = sortVars(vars);
            for (List<BooleanVariable> list : sortedLists) {
                BooleanVariable first = list.get(0);
                if (inhibition != null && !varToInhibition.containsKey(first)) {
                    varToInhibition.put(first, inhibition);
                }
                if (activation != null && !varToActivation.containsKey(first))
                    varToActivation.put(first, activation);
            }
//            if (gene.equals("BRAF"))
//                sortedLists.forEach(list -> System.out.println(list));
        }
    }
    
    private Map<String, Set<BooleanVariable>> getGeneToVars(BooleanNetwork network) {
        Map<String, Set<BooleanVariable>> geneToVars = new HashMap<>();
        network.getVariables().forEach(var -> {
            String gene = var.getProperty("gene");
            if (gene == null)
                return;
            geneToVars.compute(gene, (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(var);
                return set;
            });
        });
        return geneToVars;
    }
    
    /**
     * Sort a set of BooleanVariables into a set of ordered BooleanVariables. The first
     * variable should be injected for perturbation.
     * @param vars
     * @return
     */
    private Set<List<BooleanVariable>> sortVars(Set<BooleanVariable> vars) {
        Set<List<BooleanVariable>> rtn = new HashSet<>();
        Set<BooleanVariable> checked = new HashSet<>();
        while (vars.size() > 0) {
            // Just pick a start
            BooleanVariable first = vars.stream().findFirst().get();
            Set<BooleanVariable> set = new HashSet<>();
            checked.clear();
            set.add(first);
            while (true) {
                int preSize = set.size();
                for (BooleanVariable var : set) {
                    if (checked.contains(var))
                        continue;
                    collectVars(var, set, vars);
                    checked.add(var);
                    break; // Start another round to avoid exception
                }
                int afterSize = set.size();
                if (preSize == afterSize)
                    break;
            }
            // We may get some accessory vars in the set. Do a filtering
            set.retainAll(vars);
            vars.removeAll(set);
            List<BooleanVariable> sortedList = getSortedList(set);
            rtn.add(sortedList);
        }
        return rtn;
    }
    
    private List<BooleanVariable> getSortedList(Set<BooleanVariable> vars) {
        List<BooleanVariable> list = new ArrayList<>(vars);
        // Find if there is a BooleanVariable that has no inputs. If true, put it at the top
        BooleanVariable found = null;
        for (BooleanVariable var : list) {
            if (var.getInRelations() == null || var.getInRelations().size() == 0) {
                found = var;
                break;
            }
        }
        if (found != null) {
            list.remove(found);
            list.add(0, found);
            return list;
        }
        // Sort based on DB_IDs
        list.sort((var1, var2) -> {
            String id1 = var1.getProperty("reactomeId");
            if (id1 == null)
                id1 = "0";
            String id2 = var2.getProperty("reactomeId");
            if (id2 == null)
                id2 = "0";
            return new Long(id1).compareTo(new Long(id2));
        });
        return list;
    }
    
    private void collectVars(BooleanVariable first,
                             Set<BooleanVariable> set,
                             Set<BooleanVariable> vars) {
        if (first.getInRelations() != null) {
            first.getInRelations().forEach(rel -> {
                rel.getInputVariables().forEach(var -> {
                    set.add(var);
                });
            });
        }
        if (first.getOutRelations() != null) {
            first.getOutRelations().forEach(rel -> {
                set.add(rel.getOutputVariable());
            });
        }
    }
    
}
