/*
 * Created on Jun 25, 2014
 *
 */
package org.reactome.fi.pgm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to process a small gene regulatory network that is
 * used to learn parameters.
 * @author gwu
 *
 */
public class GeneRegulatoryNetworkHandler {
    
    /**
     * Default constructor.
     */
    public GeneRegulatoryNetworkHandler() {
    }
    
    @Test
    public void checkTFTargetInteractions() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_47_plus_i",
                                            "root", 
                                            "macmysql01");
        Collection<GKInstance> insts = dba.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
        dba.loadInstanceAttributeValues(insts,
                                        new String[]{ReactomeJavaConstants.dataSource,
                                                     ReactomeJavaConstants.definition});
        List<GKInstance> selected = new ArrayList<GKInstance>();
        for (GKInstance inst : insts) {
            GKInstance dataSource = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.dataSource);
            if (!dataSource.getDisplayName().equals("ENCODE"))
                continue;
            String definition = (String) inst.getAttributeValue(ReactomeJavaConstants.definition);
            if (definition != null && definition.contains("co-expression"))
                selected.add(inst);
        }
        System.out.println("Total selected: " + selected.size());
        Set<String> factors = new HashSet<String>();
        Set<String> targets = new HashSet<String>();
        Set<String> interactions = new HashSet<String>();
        for (GKInstance inst : selected) {
            GKInstance factor = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.factor);
            GKInstance ref = (GKInstance) factor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            String factorGene = (String) ref.getAttributeValue(ReactomeJavaConstants.geneName);
            factors.add(factorGene);
            GKInstance target = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.target);
            ref = (GKInstance) target.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            String targetGene = (String) ref.getAttributeValue(ReactomeJavaConstants.geneName);
            targets.add(targetGene);
            interactions.add(factorGene + "\t" + targetGene);
        }
        System.out.println("Factors: " + factors.size());
        System.out.println("Targets: " + targets.size());
        Set<String> shared = InteractionUtilities.getShared(factors, targets);
        System.out.println("Shared: " + shared.size());
        Set<String> total = new HashSet<String>(factors);
        total.addAll(targets);
        System.out.println("Total: " + total.size());
        for (String gene : total)
            System.out.println(gene);
        System.out.println("\nInteractions: " + interactions.size());
        for (String fi : interactions)
            System.out.println(fi);
        System.out.println("\nGraph component analysis:");
        // Get the biggest connected component
        List<Set<String>> comps = new GraphAnalyzer().calculateGraphComponents(interactions);
        for (Set<String> comp : comps) {
            System.out.println(comp.size());
        }
        Set<String> genesInComp = comps.get(0);
        System.out.println("Total genes in the biggest component: " + genesInComp.size());
        // Filter interactions to the biggest component
        Set<String> fisInComp = new HashSet<String>();
        for (String fi : interactions) {
            String[] tokens = fi.split("\t");
            if (genesInComp.contains(tokens[0]) && genesInComp.contains(tokens[1]))
                fisInComp.add(fi);
        }
        System.out.println("Total FIs in the biggest component: " + fisInComp.size());
        String fileName = FIPGMConfiguration.RESULT_DIR + "GeneRegulatoryNetwork/TF_Target_Linked_Interactions.txt";
        new FileUtility().saveInteractions(fisInComp, fileName);
    }
    
}
