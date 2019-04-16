/*
 * Created on Sep 10, 2018
 *
 */
package org.reactome.pathway.booleannetwork;

import java.util.HashMap;
import java.util.Map;

/**
 * Use this class to map a variety of interaction types specified in drug/target interactions
 * into a simple inhibition and activation. This mapping is manually generated.
 * @author wug
 *
 */
public class DrugTargetInteractionTypeMapper {
    private Map<String, ModificationType> intTypeToModType;
    
    public DrugTargetInteractionTypeMapper() {
    }
    
    public ModificationType getModificationType(String interactionType) {
        if (interactionType == null)
            return ModificationType.Inhibition; // Default
        if (intTypeToModType == null)
            intTypeToModType = loadTypeMap();
        ModificationType type = intTypeToModType.get(interactionType.toUpperCase());
        if (type == null)
            type = ModificationType.Inhibition; //This is always default
        return type;
    }
    
    private Map<String, ModificationType> loadTypeMap() {
        // Now this is hard coded
        String text = "ACTIVATOR,Activation\n" + 
                "AGONIST,Activation\n" + 
                "ALLOSTERIC ANTAGONIST,Inhibition\n" + 
                "ANTAGONIST,Inhibition\n" + 
                "BLOCKER,Inhibition\n" + 
                "FULL AGONIST,Activation\n" + 
                "INHIBITION,Inhibition\n" + 
                "INHIBITOR,Inhibition\n" + 
                "INVERSE AGONIST,Inhibition\n" + 
                "NEGATIVE,Inhibition\n" + 
                "NEGATIVE ALLOSTERIC MODULATOR,Inhibition\n" + 
                "NEGATIVE MODULATOR,Inhibition\n" + 
                "PARTIAL AGONIST,Activation\n" + 
                "POSITIVE ALLOSTERIC MODULATOR,Activation\n" + 
                "POSITIVE MODULATOR,Activation\n" + 
                "SUBSTRATE,Activation";
        String[] lines = text.split("\n");
        Map<String, ModificationType> map = new HashMap<>();
        for (String line : lines) {
            String[] tokens = line.split(",");
            map.put(tokens[0], ModificationType.valueOf(tokens[1]));
        }
        return map;
    }

}
