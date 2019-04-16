/*
 * Created on Apr 17, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.render.HyperEdge;
import org.junit.Test;
import org.reactome.booleannetwork.Attractor;
import org.reactome.booleannetwork.BooleanNetwork;
import org.reactome.booleannetwork.BooleanNetworkUtilities;
import org.reactome.booleannetwork.BooleanVariable;
import org.reactome.booleannetwork.FuzzyLogicSimulator;
import org.reactome.booleannetwork.FuzzyLogicSimulator.ANDGateMode;
import org.reactome.booleannetwork.SimulationConfiguration;
import org.reactome.r3.util.R3Constants;
import org.reactome.r3.util.ReactomeDataUtilities;

/**
 * Methods that are used to convert a Reactome pathway into a Boolean network is grouped here.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class PathwayToBooleanNetworkConverter {
    private static Logger logger = Logger.getLogger(PathwayToBooleanNetworkConverter.class);
    private BNReactionHandler reactionHandler;
    private BNVariableManager varManager;
    private BNEntityExpandHelper expandHelper;
    
    /**
     * Default constructor.
     */
    public PathwayToBooleanNetworkConverter() {
        reactionHandler = new BNReactionHandler();
        varManager = new BNVariableManager();
        expandHelper = new BNEntityExpandHelper();
    }
    
    public void setFocusedEntities(Set<String> names) {
        expandHelper.setFocusedEntities(names);
    }
    
    /**
     * Convert a Reactome pathway into a BooleanNetwork. The method is similar to the one used to
     * convert a pathway into a factor graph. The implementation here is a simple way to require
     * all annotated entities linked to reactions as required (e.g. activator somehow works as
     * a catalyst). This may be improved in the future.
     * @param pathway
     * @focusedEntities 
     * @return
     * @throws Exception
     */
    public synchronized BooleanNetwork convert(GKInstance pathway,
                                               Collection<String> focusedEntities) throws Exception {
        if (focusedEntities == null)
            setFocusedEntities(new HashSet<String>());
        else
            setFocusedEntities(new HashSet<>(focusedEntities));
        varManager.reset(); // Should not reuse pre-generated variables
        logger.info("Converting pathway " + pathway + "...");
        // Have to make sure there is a PathwayDiagram associated with the passed pathway
        Collection<GKInstance> referrers = pathway.getReferers(ReactomeJavaConstants.representedPathway);
        if (referrers == null || referrers.size() == 0)
            throw new IllegalArgumentException(pathway + " doesn't have a PathwayDiagram associated!");
        // Get all reactions
        Set<GKInstance> containedReactions = InstanceUtilities.getContainedEvents(pathway);
        List<HyperEdge> displayedEdges = ReactomeDataUtilities.getDisplayedEdges(pathway); 
        logger.info("Total displayed edges: " + displayedEdges.size());
        
        BooleanNetwork network = new BooleanNetwork();
        network.setName(pathway.getDisplayName());
        
        PersistenceAdaptor adaptor = pathway.getDbAdaptor();
        List<HyperEdge> convertedEdges = new ArrayList<>();
        // Used to mark variables as outputs
        Set<GKInstance> allOutputs = new HashSet<>();
        for (HyperEdge edge : displayedEdges) {
            if (edge.getReactomeId() == null) {
                continue;
            }
            GKInstance inst = adaptor.fetchInstance(edge.getReactomeId());
            if (inst == null) {
                // This may be possible for some disease shared pathways.
                logger.error(pathway + "'s diagram has null object for " + edge.getReactomeId());
                continue; 
            }
            if (!containedReactions.contains(inst)) { // Shared diagrams: disease reactions should be excluded.
                logger.warn(inst + " is not contained by pathway!");
                continue;
            }
            if (inst.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                // If the reaction doesn't have output, ignore it
                List<GKInstance> outputs = inst.getAttributeValuesList(ReactomeJavaConstants.output);
                if (outputs == null || outputs.size() == 0) {
                    logger.info(inst + " doesn't have output. Escape it!");
                    continue;
                }
                // Check if there is any left side entities
                Set<GKInstance> participants = InstanceUtilities.getReactionParticipants(inst);
                participants.removeAll(outputs);
                if (participants.size() == 0) {
                    logger.info(inst + " doesn't have left-side entities. Escape it!");
                    continue;
                }
                logger.info("Handling rection " + inst);
                convertedEdges.add(edge);
                reactionHandler.handleReaction(inst,
                                               varManager,
                                               network);
                allOutputs.addAll(outputs);
            }
        }
        logger.info("Expanding entities...");
        expandHelper.augumentEntities(convertedEdges,
                                      varManager, 
                                      pathway.getDbAdaptor(),
                                      network);
        network.validateVariables();
        markOutputVars(allOutputs);
        return network;
    }
    
    private void markOutputVars(Set<GKInstance> outputs) {
        outputs.forEach(output -> {
            try {
                BooleanVariable var = varManager.getVariable(output);
                if (var == null)
                    return;
                var.addProperty("output", Boolean.TRUE + "");
            }
            catch(Exception e) {
                logger.error(e.getMessage(), e);
            }
        });
    }
    
    public BooleanNetwork convert(GKInstance pathway) throws Exception {
        return convert(pathway, new HashSet<String>()); 
    }
    
    @Test
    public void testAllPathways() throws Exception {
        String fileName = R3Constants.RESULT_DIR + "ReactomePathways090217.txt";
        Set<Long> dbIds = Files.lines(Paths.get(fileName))
                               .map(line -> new Long(line.split("\t")[0]))
                               .collect(Collectors.toSet());
        logger.info("Total pathways: " + dbIds.size());
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                "reactome_59_plus_i",
                "root",
                "macmysql01");
        dbIds.forEach(dbId -> {
            try {
                GKInstance pathway = dba.fetchInstance(dbId);
                testConvert(new HashSet<>(), pathway);
            }
            catch(Exception e) {
                logger.error(e.getMessage(), e);
            }
        });
    }
    
    @Test
    public void testConvert() throws Exception {
//        FileUtility.initializeLogging();
        
        Long dbId = 2032785L; // YAP1- and WWTR1-stimulated gene expression
        dbId = 1257604L; // PIP3 activates AKT signaling
        
        Set<String> focusedEntities = new HashSet<>();
//        focusedEntities.add("EGFR");
//        setFocusedEntities(focusedEntities);
        
//        dbId = 5621481L; // C-type lectin receptors (CLRs), having some loops
        dbId = 5693567L; // HDR through Homologous Recombination (HR) or Single Strand Annealing (SSA)
                            // A pathway has some loops
        
        dbId = 400253L; // Circadian Clock
        dbId = 71387L; // [Pathway:71387] Metabolism of carbohydrates: This pathway takes a long time 
        dbId = 73923L; // Wrong with convert to a graph. [Pathway:73923] Lipid digestion, mobilization, and transport...
        
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_59_plus_i",
                                            "root",
                                            "macmysql01");
        GKInstance pathway = dba.fetchInstance(dbId);
        BooleanNetwork network = testConvert(focusedEntities,
                                             pathway);
        
//        JAXBContext context = JAXBContext.newInstance(BooleanNetwork.class);
//        Marshaller marshaller = context.createMarshaller();
//        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
//        StringWriter writer = new StringWriter();
//        marshaller.marshal(network, writer);
//        System.out.println(writer.toString());
//        
//        if (true)
//            return;
//        
//        System.out.println("\nReading...");
//        Unmarshaller unmarshaller = context.createUnmarshaller();
//        network = (BooleanNetwork) unmarshaller.unmarshal(new StringReader(writer.toString()));
//        System.out.println("Name: " + network.getName());
//        System.out.println("Total variables: " + network.getVariables().size());
//        for (BooleanVariable var1 : network.getVariables()) {
//            if (var1.getName().equals("Insulin:p-6Y-Insulin receptor [plasma membrane]")) {
//                System.out.println(var1.getName() + ": ");
//                System.out.println("inRelations: " + var1.getInRelations().size());
//                break;
//            }
//        }
//        System.out.println("Total relatioins: " + network.getRelations().size());
//        
        if (true)
            return;
        
        FuzzyLogicSimulator simulator = new FuzzyLogicSimulator();
        // Need to use product based on biochemistry mass kinetics
        simulator.setAndGateMode(ANDGateMode.PROD);
        
        BooleanVariable checkedVar = null;
        // Assign start values 0.5 for all variables
        for (BooleanVariable var : network.getVariables()) {
            if (var.getInRelations() == null || var.getInRelations().size() == 0)
                var.setValue(1.0d);
            else
                var.setValue(0.0d);
            if (var.getName().equals("p-6Y-EGFR [plasma membrane]"))
                checkedVar = var;
        }
        Map<String, Number> varNameToStimulation = new HashMap<>();
        
//        varNameToStimulation.put("RUNX2 [nucleoplasm]", 0.5);
//        varNameToStimulation.put("WWTR1 [nucleoplasm]", 1.0);
//        varNameToStimulation.put("PI(3,4,5)P3 [plasma membrane]", 1.0);
//        varNameToStimulation.put("Activator:PI3K [plasma membrane]", 0.5);
//        varNameToStimulation.put("p-6Y-EGFR [plasma membrane]", 0.1);
//        varNameToStimulation.put("ATP [cytosol]", 1.0d);
//        varNameToStimulation.put("PI(4,5)P2 [plasma membrane]", 0.5d);
        Map<BooleanVariable, Double> varToInhibition = new HashMap<>();
        varToInhibition.put(checkedVar, 0.99);
        SimulationConfiguration configuration = new SimulationConfiguration();
        configuration.setInhibition(varToInhibition);
        
        Attractor base = null;
        simulator.simulate(network, configuration);
        
        if (simulator.isAttractorReached()) {
            System.out.println("Attractor is reached: ");
            Attractor attractor = simulator.getAttractor();
//            System.out.println(attractor.outputAsText());
            base = attractor;
        }
        else
            System.out.println("Attractor cannot be reached!");
        
        for (BooleanVariable var : network.getVariables()) {
            if (var.getName().equals("Activator:PI3K [plasma membrane]"))
                System.out.println(Arrays.asList(var.getTrack()));
        }
        
        if (true)
            return;
        
        Attractor perturbed = null;
//        varNameToStimulation.put("p-6Y-EGFR [plasma membrane]", 1.0);
//        varNameToStimulation.put("Activator:PI3K [plasma membrane]", 1.0d);
//        varNameToStimulation.put("THEM4/TRIB3 [plasma membrane]", 0.1);
//        varNameToStimulation.put("ATP [cytosol]", 1.0d);
//        varNameToStimulation.put("PI(4,5)P2 [plasma membrane]", 0.5d);
        simulator.simulate(network, varNameToStimulation);
        if (simulator.isAttractorReached()) {
            System.out.println("Attractor is reached: ");
            Attractor attractor = simulator.getAttractor();
//            System.out.println(attractor.outputAsText());
            perturbed = attractor;
        }
        else
            System.out.println("Attractor cannot be reached!");
        
        System.out.println("\nCheck Perturbation:");
        calculatePerturbation(perturbed, base);
    }

    private BooleanNetwork testConvert(Set<String> focusedEntities, GKInstance pathway) throws Exception {
        BooleanNetwork network = convert(pathway, focusedEntities);
        System.out.println("\nBooleanNetwork name: " + network.getName());
        System.out.println("Total variables: " + network.getVariables().size());
//        network.getVariables().stream().forEach(var -> {
//            System.out.println(var.getId() + "\t" + var.getName() + "\t" + var.getProperty("gene"));
//        });
        
        System.out.println("Total relations: " + network.getRelations().size());
//        network.getRelations().stream().forEach(System.out::println);
        
        Set<BooleanVariable> cycleVars = BooleanNetworkUtilities.getVariablesInCycles(network);
        System.out.println("\nThe following variables are in cycles: " + cycleVars.size());
//        cycleVars.forEach(System.out::println);
        return network;
    }
    
    private void calculatePerturbation(Attractor perturbed,
                                       Attractor base) {
        Map<BooleanVariable, List<Number>> perturnedVarToValues = perturbed.getVarToValues();
        Map<BooleanVariable, List<Number>> baseVarToValues = base.getVarToValues();
        System.out.println("Variable\tPerturbed_Value\tBase_Value\tRatio");
        for (BooleanVariable var : baseVarToValues.keySet()) {
            if (var.getName().endsWith("_output"))
                continue;
            Number perturbedValue = perturnedVarToValues.get(var).get(0);
            Number baseValue = baseVarToValues.get(var).get(0);
            double ratio = perturbedValue.doubleValue() / baseValue.doubleValue();
            if (Math.abs(ratio - 1.0) < 1.0E-6)
                continue; // Don't show anything that is the same
            System.out.println(var.getName() + "\t" + 
                               perturbedValue.doubleValue() + "\t" + 
                               baseValue.doubleValue() + "\t" + 
                               ratio);
        }
    }
    
}
