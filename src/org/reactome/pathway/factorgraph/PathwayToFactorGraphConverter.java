/*
 * Created on Oct 13, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collection;
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
import org.gk.persistence.MySQLAdaptor;
import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.junit.Test;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.factorgraph.LoopyBeliefPropagation;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.SharedEMFactors;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.DiscreteObservationFactorhandler;
import org.reactome.factorgraph.common.ObservationFactorHandler;
import org.reactome.factorgraph.common.ObservationFileLoader;
import org.reactome.factorgraph.common.ObservationFileLoader.ObservationData;
import org.reactome.factorgraph.common.ObservationRandomizer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.ReactomeDataUtilities;

/**
 * This class is used to convert a Reactome pathway into a FactorGraph. The converting happens as the following
 * steps:
 * 1). Grep all reactions from the passed pathway.
 * 2). Convert individual reactions into factors.
 * 3). Link output_nodes converted from reactions into output variables. An output may be generated by multiple
 * reactions. It should be linked to these output_nodes. If more than 3 output_nodes are linked to an output,
 * add accessory node using OR (edge label max). However, there is a special case to this rule: gene regulatory 
 * reaction generates protein A, which will be expanded by adding central dogma nodes. So protein A is going to 
 * have at least two parent nodes: output from gene regulatory reaction and output from dogma node. In this case, 
 * we treat them separately and hope to get an aggregated inference value.
 * 4). To convert a reaction into a factor, use the following procedures:
 * a). If more than one input, one activator, or one inhibitor, add accessory node, respectively.
 * b). If there are more than 3 inputs, 3 activators, or 3 inhibitors, add accessory node. AND (edge label min) should
 * be used for input, OR (edge label max) for regulators.
 * c). Label edges for catalyst, input_node and activtor_node as min. Aggregate labels for them and use the summarized
 * value for vote.
 * d). There is always an output_node state regardless how many output a reaction can generate. 
 * e). The expected state of output_node is decided by vote of inhibitor_node and vote from c).
 * f). If there is only one input, activator, inhibitor, don't add accessory node.
 * 5). Find inputs to the whole pathways, and expand EntitySet, Complex. Add central dogma nodes as PARADIGM.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class PathwayToFactorGraphConverter {
    private static final Logger logger = Logger.getLogger(PathwayToFactorGraphConverter.class);
    // We only need one copy of FactorValueAssigner
    protected FactorValueAssigner valueAssigner;
    // Manage variable generation
    private PathwayVariableManager variableManager;
    // For handling a group of variables
    protected VariableSetHandler setHandler;
    // Cache in order to check the loaded data
    private ObservationFileLoader observationLoader;
    // Helper object for entity expanding by adding central dogma nodes
    protected EntityExpandHelper expandHelper;
    // For handling reactions
    protected ReactionHandler reactionHandler;
    // A flag to avoid loading observation data files so that we can 
    // use a pathway factor graph model for perturbation analysis (e.g. drug, mutation). 
    // The default is false for PARADIGM
    private boolean isForPerturbation;
    
    /**
     * Default constructor.
     */
    public PathwayToFactorGraphConverter() {
        valueAssigner = getFactorValueAssigner();
        variableManager = new PathwayVariableManager();
        setHandler = new VariableSetHandler();
        setHandler.setValueAssigner(valueAssigner);
        setHandler.setVariableManager(variableManager);
        observationLoader = new ObservationFileLoader();
        observationLoader.setPGMConfiguration(PathwayPGMConfiguration.getConfig());
        // Set up ObservationFactorHandler
        ObservationFactorHandler obsFactorHandler = new DiscreteObservationFactorhandler();
        observationLoader.setObservationFactorHandler(DataType.mRNA_EXP, obsFactorHandler);
        observationLoader.setObservationFactorHandler(DataType.CNV, obsFactorHandler);
        // This is just a test
//        observationLoader.setObservationFactorHandler(DataType.mRNA_EXP, new EmpiricalFactorHandler());
//        observationLoader.setObservationFactorHandler(DataType.CNV, new EmpiricalFactorHandler());
        configureDogmaFactorValues();
        // Configure expandHelper
        expandHelper = new EntityExpandHelper();
        expandHelper.setForPerturbation(isForPerturbation);
        expandHelper.setVariableManager(variableManager);
        expandHelper.setVariableSetHandler(setHandler);
        expandHelper.setPGMConfiguration(PathwayPGMConfiguration.getConfig());
        // For handle reactions
        reactionHandler = new ReactionHandler();
        reactionHandler.setSetHandler(setHandler);
        reactionHandler.setValueAssigner(valueAssigner);
        reactionHandler.setVariableManager(variableManager);
        reactionHandler.setExpandHelper(expandHelper);
    }
    
    private FactorValueAssigner getFactorValueAssigner() {
        if (PathwayFGConstants.NUMBER_OF_STATES == 2)
            return new PerturbationFactorValueAssigner();
        return new FactorValueAssigner(); // Default
    }

    public boolean isForPerturbation() {
        return isForPerturbation;
    }

    public void setForPerturbation(boolean isForPerturbation) {
        this.isForPerturbation = isForPerturbation;
        this.expandHelper.setForPerturbation(isForPerturbation);
    }

    public ObservationFileLoader getObservationLoader() {
        return this.observationLoader;
    }
    
    public Map<DataType, SharedEMFactors> getTypeToSharedFactors() {
        if (observationLoader != null)
            return observationLoader.getTypeToSharedFactors();
        return null;
    }
    
    private void configureDogmaFactorValues() {
        // Just two dummy nodes for generating values
        Variable output = variableManager.getVarForName("output", PathwayPGMConfiguration.getConfig().getNumberOfStates());
        Variable input = variableManager.getVarForName("input", PathwayPGMConfiguration.getConfig().getNumberOfStates());
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(output);
        variables.add(input);
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        edgeTypes.add(FactorEdgeType.OUTPUT);
        edgeTypes.add(FactorEdgeType.INPUT);
        List<Double> values = valueAssigner.generateFactorValues(variables, edgeTypes);
        double[] values1 = new double[values.size()];
        for (int i = 0; i < values.size(); i++)
            values1[i] = values.get(i);
        PathwayPGMConfiguration.getConfig().setCentralDogmaValues(values1);
    }
    
    public void setNamesForEscape(List<String> names) {
        variableManager.setNamesForEscape(names);
    }
    
    /**
     * The main method that should be used to convert a GKInstance pathway into a FactorGraph. 
     * Note: The passed pathway has to have a PathwayDiagram associated.
     * @param pathway
     * @return
     * @throws Exception
     */
    public FactorGraph convertPathway(GKInstance pathway) throws Exception {
        logger.info("Converting pathway " + pathway + "...");
        variableManager.reset();
        // Have to make sure there is a PathwayDiagram associated with the passed pathway
        Collection<GKInstance> referrers = pathway.getReferers(ReactomeJavaConstants.representedPathway);
        if (referrers == null || referrers.size() == 0)
            throw new IllegalArgumentException(pathway + " doesn't have a PathwayDiagram associated!");
        // Get all reactions
        Set<GKInstance> containedReactions = InstanceUtilities.getContainedEvents(pathway);
//        filterToReactions(reactions);
        List<HyperEdge> displayedEdges = ReactomeDataUtilities.getDisplayedEdges(pathway);
        Set<Factor> factors = new HashSet<Factor>();
        PersistenceAdaptor adaptor = pathway.getDbAdaptor();
        // Keep track output nodes for converted Reactions so that
        // we can link an output GKInstance to multiple reactions
        Map<GKInstance, Variable> reactionToOutputVar = new HashMap<GKInstance, Variable>();
        for (Iterator<HyperEdge> it = displayedEdges.iterator(); it.hasNext();) {
            HyperEdge edge = it.next();
            // In this implementation, RenderableInteraction will not be handled since
            // they are actually not annotated in Reactome pathways
            if (edge.getReactomeId() == null) {
                it.remove(); // Remove it for further analysis
                continue;
            }
            GKInstance inst = adaptor.fetchInstance(edge.getReactomeId());
            if (inst == null) {
                // This may be possible for some disease shared pathways.
                logger.error(pathway + "'s diagram has null object for " + edge.getReactomeId());
                it.remove();
                continue; 
            }
            if (!containedReactions.contains(inst)) { // Shared diagrams: disease reactions should be excluded.
                logger.warn(inst + " is not contained by pathway!");
                it.remove();
                continue;
            }
            if (inst.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                reactionHandler.handleReaction(factors, 
                                               inst,
                                               edge,
                                               reactionToOutputVar);
            }
        }
        logger.info("Total reactions converted: " + reactionToOutputVar.size());;
        if (reactionToOutputVar.size() == 0) // Nothing to be done since no outputs have been generated.
            return null;
        linkOutputs(displayedEdges,
                    reactionToOutputVar,
                    factors,
                    adaptor);
        expandHelper.augmentInputs(displayedEdges,
                                   factors, 
                                   adaptor);
        addVariableRoles(displayedEdges, adaptor);
        
        ComponentHelper compHelper = new ComponentHelper();
        compHelper.ensureOneComponent(factors);
        
        if (!isForPerturbation) { // Default we need to load observation files
            // Load data files
            Map<DataType, String> typeToFile = PathwayPGMConfiguration.getConfig().getTypeToEvidenceFile();
            if (typeToFile != null) {
                for (DataType type : typeToFile.keySet()) {
                    String file = typeToFile.get(type);
                    Map<String, Map<String, Float>> sampleToGeneToValue = observationLoader.loadObservationData(file, type);
                    observationLoader.addObservation(sampleToGeneToValue, 
                                                     type,
                                                     variableManager,
                                                     factors);
                }
            }
        }
//        else {
//            // We will create reverse relationships from catalysts to inputs so that quantitative
//            // perturbation to the upstream inputs can be considered in a way when catalysts'
//            // functions are impacted
//            for (GKInstance reaction : reactionToOutputVar.keySet()) {
//                reactionHandler.addReverseForPerturbation(reaction,
//                                                          variableManager,
//                                                          factors);
//            }
//        }
        
        FactorGraph fg = new FactorGraph();
        fg.setName(pathway.toString());
        fg.setFactors(factors);
        fg.validatVariables();
        fg.setIdsInFactors();
        return fg;
    }
    
    /**
     * Link outputs of reactions to output variables converted for reactions. An output GKInstance
     * may be linked to multiple reactions. In this method, such kind of output will be linked to
     * multiple reactions in a single factor.
     * @param edges
     * @param rxtToOuputVar
     */
    private void linkOutputs(List<HyperEdge> edges,
                             Map<GKInstance, Variable> rxtToOuputVar,
                             Set<Factor> factors,
                             PersistenceAdaptor adaptor) throws Exception {
        // Get the map from an output to its linked HyperEdges
        Map<Node, Set<HyperEdge>> outputToEdges = new HashMap<Node, Set<HyperEdge>>();
        for (HyperEdge edge : edges) {
            List<Node> outputs = edge.getOutputNodes();
            for (Node output : outputs) {
                InteractionUtilities.addElementToSet(outputToEdges, 
                                                     output,
                                                     edge);
            }
        }
        // Link outputs to converted output variables using edges labeled as member
        // since we only need one Reaction to generate output.
        List<Variable> outputVariables = new ArrayList<Variable>();
        for (Node output : outputToEdges.keySet()) {
            if (output.getReactomeId() == null)
                continue;
            GKInstance outputInst = adaptor.fetchInstance(output.getReactomeId());
            if (outputInst == null)
                continue;
            Variable outputInstVar = variableManager.getVariable(outputInst);
            if (outputInstVar == null)
                continue;
            Set<HyperEdge> linkedEdges = outputToEdges.get(output);
            outputVariables.clear();
            for (HyperEdge edge : linkedEdges) {
                if (edge.getReactomeId() == null)
                    continue;
                GKInstance rxt = adaptor.fetchInstance(edge.getReactomeId());
                if (rxt == null)
                    continue;
                Variable outputVar = rxtToOuputVar.get(rxt);
                if (outputVar != null)
                    outputVariables.add(outputVar);
            }
            setHandler.handleSetOfVariables(outputVariables,
                                            outputInstVar,
                                            FactorEdgeType.MEMBER,
                                            factors);
        }
    }

    public Map<GKInstance, Variable> getInstToVarMap() {
        return variableManager.getInstToVarMap();
    }
    
    /**
     * Attach the roles of variables into variables for future uses.
     * @param displayedEdges
     */
    private void addVariableRoles(List<HyperEdge> displayedEdges,
                                  PersistenceAdaptor adaptor) throws Exception {
        Map<Variable, Set<VariableRole>> varToRoles = new HashMap<Variable, Set<VariableRole>>();
        Map<GKInstance, Variable> instToVar = getInstToVarMap();
        for (HyperEdge edge : displayedEdges) {
            List<Node> inputs = edge.getInputNodes();
            addVariableRoles(adaptor, varToRoles, instToVar, inputs, VariableRole.INPUT);
            List<Node> outputs = edge.getOutputNodes();
            addVariableRoles(adaptor, varToRoles, instToVar, outputs, VariableRole.OUTPUT);
            List<Node> catalysts = edge.getHelperNodes();
            addVariableRoles(adaptor, varToRoles, instToVar, catalysts, VariableRole.CATALYST);
            List<Node> inhibitors = edge.getInhibitorNodes();
            addVariableRoles(adaptor, varToRoles, instToVar, inhibitors, VariableRole.INHIBITOR);
            List<Node> activators = edge.getActivatorNodes();
            addVariableRoles(adaptor, varToRoles, instToVar, activators, VariableRole.ACTIVATOR);
        }
        for (Variable var : varToRoles.keySet()) {
            Set<VariableRole> roles = varToRoles.get(var);
            var.setProperty("role", InteractionUtilities.joinStringElements(",", roles));
        }
    }

    private void addVariableRoles(PersistenceAdaptor adaptor,
                                  Map<Variable, Set<VariableRole>> varToRoles,
                                  Map<GKInstance, Variable> instToVar,
                                  List<Node> nodes,
                                  VariableRole role) throws Exception {
        for (Node node : nodes) {
            if (node.getReactomeId() == null)
                continue;
            GKInstance inst = adaptor.fetchInstance(node.getReactomeId());
            if (inst == null)
                continue;
            Variable var = instToVar.get(inst);
            if (var == null)
                continue;
            Set<VariableRole> roles = varToRoles.get(var);
            if (roles == null) {
                roles = new HashSet<VariableRole>();
                varToRoles.put(var, roles);
            }
            roles.add(role);
        }
    }
    
    @Test
    public void testConvertPathway() throws Exception {
        FileUtility.initializeLogging();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_59_plus_i", 
                                            "root", 
                                            "macmysql01");
        Long dbId = 2032785L; // YAP1- and WWTR1-stimulated gene expression
        dbId = 1257604L; // PIP3 activates AKT signaling: a pathway has a known complex loop
        // There are multiple activators of this pathway
//        dbId = 452723L; // PathwayDetailedInfo_OSTypes_NoMeta_All_Combined.xlsx
//        dbId = 112310L;
//        dbId = 69620L; // Cell cycle check points
//        dbId = 453279L; // Mitotic G1-G1/S phases
        GKInstance pathway = dba.fetchInstance(dbId);
//        System.out.println("Observations are loaded:");
//        FactorGraph fg = convertPathway(pathway);
//        System.out.println("Total factors: " + fg.getFactors().size());
//        System.out.println("Total variables: " + fg.getVariables().size());
//        long totalMemory = Runtime.getRuntime().totalMemory();
//        long freeMemory = Runtime.getRuntime().freeMemory();
//        long usedMemory = totalMemory - freeMemory;
//        System.out.println("Total used memory: " + usedMemory / (1024 * 1024.0) + " M");
        
        setForPerturbation(true);
        
        // This is only for test: a circular refrence has been created here
        setNamesForEscape(new ReactomePathwayFGRunner().getNameForEscapeList());
        
        System.out.println("\nObservations are not loaded:");
        FactorGraph fg = convertPathway(pathway);
        logger.info("Total factors: " + fg.getFactors().size());
        logger.info("Total variables: " + fg.getVariables().size());
        logger.info("Is it a tree: " + fg.isTree());
//        List<String> varNames = new ArrayList<String>();
//        for (Variable variable : fg.getVariables())
//            varNames.add(variable.getName());
//        Collections.sort(varNames);
//        for (String varName : varNames)
//            System.out.println(varName);
        
        double min = Double.MAX_VALUE;
        for (Factor factor : fg.getFactors()) {
            for (double value : factor.getValues()) {
                if (value > 0.0d && value < min)
                    min = value;
            }
        }
        System.out.println("Minimum value: " + min);
        if (true)
            return;
        
        long totalMemory = Runtime.getRuntime().totalMemory();
        long freeMemory = Runtime.getRuntime().freeMemory();
        long usedMemory = totalMemory - freeMemory;
        System.out.println("Total used memory: " + usedMemory / (1024 * 1024.0) + " M");
        
        // Just for a simple inference
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        lbp.setUseLogSpace(true);
        lbp.setFactorGraph(fg);
//        lbp.runInference();
        
        // For checking
        String fileName = "tmp/" + pathway.getAttributeValue(ReactomeJavaConstants.name) + "_020117.xml";
//        String fileName = "tmp/YAP1_FG.xml";
        FileOutputStream fos = new FileOutputStream(fileName);
        fg.exportFG(fos);
        
        if (true)
            return;
        // Check observations
        List<Observation<Number>> observations = observationLoader.getObservations();
        Map<DataType, ObservationData> typeToData = observationLoader.getLoadedData();
        ObservationRandomizer randomizer = new ObservationRandomizer();
        randomizer.setNumberOfPermutation(10);
        Map<DataType, ObservationFactorHandler> typeToHandler = observationLoader.getDataTypeToObservationFactorHandler();
        List<Observation<Number>> randomObs = randomizer.randomize(observations, 
                                                                   new ArrayList<ObservationData>(typeToData.values()),
                                                                   typeToHandler);
        Observation<Number> obs = observations.get(0);
        System.out.println("One real data: " + obs.getName());
        Observation<Number> rObs = randomObs.get(0);
        System.out.println("One random data: " + rObs.getName());
        for (Variable var : obs.getVariableToAssignment().keySet()) {
            Number state = obs.getVariableToAssignment().get(var);
            Number rState = rObs.getVariableToAssignment().get(var);
//            !!!random data is wrong: geneexp and CNV are the same for the same gene!!!
            System.out.println(var.getName() + "\t" + state + "\t" + rState);
        }
    }
    
}
