/*
 * Created on Nov 11, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.gk.schema.InvalidAttributeException;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.CentralDogmaHandler;
import org.reactome.factorgraph.common.PGMConfiguration;

/**
 * This class is used to handle reaction during converting from a GKInstance pathway to a FactorGraph object.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactionHandler {
    // We only need one copy of FactorValueAssigner
    private FactorValueAssigner valueAssigner;
    // Manage variable generation
    private PathwayVariableManager variableManager;
    // For handling a group of variables
    private VariableSetHandler setHandler;
    // for handling gene regulatory interaction
    private EntityExpandHelper expandHelper;
    
    /**
     * Default constructor.
     */
    public ReactionHandler() {
    }
    
    public FactorValueAssigner getValueAssigner() {
        return valueAssigner;
    }

    public void setValueAssigner(FactorValueAssigner valueAssigner) {
        this.valueAssigner = valueAssigner;
    }

    public PathwayVariableManager getVariableManager() {
        return variableManager;
    }

    public void setVariableManager(PathwayVariableManager variableManager) {
        this.variableManager = variableManager;
    }

    public VariableSetHandler getSetHandler() {
        return setHandler;
    }

    public void setSetHandler(VariableSetHandler setHandler) {
        this.setHandler = setHandler;
    }
    
    
    public EntityExpandHelper getExpandHelper() {
        return expandHelper;
    }

    public void setExpandHelper(EntityExpandHelper expandHelper) {
        this.expandHelper = expandHelper;
    }

    /**
     * A special case: we need to add central dogma nodes here so that regulators can be attached
     * to factor: GENE_DNA -> GENE_mRNA directly at this stage.
     * @param factors
     * @param reaction
     * @param rxtToOutputVar
     * @throws Exception
     */
    private void handleGeneRegulatoryReaction(Set<Factor> factors,
                                              GKInstance reaction, 
                                              Map<GKInstance, Variable> rxtToOutputVar) throws Exception {
        GKInstance output = (GKInstance) reaction.getAttributeValue(ReactomeJavaConstants.output);
        if (output == null || !output.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
            return; // Nothing to be done. Note: the output should be an EWAS in order to make a meaningful gene regulatory interaction.
        GKInstance refEntity = (GKInstance) output.getAttributeValue(ReactomeJavaConstants.referenceEntity);
        // Don't augument DNA
        if (refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence))
            return ;
        String name = expandHelper.getGeneName(output,
                                               refEntity);
        if (name == null)
            return;
        CentralDogmaHandler dogmaHandler = expandHelper.getDogmaHandler();
        PGMConfiguration configuration = dogmaHandler.getConfiguration();
        // Add a protein node
        int varStates = configuration.getNumberOfStates();
        String protein = name + "_" + PGMConfiguration.protein;
        Variable proteinVar = variableManager.getVarForName(protein, varStates);
        
        Variable mRNAVar = variableManager.getVarForName(name + "_" + PGMConfiguration.mRNA, varStates);
        
        Variable dnaVar = variableManager.getVarForName(name + "_" + PGMConfiguration.DNA, varStates);
        
        Variable outputVar = variableManager.getVarForName(reaction.getDBID() + "_" + FactorEdgeType.OUTPUT, varStates);
        outputVar.setProperty(ReactomeJavaConstants.DB_ID,
                              reaction.getDBID() + "");
        rxtToOutputVar.put(reaction, outputVar);
        
        // Create factors among these nodes
        // Protein to the reaction output
        dogmaHandler.createCentralDogmaFactor(proteinVar,
                                              outputVar,
                                              factors,
                                              null);
        // mRNA to protein
        dogmaHandler.createCentralDogmaFactor(mRNAVar,
                                              proteinVar,
                                              factors,
                                              null);
        // We need to attach regulators to relationship between DNA and mRNA
        // and treat them in a single Reactome reaction.
        List<Variable> rxtVars = new ArrayList<Variable>();
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        rxtVars.add(dnaVar);
        edgeTypes.add(FactorEdgeType.INPUT);
        handleReactionRegulations(reaction, 
                                  factors,
                                  rxtVars,
                                  edgeTypes);
        // Don't forget to add mRNA
        rxtVars.add(mRNAVar);
        edgeTypes.add(FactorEdgeType.OUTPUT);
        Factor factor = variableManager.createFactor(rxtVars, 
                                                     edgeTypes,
                                                     valueAssigner, 
                                                     reaction.toString());
        factor.setProperty(ReactomeJavaConstants.DB_ID, 
                           reaction.getDBID() + "");
        factors.add(factor);
    }
    
    public void addReverseForPerturbation(GKInstance reaction,
                                          PathwayVariableManager varManager,
                                          Set<Factor> factors) throws Exception {
        // Check if there is a catalyst existing
        GKInstance cas = (GKInstance) reaction.getAttributeValue(ReactomeJavaConstants.catalystActivity);
        if (cas == null)
            return; // Do nothing
        GKInstance catalyst = (GKInstance) cas.getAttributeValue(ReactomeJavaConstants.physicalEntity);
        if (catalyst == null)
            return; // Do nothing
        Variable catVar = varManager.getVariable(catalyst);
        if (catVar == null)
            throw new IllegalStateException("Cannot find converted Variable for " + catalyst);
        // Get inputs variables first. We may not be able to get all variables
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (inputs == null || inputs.size() == 0)
            return;
        Set<Variable> inputVars = new HashSet<Variable>();
        for (GKInstance input : inputs) {
            Variable inputVar = varManager.getVariable(input);
            if (inputVar == null)
                continue; // Some entities (e.g H2O) are escaped
            inputVars.add(inputVar);
        }
        if (inputVars.size() == 0)
            return; // Nothing to be done
        
        // As with the forward reaction, we will add an input variable to tune the effects a little bit
        Variable reverseInputVar = varManager.getVarForName("reverse_input_" + reaction.getDBID(),
                                                            PathwayPGMConfiguration.getConfig().getNumberOfStates());
        List<Variable> variables = new ArrayList<Variable>();
        variables.add(catVar);
        variables.add(reverseInputVar);
        List<FactorEdgeType> edgeTypes = new ArrayList<FactorEdgeType>();
        edgeTypes.add(FactorEdgeType.INPUT);
        edgeTypes.add(FactorEdgeType.OUTPUT);
        
        Factor factor = variableManager.createFactor(variables, 
                                                     edgeTypes,
                                                     valueAssigner, 
                                                     reverseInputVar.getName());
        factors.add(factor);
        
        for (Variable inputVar : inputVars) {
            variables.clear();
            variables.add(reverseInputVar);
            variables.add(inputVar);
            factor = variableManager.createFactor(variables, 
                                                  edgeTypes,
                                                  valueAssigner, 
                                                  "reverse_" + reaction.getDBID() + "_input_" + inputVar.getName());
            factors.add(factor);
        }
    }

    public void handleReaction(Set<Factor> factors,
                                GKInstance reaction,
                                HyperEdge edge,
                                Map<GKInstance, Variable> rxtToOutputVar) throws Exception {
        if (isGeneRegulatoryInteraction(edge, 
                                        reaction.getDbAdaptor())) {
            handleGeneRegulatoryReaction(factors, reaction, rxtToOutputVar);
            return;
        }
        // Variables related to this reaction instance directly
        List<Variable> rxtVariables = new ArrayList<Variable>();
        // Edge types for the above variables in the same order
        List<FactorEdgeType> rxtEdgeTypes = new ArrayList<FactorEdgeType>();
        handleInputs(reaction, 
                     factors,
                     rxtVariables,
                     rxtEdgeTypes);
        handleCatalystActivity(reaction, 
                               rxtVariables, 
                               rxtEdgeTypes);
        handleReactionRegulations(reaction,
                                  factors, 
                                  rxtVariables, 
                                  rxtEdgeTypes);
        // Need to handle outputs
        // Create a new output_node
        handleReactionOutput(reaction,
                           rxtToOutputVar,
                           rxtVariables, 
                           rxtEdgeTypes);
        Factor factor = variableManager.createFactor(rxtVariables, 
                                                     rxtEdgeTypes,
                                                     valueAssigner, 
                                                     reaction.toString());
        factor.setProperty(ReactomeJavaConstants.DB_ID, 
                           reaction.getDBID() + "");
        factors.add(factor);
    }

    private void handleInputs(GKInstance reaction,
                              Set<Factor> factors, 
                              List<Variable> rxtVariables,
                              List<FactorEdgeType> rxtEdgeTypes) throws InvalidAttributeException, Exception {
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (inputs != null && inputs.size() > 0) {
            handleReactionParticipants(reaction,
                                       inputs,
                                       FactorEdgeType.INPUT,
                                       reaction.getDBID() + "_" + FactorEdgeType.INPUT,
                                       rxtVariables,
                                       rxtEdgeTypes,
                                       factors);
        }
    }

    private void handleCatalystActivity(GKInstance reaction,
                                        List<Variable> rxtVariables,
                                        List<FactorEdgeType> rxtEdgeTypes)
            throws InvalidAttributeException, Exception {
        // Handle catalyst
        GKInstance cas = (GKInstance) reaction.getAttributeValue(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            GKInstance catalyst = (GKInstance) cas.getAttributeValue(ReactomeJavaConstants.physicalEntity);
            if (catalyst != null) {
                // There will be only one catalyst
                Variable var = variableManager.getVariable(catalyst);
                // If it is a self-catalyst, no need to add this variable again.
                if (var != null && !rxtVariables.contains(var)) {
                    rxtVariables.add(var);
                    rxtEdgeTypes.add(FactorEdgeType.CATALYST);
                }
            }
        }
    }

    private void handleReactionOutput(GKInstance reaction,
                                      Map<GKInstance, Variable> rxtToOutputVar,
                                      List<Variable> rxtVariables,
                                      List<FactorEdgeType> rxtEdgeTypes) {
        Variable outputVar = variableManager.getVarForName(reaction.getDBID() + "_" + FactorEdgeType.OUTPUT, 
                                                           PathwayPGMConfiguration.getConfig().getNumberOfStates());
        outputVar.setProperty(ReactomeJavaConstants.DB_ID, 
                              reaction.getDBID() + "");
        rxtToOutputVar.put(reaction, outputVar);
        rxtVariables.add(outputVar);
        rxtEdgeTypes.add(FactorEdgeType.OUTPUT);
    }

    private void handleReactionRegulations(GKInstance reaction,
                                           Set<Factor> factors,
                                           List<Variable> rxtVariables,
                                           List<FactorEdgeType> rxtEdgeTypes) throws Exception, InvalidAttributeException {
        // Handle regulators
        Collection<GKInstance> regulations = InstanceUtilities.getRegulations(reaction);
        List<GKInstance> activators = new ArrayList<GKInstance>();
        List<GKInstance> inhibitors = new ArrayList<GKInstance>();
        if (regulations != null && regulations.size() > 0) {
            for (GKInstance regulation : regulations) {
                GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
                if (regulator != null && regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                    if (regulation.getSchemClass().isa(ReactomeJavaConstants.PositiveRegulation))
                        activators.add(regulator);
                    else if (regulation.getSchemClass().isa(ReactomeJavaConstants.NegativeRegulation))
                        inhibitors.add(regulator);
                }
            }
        }
        // Handle activators
        if (activators.size() > 0) {
            handleReactionParticipants(reaction,
                                       activators, 
                                       FactorEdgeType.ACTIVATOR, // We are going to label edge with max, so use "MEMBER" 
                                       reaction.getDBID() + "_" + FactorEdgeType.ACTIVATOR,
                                       rxtVariables,
                                       rxtEdgeTypes, 
                                       factors);
        }
        if (inhibitors.size() > 0)
            handleReactionParticipants(reaction,
                                       inhibitors, 
                                       FactorEdgeType.INHIBITOR, // same as activators 
                                       reaction.getDBID() + "_" + FactorEdgeType.INHIBITOR,
                                       rxtVariables,
                                       rxtEdgeTypes, 
                                       factors);
    }

    private void handleReactionParticipants(GKInstance reaction,
                                            List<GKInstance> instances,
                                            FactorEdgeType edgeType,
                                            String variableName,
                                            List<Variable> rxtVariables,
                                            List<FactorEdgeType> rxtEdgeTypes,
                                            Set<Factor> factors) throws Exception {
        if (instances.size() == 1) {
            // Just use input without creating an accessory input node
            GKInstance inst = instances.get(0);
            Variable var = variableManager.getVariable(inst);
            // Don't add the same variable twice, which will break the inference algorithm.
            //TODO: If an input is used as an inhibitor, the inhibitor effect will not be considered
            // by this check!
            if (var != null && !rxtVariables.contains(var)) {
                rxtVariables.add(var);
                rxtEdgeTypes.add(edgeType);
            }
        }
        else {
            // Need to create an input node
            Variable var1 = variableManager.getVarForName(variableName,
                                                          PathwayPGMConfiguration.getConfig().getNumberOfStates());
            var1.setProperty(ReactomeJavaConstants.DB_ID, 
                             reaction.getDBID() + "");
            Variable var = var1;
            rxtVariables.add(var);
            rxtEdgeTypes.add(edgeType);
            handleListOfInstances(instances,
                                  var,
                                  FactorEdgeType.mapGroupTypeToSingleType(edgeType),
                                  factors);
        }
    }
    
    /**
     * Generate Factors between a list of GKInstances and the passed target Variable object. If
     * there are too many instance in the list (aka > configured maximum number), accessory nodes
     * will be created.
     * @param instances
     * @param targetVar
     * @param edgeType
     */
    private void handleListOfInstances(List<GKInstance> instances,
                                       Variable targetVar,
                                       FactorEdgeType edgeType,
                                       Set<Factor> factors) throws Exception {
        // Use set to avoid duplication
        Set<Variable> variables = new HashSet<Variable>();
        for (GKInstance inst : instances) {
            Variable var = variableManager.getVariable(inst);
            if (var != null)
                variables.add(var);
        }
        if (variables.size() == 0)
            return; // There is no need
        List<Variable> varList = new ArrayList<Variable>(variables);
        setHandler.handleSetOfVariables(varList, 
                                        targetVar, 
                                        edgeType, 
                                        factors);
    }
    
    
    private boolean isGeneRegulatoryInteraction(HyperEdge edge,
                                                PersistenceAdaptor dba) throws Exception {
        if (edge.getReactomeId() == null)
            return false;
        List<Node> list = edge.getInputNodes();
        if (list != null && list.size() > 0)
            return false; // Should not be considered as gene regulatory interaction.
        list = edge.getOutputNodes();
        if (list != null && list.size() != 1)
            return false; // Only one output only
        GKInstance rxt = dba.fetchInstance(edge.getReactomeId());
        if (!rxt.getSchemClass().isa(ReactomeJavaConstants.BlackBoxEvent))
            return false;
        // A gene regulatory interaction should be regulated by others
        Collection<GKInstance> regulators = InstanceUtilities.getRegulations(rxt);
        if (regulators == null || regulators.size() == 0)
            return false;
        return true;
    }
    
}
