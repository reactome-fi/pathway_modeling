/*
 * Created on Oct 15, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.render.HyperEdge;
import org.gk.render.Node;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.common.CentralDogmaHandler;
import org.reactome.factorgraph.common.PGMConfiguration;
import org.reactome.r3.util.ReactomeDataUtilities;

/**
 * This class is used to help expand PhysicalEntity instances used in a Pathway.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class EntityExpandHelper {
    private static final Logger logger = Logger.getLogger(EntityExpandHelper.class);
    private PathwayVariableManager variableManager;
    private VariableSetHandler setHandler;
    // Used to handle central dogma node
    private CentralDogmaHandler dogmaHandler;
    // A flag to indicate if this object is used for perturbation converting
    private boolean isForPerturbation;
    
    /**
     * Default constructor.
     */
    public EntityExpandHelper() {
        dogmaHandler = new CentralDogmaHandler();
        dogmaHandler.setForPerturbation(isForPerturbation);
    }
    
    public boolean isForPerturbation() {
        return isForPerturbation;
    }

    public void setForPerturbation(boolean isForPerturbation) {
        this.isForPerturbation = isForPerturbation;
        dogmaHandler.setForPerturbation(isForPerturbation);
    }

    public CentralDogmaHandler getDogmaHandler() {
        return this.dogmaHandler;
    }
    
    public void setVariableManager(PathwayVariableManager manager) {
        this.variableManager = manager;
    }
    
    public void setVariableSetHandler(VariableSetHandler setHandler) {
        this.setHandler = setHandler;
    }
    
    public void setPGMConfiguration(PGMConfiguration config) {
        this.dogmaHandler.setConfiguration(config);
    }
    
    /**
     * Expand inputs in the pathway.
     * @param edges
     * @param factors
     * @param dba
     * @throws Exception
     */
    public void augmentInputs(List<HyperEdge> edges,
                              Set<Factor> factors,
                              PersistenceAdaptor dba) throws Exception {
        Set<GKInstance> outputInstances = new HashSet<GKInstance>();
        Set<Node> inputs = ReactomeDataUtilities.getNodesForAuguemnt(edges, dba, outputInstances);
        logger.info("Total inputs (inputs, catalysts, and regulators): " + inputs.size());
        
        Set<GKInstance> processed = new HashSet<GKInstance>();
        for (Node node : inputs) {
            if (node.getReactomeId() == null)
                continue; // This should not be possible. But just in case.
            GKInstance input = dba.fetchInstance(node.getReactomeId());
            if (input == null) {
                logger.error(node.getDisplayName() + " with DB_ID = " + node.getReactomeId() + " is not in the database!");
                continue;
            }
            // In case this node has been escaped (e.g. ATP)
            if (!variableManager.isInstanceConverted(input)) {
                continue;
            }
            augumentEntity(input,
                           factors,
                           processed,
                           outputInstances);
        }
    }
    
    /**
     * a protein needs to add DNA and mRNA nodes, and a RNA should add DNA.
     * @param ewas
     * @param factors
     * @param keyToVar
     * @throws Exception
     * TODO: Need to check cases where mRNA or DNA are used directly in reactions. Currently they are not supported.
     */
    private void augmentEWAS(GKInstance ewas,
                             Set<Factor> factors) throws Exception {
//        logger.info("Augment EWAS: " + ewas);
        // Determine what type of EWAS it is
        GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
        if (refEntity == null)
            return; // Maybe during test
        // Don't augument DNA
        if (refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence))
            return ;
        String name = getGeneName(ewas, refEntity);
        if (name == null)
            return;
        // Regardless the type of EWAS (either protein or RNA), add extra EWAS node as product, mRNA, and DNA
        // Add a Protein/RNA node
        Variable ewasVar = variableManager.getVariable(ewas);
        dogmaHandler.addCentralDogmaNodes(ewasVar,
                                          name, 
                                          factors,
                                          variableManager);
    }

    public String getGeneName(GKInstance ewas,
                              GKInstance refEntity) throws Exception {
        String geneName = null;
        if (refEntity.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
            geneName = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
        }
        if (geneName == null)
            geneName = (String) ewas.getAttributeValue(ReactomeJavaConstants.name);
        if (geneName == null)
            geneName = ewas.getDisplayName();
        return geneName;
    }
    
    private void augmentEntitySet(GKInstance set,
                                  Set<Factor> factors,
                                  Set<GKInstance> processed,
                                  Set<GKInstance> outputs) throws Exception {
        logger.info("Augment EntitySet: " + set);
        Set<Variable> memberVars = new HashSet<Variable>();
        List<GKInstance> list = set.getAttributeValuesList(ReactomeJavaConstants.hasMember);
        Set<GKInstance> members = new HashSet<GKInstance>();
        if (list != null && list.size() > 0) {
            for (GKInstance inst : list) {
                Variable var = variableManager.getVariable(inst);
                if (var != null) {
                    memberVars.add(var);
                    members.add(inst);
                }
            }
        }
        if (set.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasCandidate)) {
            list = set.getAttributeValuesList(ReactomeJavaConstants.hasCandidate);
            if (list != null && list.size() > 0) {
                for (GKInstance inst : list) {
                    Variable var = variableManager.getVariable(inst);
                    if (var != null) {
                        memberVars.add(var);
                        members.add(inst);
                    }
                }
            }
        }
        if (memberVars.size() == 0) // Nothing to be done.
            return;
        Variable setVar = variableManager.getVariable(set);
        List<Variable> varList = new ArrayList<Variable>(memberVars);
        setHandler.handleSetOfVariables(varList,
                                        setVar,
                                        FactorEdgeType.MEMBER, 
                                        factors);
        // Need to call recursive
        for (GKInstance inst : members) {
            augumentEntity(inst,
                           factors,
                           processed,
                           outputs);
        }
    }
    
    private void augmentComplex(GKInstance complex,
                                Set<Factor> factors,
                                Set<GKInstance> processed,
                                Set<GKInstance> outputs) throws Exception {
//        logger.info("Augment EntitySet: " + complex);
        Set<Variable> compVars = new HashSet<Variable>();
        List<GKInstance> list = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        Set<GKInstance> comps = new HashSet<GKInstance>();
        if (list != null && list.size() > 0) {
            for (GKInstance inst : list) {
                Variable var = variableManager.getVariable(inst);
                if (var != null) {
                    compVars.add(var);
                    comps.add(inst);
                }
            }
        }
        if (compVars.size() == 0) // Nothing to be done.
            return;
        Variable setVar = variableManager.getVariable(complex);
        List<Variable> varList = new ArrayList<Variable>(compVars);
        setHandler.handleSetOfVariables(varList,
                                        setVar,
                                        FactorEdgeType.COMPLEX, 
                                        factors);
        // Need to call recursive
        for (GKInstance inst : comps) {
            augumentEntity(inst,
                           factors,
                           processed,
                           outputs);
        }
    }
    
    private void augumentEntity(GKInstance inst, 
                              Set<Factor> factors,
                              Set<GKInstance> processed,
                              Set<GKInstance> outputs) throws Exception {
        if (processed.contains(inst) || // This has been processed
            outputs.contains(inst))  // It is used as an output, and should not be processed
            return;
        processed.add(inst);
        if (inst.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
            augmentEWAS(inst, 
                        factors);
        else if (inst.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            augmentEntitySet(inst, 
                             factors,
                             processed,
                             outputs);
        else if (inst.getSchemClass().isa(ReactomeJavaConstants.Complex))
            augmentComplex(inst, 
                           factors,
                           processed,
                           outputs);
    }
    
}
