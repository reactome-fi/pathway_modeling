/*
 * Created on May 27, 2014
 *
 */
package org.reactome.factorgraph.common;

import java.util.Collection;
import java.util.Map;

import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.SharedEMFactors;
import org.reactome.factorgraph.Variable;

/**
 * This class is used to add central dogma nodes.
 * @author gwu
 *
 */
public class CentralDogmaHandler {    
    private PGMConfiguration configuration;
    // Check to indicate this is for perturbation studies.
    private boolean isForPerturbation;
    
    /**
     * Default constructor.
     */
    public CentralDogmaHandler() {
    }
    
    public boolean isForPerturbation() {
        return isForPerturbation;
    }

    public void setForPerturbation(boolean isForPerturbation) {
        this.isForPerturbation = isForPerturbation;
    }

    public PGMConfiguration getConfiguration() {
        return configuration;
    }

    public void setConfiguration(PGMConfiguration configuration) {
        this.configuration = configuration;
    }

    /**
     * This method is used to create three central dogma nodes.
     * @param geneVar
     * @param factors
     */
    public void addCentralDogmaNodes(Variable targetVar,
                                     String name,
                                     Collection<Factor> factors,
                                     VariableManager varManager,
                                     Map<String, SharedEMFactors> typeToSharedFactors) {
        Map<String, Variable> nameToVar = varManager.getNameToVar();
        Variable proteinVar = getProteinVar(name, nameToVar);
        if (proteinVar != null) {
            // Central dogma for this gene has been added, there is no need to do.
            // However, we may have to link this proteinVar to targetVar
            // The relationship between fiNode and protein is: protein -> fiNode
            createCentralDogmaFactor(proteinVar, 
                                     targetVar,
                                     factors,
                                     getSharedEMFactors(typeToSharedFactors, 
                                                        PGMConfiguration.protein));
            return;
        }
        // Add a protein node
        String protein = name + "_" + PGMConfiguration.protein;
        proteinVar = varManager.getVarForName(protein, configuration.getNumberOfStates());
        // The relationship between fiNode and protein is: protein -> fiNode
        createCentralDogmaFactor(proteinVar, 
                                 targetVar,
                                 factors,
                                 getSharedEMFactors(typeToSharedFactors, 
                                                    PGMConfiguration.protein));
        // Only protein variables are needed for protein-based perturbation studies
        if (isForPerturbation)
            return;
        
        Variable mRNAVar = varManager.getVarForName(name + "_" + PGMConfiguration.mRNA,
                                                    configuration.getNumberOfStates());
        createCentralDogmaFactor(mRNAVar, 
                                 proteinVar,
                                 factors,
                                 getSharedEMFactors(typeToSharedFactors, PGMConfiguration.mRNA));
        
        Variable dnaVar = varManager.getVarForName(name + "_" + PGMConfiguration.DNA,
                                                   configuration.getNumberOfStates());
        createCentralDogmaFactor(dnaVar,
                                 mRNAVar,
                                 factors,
                                 getSharedEMFactors(typeToSharedFactors, PGMConfiguration.DNA));
    }
    
    private SharedEMFactors getSharedEMFactors(Map<String, SharedEMFactors> typeToSharedFactors,
                                               String type) {
        if (typeToSharedFactors == null)
            return null;
        SharedEMFactors sharedFactors = typeToSharedFactors.get(type);
        if (sharedFactors == null) {
            sharedFactors = new SharedEMFactors();
            typeToSharedFactors.put(type, sharedFactors);
        }
        return sharedFactors;
    }
    
    public void addCentralDogmaNodes(Variable targetVar,
                                     String name,
                                     Collection<Factor> factors,
                                     VariableManager varManager) {
        addCentralDogmaNodes(targetVar, 
                             name, 
                             factors, 
                             varManager, 
                             null);
    }
    
    public Variable getAnchorVar(String name,
                                 VariableManager varManager,
                                 DataType dataType) {
        Variable parent = null;
        Map<String, Variable> nameToVar = varManager.getNameToVar();
        switch (dataType) {
            case CNV :
                parent = getDNAVar(name, nameToVar);
                break;
            case mRNA_EXP :
                parent = getmRNAVar(name, nameToVar);
                break;
            case Methylation :
                parent = getDNAVar(name, nameToVar);
                break;
            case miRNA :
                parent = getmRNAVar(name, nameToVar);
                break;
            case Mutation :
                parent = getProteinVar(name, nameToVar);
                break;
        }
        if (parent == null)
            throw new IllegalArgumentException(dataType + " has not been supported yet: no central dogma node can be assigned for " + name);
        return parent;
    }
    
    private Variable getDNAVar(String name,
                              Map<String, Variable> nameToVar) {
        return nameToVar.get(name + "_" + PGMConfiguration.DNA);
    }
    
    private Variable getmRNAVar(String name,
                               Map<String, Variable> nameToVar) {
        return nameToVar.get(name + "_" + PGMConfiguration.mRNA);
    }
    
    private Variable getProteinVar(String name, 
                                  Map<String, Variable> nameToVar) {
        return nameToVar.get(name + "_" + PGMConfiguration.protein);
    }
    
    /**
     * Check if a central dogma node has been added.
     * @param name
     * @param nameToVar
     * @return
     */
    public boolean isCentralDogmaAdded(String name,
                                       VariableManager varManager) {
        // Use DNA to check if central dogma nodes have been added since proteinVar
        // may not be needed for miRNA
        Variable var = getDNAVar(name, varManager.getNameToVar());
        return var != null;
    }
    
    public void createCentralDogmaFactor(Variable from,
                                         Variable to,
                                         Collection<Factor> factors,
                                         SharedEMFactors sharedEMFactors) {
        double[] values = configuration.getCentralDogmaValues();
        if (values == null)
            throw new IllegalStateException("Values for central dogma factors have not be configured!");
        Factor factor = new Factor(from, to, values);
        factors.add(factor);
    }
    
}
