/*
 * Created on Jan 29, 2015
 *
 */
package org.reactome.factorgraph.common;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;
import org.reactome.factorgraph.common.ObservationFileLoader.ObservationData;

/**
 * This class is used to generate a list of randomized Observation<Float> objects from
 * a passed list of Observation<Float> objects.
 * @author gwu
 *
 */
public class ObservationRandomizer {
    // Used to label the generated random samples
    private String randomSamplePrefix = "Random";
    private int numberOfPermutation = 1000;
    
    /**
     * Default constructor.
     */
    public ObservationRandomizer() {
    }
    
    public String getRandomSamplePrefix() {
        return randomSamplePrefix;
    }

    public void setRandomSamplePrefix(String randomSamplePrefix) {
        this.randomSamplePrefix = randomSamplePrefix;
    }

    public int getNumberOfPermutation() {
        return numberOfPermutation;
    }

    public void setNumberOfPermutation(int numberOfPermutation) {
        this.numberOfPermutation = numberOfPermutation;
    }
    
    /**
     * Generate a list of randomized Observtaion objects from the passed observations object.
     * This method should be used for Observation objects containing a large number of observation
     * variables to avoid a large bias (e.g. Observations generated for a FI PGM model).
     * @param observations
     * @return
     */
    public List<Observation<Number>> createRandomObservations(List<Observation<Number>> observations) {
        List<Observation<Number>> randomized = new ArrayList<Observation<Number>>();
        // Create a random gene list 
        // We will use samples having all data types to avoid bias in the final running
        // For example, if two data types, mutation and CNV, are used, but CNV samples have one
        // time more than mutation. So some of random genes may not get mutation observation in
        // the same random sample, which increase the final impact score since it is treated as
        // same probability of two states.
        ObservationHelper helper = new ObservationHelper();
        // Make a copy since we are going to change its content
        observations = new ArrayList<Observation<Number>>(observations);
        helper.filterObservationsToHaveSharedDataTypes(observations);
        Set<DataType> dataTypes = helper.getDataTypesFromObservations(observations);
        Set<Variable> allVariables = getAllVariables(observations);
        Set<String> allGenes = getGenesFromVariables(allVariables, dataTypes);
        // Make sure null should not be returned during random drawing
        Set<String> sharedGenes = new HashSet<String>(allGenes);
        // Will keep genes having observed data in all types
        // This may bias the distribution a little bit. However, it seems
        // the best we can do right now!
        filterGenesToAllDataTypes(sharedGenes, observations, dataTypes);
        List<String> sharedGeneList = new ArrayList<String>(sharedGenes);
        RandomData randomizer = new RandomDataImpl();
        int randomIndex;
        
        // For quick performance
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (Variable var : allVariables)
            nameToVar.put(var.getName(), var);
        
        for (int i = 0; i < numberOfPermutation; i++) {
            Observation<Number> randomObs = new Observation<Number>();
            String randomSample = randomSamplePrefix + i;
            randomObs.setName(randomSample);
            List<VariableAssignment<Number>> randomVarAssgns= new ArrayList<VariableAssignment<Number>>();
            for (String randomGene : allGenes) {
                // Pick up a random sample first
                randomIndex = randomizer.nextInt(0, observations.size() - 1);
                Observation<Number> realObs = observations.get(randomIndex);
                // Pick up a random variable
                randomIndex = randomizer.nextInt(0, sharedGeneList.size() - 1);
                String realGene = sharedGeneList.get(randomIndex);
                // We will get whatever there for realGene in realObs
                for (DataType dataType : dataTypes) {
                    Variable randomVar = nameToVar.get(randomGene + "_" + dataType);
                    if (randomVar == null)
                        continue; // Cannot find this variable in the observation. Just escape it!
                    String varName = realGene + "_" + dataType;
                    Variable var = nameToVar.get(varName);
                    if  (var == null)
                        throw new IllegalStateException("Cannot find variable, " + varName);
                    VariableAssignment<Number> varAssgn = realObs.getVariableAssignment(var);
                    if (varAssgn == null)
                        throw new IllegalStateException("Cannot find value for variable, " + var.getName());
                    // We need to make a copy of this value
                    VariableAssignment<Number> randomVarAssgn = new VariableAssignment<Number>();
                    randomVarAssgn.setVariable(randomVar);
                    randomVarAssgn.setAssignment(varAssgn.getAssignment());
                    randomVarAssgn.setDistribution(varAssgn.getDistribution());
                    randomVarAssgns.add(randomVarAssgn);
                }
            }
            randomObs.setVariableAssignments(randomVarAssgns);
            randomized.add(randomObs);
        }
        return randomized;
    }
    
    private Set<String> getGenesFromVariables(Collection<Variable> variables,
                                              Collection<DataType> dataTypes) {
        Set<String> genes = new HashSet<String>();
        // Get a little bit gain of performance
        Set<String> dataTypePostfix = new HashSet<String>();
        for (DataType dataType : dataTypes)
            dataTypePostfix.add("_" + dataType);
        for (Variable var : variables) {
            String name = var.getName();
            for (String dataType : dataTypePostfix) {
                if (name.endsWith(dataType)) {
                    // Should avoid case like KCNV1 for CNV
                    int index = name.lastIndexOf(dataType.toString());
                    genes.add(name.substring(0, index));
                    break;
                }
            }
        }
        return genes;
    }
    
    private Set<Variable> getAllVariables(List<Observation<Number>> observations) {
        Set<Variable> variables = new HashSet<Variable>();
        for (Observation<Number> obs : observations) {
            List<VariableAssignment<Number>> varAssgns = obs.getVariableAssignments();
            for (VariableAssignment<Number> varAssgn : varAssgns)
                variables.add(varAssgn.getVariable());
        }
        return variables;
    }
    
    private void filterGenesToAllDataTypes(Collection<String> genes,
                                           List<Observation<Number>> observations,
                                           Collection<DataType> dataTypes) {
        Map<String, Variable> nameToVar = generateNameToVariableMap(observations);
        boolean isRemoved = false;
        for (Iterator<String> it = genes.iterator(); it.hasNext();) {
            String gene = it.next();
            for (DataType dataType : dataTypes) {
                Variable var = nameToVar.get(gene + "_" + dataType);
                if (var == null) {
                    it.remove();
                    break;
                }
                // Make sure this variable is in all observation
                isRemoved = false;
                for (Observation<Number> obs : observations) {
                    VariableAssignment<Number> assgn = obs.getVariableAssignment(var);
                    if (assgn == null) {
                        it.remove();
                        isRemoved = true;
                        break;
                    }
                }
                if (isRemoved)
                    break;
            }
        }
    }

    /**
     * The main method for applying randomization for a list of Observation<Float> objects.
     * @param observations a list of variables
     * @param observationData the real observation data.
     * @return
     */
    public List<Observation<Number>> randomize(List<Observation<Number>> observations,
                                               List<ObservationData> observationData,
                                               Map<DataType, ObservationFactorHandler> typeToFactorHandler) {
        List<ObservationData> randomData = randomize(observationData);
        return createRandomObservations(observations,
                                        randomData,
                                        typeToFactorHandler);
    }

    /**
     * Create a list of Observation<Number> objects from the passed randomized observation data.
     * @param observations a list of Observation<Number> objects used as templates for randomized observations.
     * @param randomData a list of ObservationData.
     * @param dataTypeToHandler a map from dataType to factorHandler, which should be used during real data loading.
     * @return
     */
    public List<Observation<Number>> createRandomObservations(List<Observation<Number>> observations,
                                                              List<ObservationData> randomData,
                                                              Map<DataType, ObservationFactorHandler> dataTypeToHandler) {
        Map<String, Variable> nameToVar = generateNameToVariableMap(observations);
        List<Observation<Number>> randomized = new ArrayList<Observation<Number>>();
        List<String> randomSamples = getSamplesFromObservations(randomData);
        int index = 0;
        String geneName = null;
        for (String randomSample : randomSamples) {
            Observation<Number> randomObs = new Observation<Number>();
            randomObs.setName(randomSample);
            randomized.add(randomObs);
            List<VariableAssignment<Number>> varAssgns = new ArrayList<VariableAssignment<Number>>();
            for (int i = 0; i < randomData.size(); i++) {
                ObservationData randomData1 = randomData.get(i);
                DataType dataType = randomData1.getDataType();
                ObservationFactorHandler handler = dataTypeToHandler.get(dataType);
                Map<String, Map<String, Float>> randomSampleToGeneToState = randomData1.getSampleToGeneToValue();
                Map<String, Float> randomGeneToState = randomSampleToGeneToState.get(randomSample);
                // Use a smaller size of nameToVar to increase the performance.
                for (String varName : nameToVar.keySet()) {
                    // Get gene
                    index = varName.indexOf("_");
                    geneName = varName.substring(0, index);
                    String typeName = varName.substring(index + 1);
                    if (!typeName.equals(dataType.toString()))
                        continue;
                    Float value = randomGeneToState.get(geneName);
                    if (value == null)
                        continue;
                    // Distribution for EmpiricalFactorHandler should be handled automatically
                    VariableAssignment<Number> varAssgn = handler.parseValue(value.doubleValue(),
                                                                             dataType,
                                                                             nameToVar.get(varName));
                    varAssgns.add(varAssgn);
                }
            }
            randomObs.setVariableAssignments(varAssgns);
        }
        return randomized;
    }
    
    private Map<String, Variable> generateNameToVariableMap(List<Observation<Number>> observations) {
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (Observation<Number> observation : observations) {
            for (Variable var : observation.getVariableToAssignment().keySet()) {
                if (nameToVar.containsKey(var.getName()))
                    continue;
                nameToVar.put(var.getName(), var);
            }
        }
        return nameToVar;
    }
    
    /**
     * Randomize a list of processed observation data in a tuple of genes. So that all data types for one gene
     * is kept and randomly assigned to another gene to generate a randomized data set.
     * @param sampleToGeneToStateList
     * @return
     */
    public List<ObservationData> randomize(List<ObservationData> observationData) {
        List<ObservationData> randomData = new ArrayList<ObservationFileLoader.ObservationData>();
        // Create a random gene list 
        List<String> allGenes = getAllGenesFromObservations(observationData);
        List<String> sharedGenes = getSharedGenesFromObservation(observationData);
        if (sharedGenes.size() == 0)
            return randomData;
        // Create a random sample list
        List<String> sampleList = getSamplesFromObservations(observationData);
        if (sampleList.size() == 0)
            return randomData;
        // Create random samples
        // Empty value first
        for (ObservationData data : observationData) {
            ObservationData randomData1 = new ObservationData();
            randomData1.setDataType(data.getDataType());
            Map<String, Map<String, Float>> sampleToGeneToState = new HashMap<String, Map<String, Float>>();
            randomData1.setSampleToGeneToValue(sampleToGeneToState);
            randomData.add(randomData1);
        }
        RandomData randomizer = new RandomDataImpl();
        int randomIndex;
        // Starting generating random values in a linked tuple way
        for (int i = 0; i < numberOfPermutation; i++) {
            String sampleName = randomSamplePrefix + i; // A sample in the generated random data set
            for (int k = 0; k < allGenes.size(); k++) {
                String gene = allGenes.get(k); // Gene in the random sample
                // Get a random gene from the shared gene list in order to maximize
                // the chance to get the data set.
                randomIndex = randomizer.nextInt(0, sharedGenes.size() - 1);
                String randomGene = sharedGenes.get(randomIndex);
                // Get a random sample
                randomIndex = randomizer.nextInt(0, sampleList.size() - 1);
                String randomSample = sampleList.get(randomIndex);
                // Assign states for randomGene in randomSample 
                for (int j = 0; j < observationData.size(); j++) {
                    ObservationData realData = observationData.get(j);
                    Map<String, Map<String, Float>> realSampleToGeneToState = realData.getSampleToGeneToValue();
                    Map<String, Float> realGeneToState = realSampleToGeneToState.get(randomSample);
                    Float realState = realGeneToState.get(randomGene);
                    ObservationData randomData1 = randomData.get(j);
                    Map<String, Map<String, Float>> randomSampleToGeneToState = randomData1.getSampleToGeneToValue();
                    Map<String, Float> randomGeneToState = randomSampleToGeneToState.get(sampleName);
                    if (randomGeneToState == null) {
                        randomGeneToState = new HashMap<String, Float>();
                        randomSampleToGeneToState.put(sampleName, randomGeneToState);
                    }
                    randomGeneToState.put(gene, realState);
                    // Since we used the same random gene and same random sample for all evidence maps,
                    // the correlation among these evidences should be kept..
                }
            }
        }
        return randomData;
    }
    
    /**
     * Get shared samples in all data set.
     * @param sampleToGeneToStateList
     * @return
     */
    private List<String> getSamplesFromObservations(List<ObservationData> observationData) {
        Set<String> samples = null;
        for (ObservationData data : observationData) {
            Map<String, Map<String, Float>> sampleToGeneToState = data.getSampleToGeneToValue();
            if (samples == null)
                samples = new HashSet<String>(sampleToGeneToState.keySet());
            else
                samples.retainAll(sampleToGeneToState.keySet());
        }
        if (samples == null)
            return new ArrayList<String>();
        return new ArrayList<String>(samples);
    }
    
    /**
     * Get all genes touched in the observation data.
     * @param observations
     * @return
     */
    private List<String> getAllGenesFromObservations(List<ObservationData> observationData) {
        Set<String> genes = new HashSet<String>();
        for (ObservationData data : observationData) {
            Map<String, Map<String, Float>> sampleToGeneToState = data.getSampleToGeneToValue();
            for (String sample : sampleToGeneToState.keySet()) {
                Map<String, Float> geneToState = sampleToGeneToState.get(sample);
                genes.addAll(geneToState.keySet());
            }
        }
        return new ArrayList<String>(genes);
    }
    
    /**
     * Get shared genes from all observations.
     * @param observationData
     * @return
     */
    private List<String> getSharedGenesFromObservation(List<ObservationData> observationData) {
        Set<String> genes = null;
        for (ObservationData data : observationData) {
            Map<String, Map<String, Float>> sampleToGeneToState = data.getSampleToGeneToValue();
            for (String sample : sampleToGeneToState.keySet()) {
                Map<String, Float> geneToState = sampleToGeneToState.get(sample);
                if (genes == null)
                    genes = new HashSet<String>(geneToState.keySet());
                else
                    genes.retainAll(geneToState.keySet());
            }
        }
        if (genes == null)
            return new ArrayList<String>();
        return new ArrayList<String>(genes);
    }
    
}
