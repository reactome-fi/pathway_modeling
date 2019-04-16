/*
 * Created on May 27, 2014
 *
 */
package org.reactome.factorgraph.common;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.factorgraph.EMFactor;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Observation;
import org.reactome.factorgraph.SharedEMFactors;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;
import org.reactome.r3.util.FileUtility;

/**
 * A class that is used to load an observation data file.
 * @author gwu
 *
 */
public class ObservationFileLoader {
    // Loaded evidences
    // Use Number to support multiple types
    private Map<String, Observation<Number>> sampleToObservation;
    // A helper object
    protected CentralDogmaHandler dogmaHandler;
    // Cached loaded values to avoid multiple loading, which is very slow for large data sets
    private Map<DataType, ObservationData> typeToData;
    // To handle observation file factor for different data types
    private Map<DataType, ObservationFactorHandler> typeToFactorHandler;
    // To support EM learning
    private Map<DataType, SharedEMFactors> typeToSharedFactors;
    
    /**
     * The default constructor. Since the loaded evidences are cached in the memory, an object of
     * this class should be used for one calculation only.
     */
    public ObservationFileLoader() {
        dogmaHandler = new CentralDogmaHandler();
        sampleToObservation = new HashMap<String, Observation<Number>>();
        typeToData = new HashMap<DataType, ObservationData>();
    }
    
    public void setObservationFactorHandler(DataType dataType,
                                            ObservationFactorHandler factorHandler) {
        if (typeToFactorHandler == null)
            typeToFactorHandler = new HashMap<DataType, ObservationFactorHandler>();
        typeToFactorHandler.put(dataType, 
                                factorHandler);
        factorHandler.setConfiguration(dogmaHandler.getConfiguration());
    }
    
    public ObservationFactorHandler getObservationFactorHandler(DataType dataType) {
        if (typeToFactorHandler != null)
            return typeToFactorHandler.get(dataType);
        return null;
    }
    
    /**
     * Get the cached map from DataType to ObservationFactorHandler. In the current implementation,
     * there is at most one ObservationFactorHandler for one DataType. However, ObservationFactorHandler
     * may be shared among multiple DataTypes.
     * @return
     */
    public Map<DataType, ObservationFactorHandler> getDataTypeToObservationFactorHandler() {
        return this.typeToFactorHandler;
    }
    
    /**
     * Set the configuration used for this running. There should be only one copy of PGMConfiguration object
     * in the whole application running.
     * @param config
     */
    public void setPGMConfiguration(PGMConfiguration config) {
        dogmaHandler.setConfiguration(config);
        if (typeToFactorHandler != null) {
            for (ObservationFactorHandler factorHandler : typeToFactorHandler.values())
                factorHandler.setConfiguration(config);
        }
    }
    
    /**
     * Make a copy of Observation<Float> objects in this method so that the loaded Observations
     * can be used in multiple threads by sharing a single ObservationFileLoader object.
     * @return
     */
    public List<Observation<Number>> getObservations() {
        List<Observation<Number>> rtn = new ArrayList<Observation<Number>>();
        for (Observation<Number> observation : sampleToObservation.values()) {
//            Observation<Number> copy = observation.copy(); // TODO: This needs to be checked if it is needed.
            rtn.add(observation);
        }
        return rtn;
    }
    
    /**
     * This method combined with setLoadedData() is used for multiple threading
     * environment to avoid multiple loading of a big data set multiple times 
     * by sharing the loaded typeToData.
     * @return
     */
    public Map<DataType, ObservationData> getLoadedData() {
        return this.typeToData;
    }
    
    /**
     * This method combined with getLoadedData() is used for multiple threading
     * environment to avoid multiple loading of a big data set multiple times 
     * by sharing the loaded typeToData.
     */
    public void setLoadedData(Map<DataType, ObservationData> typeToData) {
        this.typeToData = typeToData;
    }
    
    public Map<DataType, SharedEMFactors> getTypeToSharedFactors() {
        return this.typeToSharedFactors;
    }
    
    /**
     * Create accessory nodes and load observation data. The client should call this method if
     * the entities used in a factor graph are gene names. However, the factor graph converted
     * from Reactome pathway should not call this method since the entity names are PE's DB_IDs.
     * @param fileName
     * @param dataType
     * @param nameToVar
     * @throws Exception
     */
    public void processDataFile(String fileName,
                                DataType dataType,
                                Collection<Factor> factors,
                                VariableManager varManager) throws IOException {
        ObservationFactorHandler factorHandler = getObservationFactorHandler(dataType);
        if (factorHandler == null) {
            throw new IllegalStateException("ObservationFactorHandler is not set for " + dataType + "!");
        }
        if (dataType == DataType.Mutation) {
            loadMutationFile(fileName, factors, varManager);
            return;
        }
        parseNoMutationDataFile(fileName,
                                dataType,
                                factors,
                                varManager);
    }

    protected void parseNoMutationDataFile(String fileName, 
                                           DataType dataType,
                                           Collection<Factor> factors,
                                           VariableManager varManager) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] samples = line.split("\t");
        // Need to finish the loading by some post-processing steps.
        ObservationFactorHandler factorHandler = getObservationFactorHandler(dataType);
        List<Double> values = null;
        if (factorHandler instanceof EmpiricalFactorHandler)
            values = new ArrayList<Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Variable var = varManager.getVarForName(tokens[0]); // The first token should be gene name
            if (var == null)
                continue;
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].length() == 0 || tokens[i].equalsIgnoreCase("na"))
                    continue; // Unknown values
                ensureCentralDogmaNodes(tokens[0],
                                        factors,
                                        varManager);
                double value = new Double(tokens[i]);
                addObservation(value, 
                               tokens[0],
                               samples[i],
                               dataType,
                               varManager,
                               factors);
                if (values != null)
                    values.add(value);
            }
        }
        fu.close();
        if (factorHandler instanceof EmpiricalFactorHandler) {
            EmpiricalFactorHandler eFactorHandler = (EmpiricalFactorHandler) factorHandler;
            eFactorHandler.setDataForDistribution(dataType,
                                                  convertToDoubleArray(values));
        }
        factorHandler.finish();
    }
    
    /**
     * Mutation is not supported as default.
     * @param fileName
     * @param factors
     * @param nameToVar
     * @throws IOException
     */
    protected void loadMutationFile(String fileName,
                                    Collection<Factor> factors,
                                    VariableManager varManager) throws IOException {
        throw new IllegalArgumentException("Mutation data is not supported.");
    }
    
    /**
     * Load an observation file into a map as: Sample -> Gene -> Value
     * @param fileName
     * @param dataType
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, Float>> loadObservationData(String fileName,
                                                               DataType dataType) throws IOException {
        ObservationData data = typeToData.get(dataType);
        if (data != null)
            return data.sampleToGeneToValue;
        Map<String, Map<String, Float>> sampleToGeneToAssgn = _loadObservationData(fileName, dataType);
        data = new ObservationData();
        data.dataType = dataType;
        data.sampleToGeneToValue = sampleToGeneToAssgn;
        typeToData.put(data.dataType, data);
        return sampleToGeneToAssgn;
    }
    
    /**
     * Load an observation file into a map as: Sample -> Gene -> Assignment. The cached values
     * will not be used in this method if any.
     * @param fileName
     * @param dataType
     * @param threshold used to discretize continuous values.
     * @return
     * @throws IOException
     */
    private Map<String, Map<String, Float>> _loadObservationData(String fileName,
                                                                 DataType dataType) throws IOException {
        if (dataType == DataType.Mutation)
            return loadMAFFile(fileName);
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] samples = line.split("\t");
        Map<String, Map<String, Float>> sampleToGeneToValue = new HashMap<String, Map<String, Float>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].length() == 0 || tokens[i].equalsIgnoreCase("na"))
                    continue; // Unknown values
                Map<String, Float> geneToValue = sampleToGeneToValue.get(samples[i]);
                if (geneToValue == null) {
                    geneToValue = new HashMap<String, Float>();
                    sampleToGeneToValue.put(samples[i], geneToValue);
                }
                Float value = new Float(tokens[i]);
                geneToValue.put(gene, value);
            }
        }
        fu.close();
        return sampleToGeneToValue;
    }
    
    /**
     * In the default implementation, the mutation data type is not supported.
     * @param fileName
     * @return
     * @throws IOException
     */
    protected Map<String, Map<String, Float>> loadMAFFile(String fileName) throws IOException {
        throw new IllegalArgumentException("Mutation data is not supported.");
    }

    /**
     * Add observation data to a set of variables. Only genes in the passed targetGenes should be
     * handled.
     * @param sampleToGeneToValue
     * @param dataType
     * @param targetGenes
     * @param nameToVarInDogma
     * @param factors
     */
    public void addObservation(Map<String, Map<String, Float>> sampleToGeneToValue,
                               DataType dataType,
                               VariableManager varManager,
                               Collection<Factor> factors) {
        Set<String> targetGenes = grepTargetGenes(varManager.getNameToVar());
        for (String sample : sampleToGeneToValue.keySet()) {
            Map<String, Float> geneToValue = sampleToGeneToValue.get(sample);
            for (String gene : geneToValue.keySet()) {
                if (!targetGenes.contains(gene))
                    continue;
                addObservation(geneToValue.get(gene),
                               gene,
                               sample,
                               dataType,
                               varManager,
                               factors);
            }
        }
        // Need to finish the loading by some post-processing steps.
        ObservationFactorHandler factorHandler = getObservationFactorHandler(dataType);
        if (factorHandler instanceof EmpiricalFactorHandler) {
            EmpiricalFactorHandler eFactorHandler = (EmpiricalFactorHandler) factorHandler;
            double[] data = pulloutAllData(sampleToGeneToValue);
            eFactorHandler.setDataForDistribution(dataType, data);
        }
        factorHandler.finish();
    }
    
    private double[] pulloutAllData(Map<String, Map<String, Float>> sampleToGeneToValue) {
        List<Double> values = new ArrayList<Double>();
        for (String sample : sampleToGeneToValue.keySet()) {
            Map<String, Float> geneToValue = sampleToGeneToValue.get(sample);
            for (String gene : geneToValue.keySet()) {
                Float value = geneToValue.get(gene);
                values.add(value.doubleValue());
            }
        }
        double[] rtn = convertToDoubleArray(values);
        return rtn;
    }
    
    private double[] convertToDoubleArray(List<Double> values) {
        double[] rtn = new double[values.size()];
        for (int i = 0; i < values.size(); i++)
            rtn[i] = values.get(i);
        return rtn;
    }
    
    /**
     * Get a list of target genes that can be mapped to observation files.
     * @param nameToVarInDogma
     * @return
     */
    private Set<String> grepTargetGenes(Map<String, Variable> nameToVarInDogma) {
        Set<String> genes = new HashSet<String>();
        for (String name : nameToVarInDogma.keySet()) {
            if (name.endsWith("_" + PGMConfiguration.protein)) {
                int index = name.lastIndexOf("_");
                genes.add(name.substring(0, index));
            }
        }
        return genes;
    }

    protected void addObservation(double value, 
                                  String gene, 
                                  String sample,
                                  DataType dataType,
                                  VariableManager varManager,
                                  Collection<Factor> factors) {
        Variable anchorVar = dogmaHandler.getAnchorVar(gene,
                                                       varManager,
                                                       dataType);
        Variable valueVar = varManager.getVarForName(gene + "_" + dataType.toString());
        if (valueVar == null) {
            addObservationNode(anchorVar,
                               dataType,
                               gene,
                               factors,
                               varManager);
            valueVar = varManager.getVarForName(gene + "_" + dataType.toString());
        }
        addObservation(sample,
                       valueVar,
                       value,
                       dataType);
    }
    
    private void addObservation(String sample,
                                Variable var, 
                                Double value,
                                DataType dataType) {
        Observation<Number> observation = sampleToObservation.get(sample);
        if (observation == null) {
            observation = new Observation<Number>();
            observation.setName(sample);
            sampleToObservation.put(sample, observation);
        }
        ObservationFactorHandler factorHandler = getObservationFactorHandler(dataType);
        VariableAssignment<Number> varAssign = factorHandler.parseValue(value, dataType, var);
        observation.addAssignment(varAssign);
    }
    
    private void addObservationNode(Variable anchorVar,
                                    DataType type,
                                    String gene,
                                    Collection<Factor> factors,
                                    VariableManager varManager) {
        ObservationFactorHandler factorHandler = getObservationFactorHandler(type);
        Variable obsVar = factorHandler.createObservationVar(gene, 
                                                             type, 
                                                             varManager);
        Factor factor = factorHandler.createObservationFactor(anchorVar, 
                                                              obsVar,
                                                              type);
        factors.add(factor);
        if (factor instanceof EMFactor) {
            SharedEMFactors sharedFactors = getSharedEMFactors(type);
            sharedFactors.addSharedFactor((EMFactor)factor);
        }
    }

    private SharedEMFactors getSharedEMFactors(DataType type) {
        if (typeToSharedFactors == null) 
            typeToSharedFactors = new HashMap<DataType, SharedEMFactors>();
        SharedEMFactors sharedFactors = typeToSharedFactors.get(type);
        if (sharedFactors == null) {
            sharedFactors = new SharedEMFactors();
            typeToSharedFactors.put(type, sharedFactors);
        }
        return sharedFactors;
    }

    /**
     * Call this method to make sure the central dogma nodes for a specific sample has been added already.
     * @param gene
     * @param sample
     * @param nameToVar
     */
    protected void ensureCentralDogmaNodes(String gene,
                                           Collection<Factor> factors,
                                           VariableManager varManager) {
        if (dogmaHandler.isCentralDogmaAdded(gene, varManager))
            return;
        Variable var = varManager.getVarForName(gene);
        dogmaHandler.addCentralDogmaNodes(var, 
                                          gene, 
                                          factors,
                                          varManager);
    }
    
    /**
     * A very simple data structure to avoid a complicated map.
     * @author gwu
     *
     */
    public static class ObservationData {
        
        DataType dataType;
        // Use float to support both discrete and continuous variables.
        Map<String, Map<String, Float>> sampleToGeneToValue;
        
        public ObservationData() {
            
        }

        public DataType getDataType() {
            return dataType;
        }

        public void setDataType(DataType dataType) {
            this.dataType = dataType;
        }

        public Map<String, Map<String, Float>> getSampleToGeneToValue() {
            return sampleToGeneToValue;
        }

        public void setSampleToGeneToValue(Map<String, Map<String, Float>> sampleToGeneToValue) {
            this.sampleToGeneToValue = sampleToGeneToValue;
        }
    }
}
