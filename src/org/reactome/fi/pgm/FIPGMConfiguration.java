/*
 * Created on May 28, 2014
 *
 */
package org.reactome.fi.pgm;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.PGMConfiguration;
import org.reactome.r3.util.FileUtility;


/**
 * Customized configuration file for the PGM converted from the FI network.
 * @author gwu
 *
 */
public class FIPGMConfiguration extends PGMConfiguration {
    private static final Logger logger = Logger.getLogger(FIPGMConfiguration.class);
    
    private static final int NUMBER_OF_STATE = 2;
    // Single gene in the normal state vs. the abnormal state
    public static final double ALPHA = 0.01d;
    public static final double THETA = 0.1d;
    // Background correlation.
    public static final double BETA = 0.01d;
    // Define a background noise
    public static final double GAMMA = 0.01d;
    public static final String RESULT_DIR = "results/FI_PGM/";
    public static final String DATA_DIR = "test_data/";
    
    // This is a singleton.
    private static FIPGMConfiguration config;
    // Check if miRNA-target interactions are needed
    private boolean needmiRNA = false; // Default should not use miRNA
    // Cached loaded FIs
    private Set<String> loadedFIs;
    
    private FIPGMConfiguration() {
        setNumberOfStates(NUMBER_OF_STATE);
        double[] values = new double[4];
        values[0] = 1.0d - FIPGMConfiguration.BETA;
        values[1] = FIPGMConfiguration.BETA;
        values[2] = FIPGMConfiguration.BETA;
        values[3] = 1.0d - FIPGMConfiguration.BETA;
        setCentralDogmaValues(values);
        File file = new File("resources/PGM_FI_Config.xml");
        if (file.exists()) {
            try {
                configure("resources/PGM_FI_Config.xml");
                debug();
            }
            catch(Exception e) {
                logger.error("FIPGMConfiguration: " + e.getMessage(), e);
            }
        }
    }
    
    public void setFIs(Set<String> fis) {
        this.loadedFIs = fis;
    }
    
    public boolean isNeedmiRNA() {
        return needmiRNA;
    }

    public void setNeedmiRNA(boolean needmiRNA) {
        this.needmiRNA = needmiRNA;
    }

    public static final FIPGMConfiguration getConfig() {
        if (config == null)
            config = new FIPGMConfiguration();
        return config;
    }
    
    public Set<String> getFIs() throws IOException {
        if (loadedFIs == null)
            loadedFIs = loadFIs();
        return loadedFIs;
    }
    
    private Set<String> loadFIs() throws IOException {
        FileUtility fu = new FileUtility();
        String fiFileName = getProperties().get("fiFile");
        if (fiFileName == null) {
            logger.error("fiFile is not set!");
            throw new IllegalStateException("fiFile is not set!");
        }
        Set<String> fis = fu.loadInteractions(fiFileName);
        if (!needmiRNA)
            return fis;
        String miRNATargetFile = getProperties().get("miRNATargetFile");
        if (miRNATargetFile == null) {
            logger.warn("miRNATargetFile is not set though miRNATarget interactions are required!");
            return fis;
        }
        fu.setInput(miRNATargetFile);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fis.add(tokens[0] + "\t" + tokens[1]);
        }
        fu.close();
        return fis;
    }
    
    /**
     * Print out loaded information for debugging purpose.
     */
    private void debug() {
        if(logger.getLevel() == null || logger.getLevel().isGreaterOrEqual(Level.WARN))
            return;
        logger.info("Loaded configurations:");
        // Output for INFO or Lower
        if (properties != null) {
            logger.info("Properties:");
            for (String key : properties.keySet()) {
                logger.info(key + ": " + properties.get(key));
            }
        }
        StringBuilder builder = new StringBuilder();
        if (typeToFactorValues != null) {
            logger.info("typeToFactorValues:");
            for (DataType type : typeToFactorValues.keySet()) {
                builder.setLength(0);
                double[] values = typeToFactorValues.get(type);
                for (double value : values)
                    builder.append(value).append(",");
                logger.info(type + ": " + builder.toString());
            }
        }
        if (typeToFileName != null) {
            logger.info("typeToFileName: ");
            for (DataType type : typeToFileName.keySet()) {
                logger.info(type + ": " + typeToFileName.get(type));
            }
        }
        if (typeToThreshold != null) {
            logger.info("typeToThreshold:");
            for (DataType type : typeToThreshold.keySet()) {
                logger.info(type + ": " + typeToThreshold.get(type));
            }
        }
    }
}
