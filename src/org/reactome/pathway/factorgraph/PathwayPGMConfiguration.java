/*
 * Created on Oct 31, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.File;

import org.apache.log4j.Logger;
import org.reactome.factorgraph.common.PGMConfiguration;

/**
 * Customized PGMConfiguration to provide some specific information for PGMs converted
 * from Reactome pathways. This is a singleton so that it can be shared in different places.
 * @author gwu
 *
 */
public class PathwayPGMConfiguration extends PGMConfiguration {
    private final Logger logger = Logger.getLogger(PathwayPGMConfiguration.class);
    private static PathwayPGMConfiguration config;
    
    private PathwayPGMConfiguration() {
        setNumberOfStates(PathwayFGConstants.NUMBER_OF_STATES);
        // This is hard-coded
        String configFileName = "resources/PGM_Pathway_Config.xml";
        File file = new File(configFileName);
        if (!file.exists())
            return;
        try {
            configure(configFileName);
        }
        catch(Exception e) {
            logger.error("PathwayPGMConfiguration: " + e, e);
        }
    }
    
    public static PathwayPGMConfiguration getConfig() {
        if (config == null)
            config = new PathwayPGMConfiguration();
        return config;
    }
}
