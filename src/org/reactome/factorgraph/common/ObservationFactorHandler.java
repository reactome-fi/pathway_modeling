/*
 * Created on Jul 7, 2015
 *
 */
package org.reactome.factorgraph.common;

import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;

/**
 * This class is used to handle the observed values. This is an Abstract class, an subclass to
 * this class should provide implement to handle values.
 * @author gwu
 *
 */
public abstract class ObservationFactorHandler {
    protected PGMConfiguration configuration;
    
    public PGMConfiguration getConfiguration() {
        return configuration;
    }

    public void setConfiguration(PGMConfiguration configuration) {
        this.configuration = configuration;
    }

    /**
     * Parse a value token loaded from a file.
     * @param dataType
     * @param valueToken
     * @param configuration
     */
    public abstract VariableAssignment<Number> parseValue(Double value,
                                                          DataType dataType,
                                                          Variable obsVar);
    
    /**
     * Create a Variable for an observation data item. The created Variable may be a discrete
     * or continuous Variable object.
     * @param name
     * @param dataType
     * @return
     */
    public abstract Variable createObservationVar(String geneName, 
                                         DataType dataType,
                                         VariableManager variableManager);
    
    /**
     * Create a Factor object for an observation data item.
     * @param anchorVar
     * @param obsVar
     * @param dataType
     */
    public abstract Factor createObservationFactor(Variable anchorVar,
                                                   Variable obsVar,
                                                   DataType dataType);
    
    /**
     * Perform some post-processing work. The default is doing nothing.
     */
    public void finish() {}
}
