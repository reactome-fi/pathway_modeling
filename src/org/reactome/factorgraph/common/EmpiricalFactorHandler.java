/*
 * Created on Jul 9, 2015
 *
 */
package org.reactome.factorgraph.common;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.reactome.factorgraph.ContinuousVariable;
import org.reactome.factorgraph.EmpiricalFactor;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;
import org.reactome.factorgraph.VariableAssignment;

/**
 * This class is used to handle factors and variables related to empirical distribution.
 * @author gwu
 *
 */
public class EmpiricalFactorHandler extends ContinuousObservationFactorHandler {
    // Cached distribution for each data type
    protected Map<DataType, EmpiricalDistribution> typeToDist;
    
    /**
     * Default constructor.
     */
    public EmpiricalFactorHandler() {
    }
    
    /**
     * Get the cached EmpiricalDistribution for a DataType.
     * @param dataType
     * @return
     */
    public EmpiricalDistribution getDistribution(DataType dataType) {
        if (typeToDist == null)
            typeToDist = new HashMap<DataType, EmpiricalDistribution>();
        EmpiricalDistribution dist = typeToDist.get(dataType);
        if (dist == null) {
            dist = new EmpiricalDistribution();
            typeToDist.put(dataType, dist);
        }
        return dist;
    }
    
    /**
     * Set the data for the distribution specified by its DataType object.
     * @param dataType
     * @param data
     */
    public void setDataForDistribution(DataType dataType,
                                       double[] data) {
        EmpiricalDistribution dist = getDistribution(dataType);
        dist.load(data);
    }
    
    @Override
    public VariableAssignment<Number> parseValue(Double value,
                                                 DataType dataType,
                                                 Variable obsVar) {
        VariableAssignment<Number> varAssgn = super.parseValue(value, dataType, obsVar);
        // Need to inject an EmpricalDistribution object for it.
        EmpiricalDistribution dist = getDistribution(dataType);
        varAssgn.setDistribution(dist);
        return varAssgn;
    }

    @Override
    public Variable createObservationVar(String geneName,
                                         DataType dataType,
                                         VariableManager variableManager) {
        Variable var = super.createObservationVar(geneName, dataType, variableManager);
        ((ContinuousVariable)var).setDistributionType(dataType.getDistributionType());
        return var;
    }

    @Override
    public Factor createObservationFactor(Variable anchorVar,
                                          Variable obsVar,
                                          DataType dataType) {
        if (!(obsVar instanceof ContinuousVariable))
            throw new IllegalArgumentException("The passed obsVar should be a ContinousVariable.");
        ContinuousVariable cVar = (ContinuousVariable) obsVar;
        EmpiricalFactor factor = new EmpiricalFactor();
        factor.setDiscreteVariable(anchorVar);
        factor.setContinuousVariable(cVar);
        return factor;
    }

    /**
     * Provide a check to make sure the data in distributions are not empty.
     */
    @Override
    public void finish() {
        if (typeToDist == null || typeToDist.size() == 0)
            return;
        for (DataType dataType : typeToDist.keySet()) {
            EmpiricalDistribution dist = typeToDist.get(dataType);
            if (!dist.isLoaded())
                throw new IllegalStateException("The EmpricalDistribution object for " + dataType + " is not loaded.");
        }
    }
    
}
