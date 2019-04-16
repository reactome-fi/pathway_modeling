/*
 * Created on Mar 30, 2015
 *
 */
package org.reactome.factorgraph;


/**
 * This class is used to model continuous variable. In the implementation of this PGM package, such kind of variable 
 * will be used only as a leaf node in a FactorGraph to support discrete inference, and shouldn't be included in 
 * inference. Since they are leaf nodes, they will be integrated into their related factors.
 * @author gwu
 *
 */
public class ContinuousVariable extends Variable {
    // Specify the value distribution of this Continuous should be 
    // regarded as two-sided (e.g. gene expression) or one-sided
    // (e.g. mutation). The default will be ONE_SIDED.
    private DistributionType distributionType = DistributionType.ONE_SIDED;
    
    /**
     * Default constructor.
     */
    public ContinuousVariable() {
    }
    
    /**
     * @param states
     */
    public ContinuousVariable(int states) {
        super(states);
    }
    
    public DistributionType getDistributionType() {
        return distributionType;
    }

    public void setDistributionType(DistributionType distributionType) {
        this.distributionType = distributionType;
    }

    @Override
    protected double[] sendMessage(FGNode target, 
                                   InferenceType type,
                                   boolean logSpace) {
        throw new IllegalStateException("This method is not supported in class ContinuousVariable");
    }

    @Override
    protected void updateBelief(boolean logSpace) {
        throw new IllegalStateException("This method is not supported in class ContinuousVariable");
    }

    /**
     * For the time being, this type of query is not supported.
     */
    @Override
    public double[] getBelief() {
        throw new IllegalStateException("This method is not supported in class ContinuousVariable");
    }
    
    /**
     * The distribution type of the values modeled by a ContinuousVariable. For example,
     * a gene expression value should be TWO_SIDED, and a mutation should be ONE_SIDED.
     * @author gwu
     *
     */
    public static enum DistributionType {
        ONE_SIDED,
        TWO_SIDED
    }
    
}
