/*
 * Created on May 28, 2014
 *
 */
package org.reactome.factorgraph.common;

import org.reactome.factorgraph.ContinuousVariable.DistributionType;


public enum DataType { 
    CNV,
    mRNA_EXP,
    Methylation,
    miRNA,
    Mutation;

    /**
     * Call this method so that a unique key can be generated for a set of data types. The
     * returned stride is used to generate a key for a HashMap. For example, a data containing
     * CNV = 1, mRNA_Exp = 0, Methylation = 1, miRNA = null, and Mutation = 1 should have the 
     * following key: 1 * 2 (for CNV) + 10 * 1 (mRNA_Exp) + 100 * 2 (Methylation) + 1000 * 0 (miRNA) + 
     * 10000 * 2 (for Mutation) = 20212
     * @param type
     * @return
     */
    public static int getKeyStride(DataType type) {
        switch (type) {
            case CNV : return 1;
            case mRNA_EXP : return 10;
            case Methylation : return 100;
            case miRNA : return 1000;
            case Mutation : return 10000;
        }
        throw new IllegalArgumentException("Unknown DataType: " + type);
    }
    
    public DistributionType getDistributionType() {
        if (this == Mutation)
            return DistributionType.ONE_SIDED;
        else
            return DistributionType.TWO_SIDED;
    }
    
}