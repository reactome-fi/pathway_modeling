/*
 * Created on Sep 7, 2017
 *
 */
package org.reactome.pathway.booleannetwork;

public class DefaultAffinityToModificationMap implements AffinityToModificationMap {
    
    public DefaultAffinityToModificationMap() {
    }
    
    @Override
    public double getModificationStrenth(double value) {
        if (value < 1) // Less than 1 nM
            return 0.999d;
        if (value < 10)
            return 0.99d;
        if (value <= 10e6) {
            // Linear relation between log[10, 10e6] ~ [0.99, 0.01]
            double k = (0.01d - 0.99d) / (6.0d - 1.0d);
            double rtn = k * (Math.log10(value) - 1) + 0.99d;
            return rtn;
        }
        return 0.0d;
    }

}
