/*
 * Created on Sep 6, 2018
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

/**
 * An implementation of transfer function based on constrained Hill function.
 * @author wug
 *
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class HillFunction implements TransferFunction {
    // Currently hard-coded for the time being
    private int n = 3;
    private double g = 1.0d;
    private double k = 0.5503d; // Used as EC50. This may need to be changed.
    
    public HillFunction() {
    }
    
    public void setParameter_g(double g) {
    	this.g = g;
    }
    
    /**
     * Configure parameters for Hill functions. See the original paper for the
     * meanings of these parameters (Figure 1): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3292705/.
     */
    public void setParameters(int n,
                              double k,
                              double g) {
        this.n = n;
        this.k = k;
        this.g = g;
    }
    
    public Map<String, Number> getParameters() {
        Map<String, Number> paraToValue = new HashMap<>();
        paraToValue.put("n", n);
        paraToValue.put("k", k);
        paraToValue.put("g", g);
        return paraToValue;
    }

    @Override
    public Double transfer(Number input) {
        if (input == null)
            return null;
        double kpowern = Math.pow(k, n);
        double vpowern = Math.pow(input.doubleValue(), n);
        return g * (1 + kpowern) * vpowern / (vpowern + kpowern);
    }
    
    @Override
    public String toString() {
        return "Hill function (n: " + n + ", k: " + k + ", g: " + g + ")";
    }

}
