/*
 * Created on Sep 6, 2018
 *
 */
package org.reactome.booleannetwork;

/**
 * A simple identity function.
 * @author wug
 *
 */
public class IdentityFunction implements TransferFunction {
    
    public IdentityFunction() {
        
    }

    @Override
    public Double transfer(Number value) {
        if (value == null)
            return null;
        return value.doubleValue();
    }
    
    @Override
    public String toString() {
        return "Identity function";
    }

}
