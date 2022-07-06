/*
 * Created on Sep 6, 2018
 *
 */
package org.reactome.booleannetwork;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

/**
 * A simple identity function.
 * @author wug
 *
 */
@XmlAccessorType(XmlAccessType.FIELD)
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
