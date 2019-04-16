/*
 * Created on Apr 7, 2011
 *
 */
package org.reactome.r3.util;

/**
 * A simple class to store results from ProcessCaller.
 * @author wgm
 *
 */
public class ProcessCallResult {
    private String output;
    private String error;
    
    public ProcessCallResult() {
    }

    public String getOutput() {
        return output;
    }

    public void setOutput(String output) {
        this.output = output;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }
    
    
    
}
