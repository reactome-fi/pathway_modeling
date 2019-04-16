/*
 * Created on Apr 17, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A simple class to model an attractor, which may be a single value or a cycle.
 * @author gwu
 *
 */
public class Attractor {
    private List<BooleanVariable> variables;
    // This may be a single element or multiple elements
    private List<Number[]> values; 
    
    /**
     * Default constructor.
     */
    public Attractor() {
    }

    public List<BooleanVariable> getVariables() {
        return variables;
    }

    public void setVariables(List<BooleanVariable> variables) {
        this.variables = variables;
    }

    public List<Number[]> getValues() {
        return values;
    }

    public void setValues(List<Number[]> values) {
        this.values = values;
    }
    
    public Map<BooleanVariable, List<Number>> getVarToValues() {
        Map<BooleanVariable, List<Number>> varToValues = new HashMap<>();
        for (int i = 0; i < variables.size(); i++) {
            BooleanVariable var = variables.get(i);
            List<Number> varValues = new ArrayList<>();
            for (int j = 0; j < values.size(); j++) {
                varValues.add(values.get(j)[i]);
            }
            varToValues.put(var, varValues);
        }
        return varToValues;
    }
    
    public String outputAsText() {
        if (variables == null || values == null)
            return "No information available!";
        StringBuilder builder = new StringBuilder();
        builder.append("Variable");
        for (int i = 0; i < values.size(); i++) {
            builder.append("\t").append(i + 1);
        }
        builder.append("\n");
        for (int i = 0; i < variables.size(); i++) {
            builder.append(variables.get(i).getName());
            for (int j = 0; j < values.size(); j++) {
                builder.append("\t").append(values.get(j)[i]);
            }
            builder.append("\n");
        }
        return builder.toString();
    }
    
}
