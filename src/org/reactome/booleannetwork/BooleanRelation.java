/*
 * A BooleanRelation models a Boolean relation between two or more BooleanVariable objects.
 * If more than two BooleanVariables are included, make sure there is only one output and all
 * others are inputs that have "AND" relations.
 * Created on Mar 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlIDREF;
import javax.xml.bind.annotation.XmlTransient;

@XmlAccessorType(XmlAccessType.FIELD)
public class BooleanRelation extends BooleanObject {
	
    @XmlElement(name="inputToIsNegated")
    private List<InputToIsNegated> inputToIsNegateds;
	@XmlElement(name="output")
	@XmlIDREF
	private BooleanVariable outputVariable;
	// A relation can have its own transferFunction. This is optional.
	// For the time being, this is not supported at the server-side for ReactomeFIViz
	@XmlTransient
	private TransferFunction transferFunction;
	
	public BooleanRelation() {
	}

	@XmlTransient
	public TransferFunction getTransferFunction() {
		return transferFunction;
	}

	public void setTransferFunction(TransferFunction transferFunction) {
		this.transferFunction = transferFunction;
	}

	/**
	 * Don't try to modify the returned map.
	 * @return
	 */
	public Map<BooleanVariable, Boolean> getInputVarToIsNegated() {
	    Map<BooleanVariable, Boolean> rtn = new HashMap<>();
	    if (inputToIsNegateds != null) {
	        for (InputToIsNegated tmp : inputToIsNegateds)
	            rtn.put(tmp.var, tmp.isNegated);
	    }
	    return rtn;
	}
	
	public Set<BooleanVariable> getInputVariables() {
		if (inputToIsNegateds == null)
			return null;
		Set<BooleanVariable> vars = new HashSet<>();
		for (InputToIsNegated tmp : inputToIsNegateds) 
		    vars.add(tmp.var);
		return vars;
	}

	public void addInputVariable(BooleanVariable variable,
	                             Boolean isNegated) {
	    if (variable == null)
	        return;
	    if (inputToIsNegateds == null)
	        inputToIsNegateds = new ArrayList<>();
	    InputToIsNegated inputToIsNegated = new InputToIsNegated();
	    inputToIsNegated.var = variable;
	    inputToIsNegated.isNegated = isNegated;
	    inputToIsNegateds.add(inputToIsNegated);
	    variable.addOutRelation(this);
	}

	public BooleanVariable getOutputVariable() {
		return outputVariable;
	}

	public void setOutputVariable(BooleanVariable outputVariable) {
	    if (outputVariable == null)
	        return;
	    this.outputVariable = outputVariable;
	    outputVariable.addInRelation(this);
	}
	
	public String toString() {
		return inputToIsNegateds + " -> " + outputVariable;
	}
	
	/**
	 * A helper class for JAXB
	 */
	@XmlAccessorType(XmlAccessType.FIELD)
	private static class InputToIsNegated {
	    @XmlElement(name="input")
	    @XmlIDREF
	    private BooleanVariable var;
	    private Boolean isNegated;
	    
	    public InputToIsNegated() {
	    }
	    
	    public String toString() {
	        return var.getName() + ", " + isNegated;
	    }
	    
	}

}
