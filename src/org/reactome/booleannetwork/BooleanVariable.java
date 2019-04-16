/*
 * Variable is used as a node in a BooleanNetwork object.
 * Created on Mar 23, 2017
 */
package org.reactome.booleannetwork;

import java.util.HashSet;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlIDREF;
import javax.xml.bind.annotation.XmlTransient;

@XmlAccessorType(XmlAccessType.FIELD)
public class BooleanVariable extends BooleanObject {
	
    // Kept values from simulation for post analysis
    @XmlTransient
    private Number[] track;
	private Number value; // Since we want to use fuzzy logic, the value may not be just 0 or 1.
	@XmlElement(name="inRelation")
	@XmlIDREF
	private Set<BooleanRelation> inRelations;
	@XmlElement(name="outRelation")
	@XmlIDREF
	private Set<BooleanRelation> outRelations;
	// If true, an identical transfer function will be used (i.e. input copied to output)
	private boolean useIdenticalTransfer; 

	public BooleanVariable() {
	}
	
	public Number[] getTrack() {
        return track;
    }

	/**
	 * Values generated from simulation. The first number should be the initial
	 * value (ie value at time = 0).
	 * @param track
	 */
    public void setTrack(Number[] track) {
        this.track = track;
    }

    public boolean isUseIdenticalTransfer() {
        return useIdenticalTransfer;
    }
    
    /**
     * If a BooleanVariable doesn't have anything connected into it, it is regarded as
     * an input variable.
     * @return
     */
    public boolean isInputVariable() {
        if (inRelations == null || inRelations.size() == 0)
            return true;
        return false;
    }

    public void setUseIdenticalTransfer(boolean useIdenticalTransfer) {
        this.useIdenticalTransfer = useIdenticalTransfer;
    }

    public Number getValue() {
		return value;
	}

	public void setValue(Number value) {
		this.value = value;
	}

	public Set<BooleanRelation> getInRelations() {
		return inRelations;
	}

	public void addInRelation(BooleanRelation inRelation) {
		if (inRelations == null)
			inRelations = new HashSet<BooleanRelation>();
		inRelations.add(inRelation);
	}

	public Set<BooleanRelation> getOutRelations() {
		return outRelations;
	}
	
	public void addOutRelation(BooleanRelation outRelation) {
		if (outRelations == null)
			outRelations = new HashSet<BooleanRelation>();
		outRelations.add(outRelation);
	}

	public String toString() {
		return "Variable: " + name;
	}
	
}
