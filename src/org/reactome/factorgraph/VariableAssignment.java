/*
 * Created on Jul 2, 2015
 *
 */
package org.reactome.factorgraph;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlIDREF;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.commons.math3.random.EmpiricalDistribution;

/**
 * This class is used to describe an assignment to a Variable object.
 * @author gwu
 *
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class VariableAssignment<T extends Number> {
    // The Variable whose value is assigned
    @XmlIDREF
    private Variable variable;
    // The assigned value: this may be a state in an Integer or a continuous value
    private T assignment;
    // A description about this value
    @XmlTransient
    private EmpiricalDistribution distribution;
    
    /**
     * Default constructor.
     */
    public VariableAssignment() {
    }

    public Variable getVariable() {
        return variable;
    }

    public void setVariable(Variable variable) {
        this.variable = variable;
    }

    public T getAssignment() {
        return assignment;
    }

    public void setAssignment(T assignment) {
        this.assignment = assignment;
    }

    public EmpiricalDistribution getDistribution() {
        return distribution;
    }

    public void setDistribution(EmpiricalDistribution distribution) {
        this.distribution = distribution;
    }
    
}
