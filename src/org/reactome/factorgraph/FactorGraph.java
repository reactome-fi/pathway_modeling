/*
 * Created on Feb 26, 2014
 *
 */
package org.reactome.factorgraph;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.jgrapht.GraphTests;
import org.jgrapht.graph.AsUndirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;

/**
 * A simple Factor graph for XML/JSON converting purpose.
 * @author gwu
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlRootElement
public class FactorGraph {
    // Add a name for debugging purpose
    private String name;
    // Use set for quick performance
    @XmlElement(name="factor")
    private Set<Factor> factors;
    // Keep all variables in the factor graph so that JAXB marshall/unmarshall
    // can create shared instances of variables. 
    @XmlElement(name="variable")
    private Set<Variable> variables;
    
    /**
     * Default constructor.
     */
    public FactorGraph() {
    }

    public Set<Factor> getFactors() {
        return factors;
    }

    public void setFactors(Set<Factor> factors) {
        this.factors = factors;
    }
    
    public void addFactor(Factor factor) {
        if (factors == null)
            factors = new HashSet<Factor>();
        factors.add(factor);
    }
    
    /**
     * Use this method to make sure cached variables are correct.
     * This method should be called as the last step when a new
     * PGMFactorGraph is created.
     */
    public void validatVariables() {
        if (variables == null)
            variables = new HashSet<Variable>();
        else
            variables.clear();
        if (factors == null) {
            return;
        }
        for (Factor factor : factors) {
            if (factor.getVariables() != null) {// Just in case
                variables.addAll(factor.getVariables());
                for (Variable var : factor.getVariables())
                    var.addFactor(factor);
            }
        }
    }
    
    /**
     * Set FGNode ids so that they are unique. Variable and Factor objects should have
     * different ids. If an Variable and Factor have the same id, they may be messed up
     * during unmarshalling.
     * @param fg
     */
    public void setIdsInFactors() {
        int id = 0;
        if (getVariables() != null) {
            for (Variable var : getVariables()) {
                var.setId(id ++);
            }
        }
        if (getFactors() != null) {
            for (Factor factor : getFactors())
                factor.setId(id ++);
        }
    }
    
    /**
     * Calculate the log likelihood of an assignment of all variables contained by this 
     * FactorGraph object. The value actually is the log of products of all factors determined
     * by the passed assignment.
     * @param varToState
     * @return
     */
    public double getLogLikelihood(Map<Variable, Integer> varToState) {
        // Check to make sure all variables have been assigned
        Set<Variable> shared = new HashSet<Variable>(getVariables());
        shared.removeAll(varToState.keySet());
        if (shared.size() > 0) {
            throw new IllegalArgumentException("States of one or more variables contained by this FactorGraph are not specified.");
        }
        // Use log-space to calculation
        double loglikelihood = 0.0d;
        for (Factor factor : factors) {
            double value = factor.getValue(varToState);
            loglikelihood += Math.log(value);
        }
        return loglikelihood;
    }
    
    public void setVariables(Set<Variable> variables) {
        this.variables = variables;
    }
    
    public Set<Variable> getVariables() {
        return variables;
    }
    
    public Variable getVariable(String varName) {
        if (variables == null)
            return null;
        for (Variable var : variables) {
            if (var.getName().equals(varName))
                return var;
        }
        return null;
    }
    
    /**
     * Export this object using JAXB.
     * @param os
     * @throws JAXBException
     */
    public void exportFG(OutputStream os) throws JAXBException {
        JAXBContext context = JAXBContext.newInstance(FactorGraph.class);
        Marshaller marshaller = context.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        marshaller.marshal(this, os);
    }
    
    /**
     * Import a serialized FactorGraph from an inputstream. 
     */
    public void importFG(InputStream is) throws JAXBException {
        JAXBContext context = JAXBContext.newInstance(FactorGraph.class);
        Unmarshaller unmarshaller = context.createUnmarshaller();
        FactorGraph fg = (FactorGraph) unmarshaller.unmarshal(is);
        this.factors = fg.getFactors();
        this.variables = fg.getVariables();
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    /**
     * Check if this FactorGraph is a tree, which doesn't have any loop. This method uses
     * the data structures in jGraphT.
     * @return
     */
    public boolean isTree() {
        if (factors == null || factors.size() == 0)
            return true; // Just an empty graph.
        // Convert it as a JGraphT to check
        DefaultDirectedGraph<FGNode, DefaultEdge> graph = new DefaultDirectedGraph<FGNode, DefaultEdge>(DefaultEdge.class);
        for (Factor factor : factors) {
            graph.addVertex(factor);
            List<Variable> variables = factor.getVariables();
            for (Variable var : variables) {
                if (!graph.containsVertex(var))
                    graph.addVertex(var);
                graph.addEdge(factor, var);
            }
        }
        // There is no need to consider directions
        AsUndirectedGraph<FGNode, DefaultEdge> unDirGraph = new AsUndirectedGraph<FGNode, DefaultEdge>(graph);
        return GraphTests.isTree(unDirGraph);
    }
    
    /**
     * Check if this FactorGraph can be used for inference. Currently the Inference algorithms implemented
     * can be used for discrete variables only FactorGraphs or FactorGraphs containing leaf ContinuousVariable and
     * ContinuousFactor with discrete variables.
     * @return
     */
    public boolean isInferreable() {
        List<ContinuousVariable> continuousVars = new ArrayList<ContinuousVariable>();
        for (Variable var : getVariables()) {
            if (var instanceof ContinuousVariable)
                continuousVars.add((ContinuousVariable)var);
        }
        if (continuousVars.size() == 0)
            return true;
        // Make sure ContinousVariable is contained by a single ContinuousFactor only
        for (ContinuousVariable var : continuousVars) {
            List<Factor> factors = var.getFactors();
            if (factors.size() > 1)
                return false;
            if (factors.size() == 0)
                continue;
            Factor factor = factors.get(0);
            if (!(factor instanceof ContinuousFactor)) {
                return false;
            }
            // Make sure there are only two variables: one this var and another is a discrete Variable
            List<Variable> vars = factor.getVariables();
            if (vars.size() > 2)
                return false;
            for (Variable var1 : vars) {
                if (var1 == var)
                    continue;
                if (var1 instanceof ContinuousVariable)
                    return false; // There are two CLGVariable here
            }
        }
        return true;
    }
    
    @Override
    public String toString() {
        if (name != null)
            return name;
        return super.toString();
    }
    
}
