/*
 * Created on Mar 7, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlID;
import javax.xml.bind.annotation.XmlTransient;

/**
 * An abstract class for both factor and variables.
 * @author gwu
 */
@XmlAccessorType(XmlAccessType.FIELD)
public abstract class FGNode {
    // For some customized information
    // Don't annotate this map as element. Otherwise, nothing will be
    // output in JAXB
    private HashMap<String, String> prop;
    // A unique id in the containing factor graph.
    @XmlID
    private String id;
    // A human readable name. Two PGMNode may have the same name.
    @XmlAttribute
    private String name;
    // edges for messages sent in
    @XmlTransient // We don't want to serialize these values
    protected List<Edge> inEdges;
    // edges for messages sent out
    @XmlTransient
    protected List<Edge> outEdges;
    @XmlTransient
    protected double[] belief;
    // Used to hold messages to avoid GC 
    @XmlTransient
    protected double[] message;
    
    public FGNode() {
    }
    
    public String getProperty(String key) {
        if (prop == null)
            return null;
        return prop.get(key);
    }

    public void setProperty(String key, String value) {
        // We want to keep the use of memory at minimum
        if (key == null) // Do nothing
            return;
        if (value == null) {
            if (prop != null)
                prop.remove(key); // Don't store a null value
        }
        else {
            if (prop == null)
                prop = new HashMap<String, String>();
            prop.put(key, value);
        }
    }
    
    /**
     * Set some customized information here for future reference.
     * @param map
     */
    public void setProperties(Map<String, String> map) {
        if (map == null)
            prop = null;
        else
            prop = new HashMap<String, String>(map);
    }
    
    public Map<String, String> getProperties() {
        return this.prop;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }
    
    public void setId(Integer id) {
        this.id = id + "";
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    public String toString() {
        return getName() + " (" + getId() + ")";
    }
    
    public void resetEdges() {
        if (inEdges != null)
            inEdges.clear();
        if (outEdges != null)
            outEdges.clear();
    }

    public List<Edge> getInEdges() {
        return inEdges;
    }

    public void setInEdges(List<Edge> inEdges) {
        this.inEdges = inEdges;
    }
    
    public void addInEdge(Edge edge) {
        if (inEdges == null)
            inEdges = new ArrayList<Edge>();
        inEdges.add(edge);
    }
    
    public void removeInEdge(Edge edge) {
        if (inEdges != null)
            inEdges.remove(edge);
    }

    public List<Edge> getOutEdges() {
        return outEdges;
    }

    public void setOutEdges(List<Edge> outEdges) {
        this.outEdges = outEdges;
    }
    
    public void addOutEdge(Edge edge) {
        if (outEdges == null)
            outEdges = new ArrayList<Edge>();
        outEdges.add(edge);
    }
    
    /**
     * Remove a registered outward edge pointing to the targeted
     * target.
     * @param target
     */
    public void removeOutEdge(FGNode target) {
        if (outEdges == null || outEdges.size() == 0)
            return;
        for (Iterator<Edge> it = outEdges.iterator(); it.hasNext();) {
            Edge edge = it.next();
            if (edge.getToNode() == target) {
                it.remove();
                break;
            }
        }
    }
    
    public void removeOutEdge(Edge edge) {
        if (outEdges != null)
            outEdges.remove(edge);
    }
    
    /**
     * Send message from this FGNode to a target FGNode. 
     * @param target
     * @param inferenceType SUM_PROD or MAX_PROD
     * @param logSpace true to use log-space for multiplication.
     * @return
     */
    protected abstract double[] sendMessage(FGNode target,
                                            InferenceType inferenceType,
                                            boolean logSpace);

    protected void normalize(double[] belief,
                             boolean logSpace,
                             boolean needBackToLogSpace) {
        if (logSpace)
            convertLogToProb(belief);
        // Do a normalization
        double sum = 0.0d;
        for (int i = 0; i < belief.length; i++)
            sum += belief[i];
        for (int i = 0; i < belief.length; i++)
            belief[i] /= sum;
        if (needBackToLogSpace)
            convertProbToLog(belief);
    }
    
    protected String convertToString(double[] belief) {
        StringBuilder builder = new StringBuilder();
        for (double v : belief)
            builder.append(v).append("\t");
        return builder.toString();
    }

    public double[] getBelief() {
        return belief;
    }

    public void setBelief(double[] belief) {
        this.belief = belief;
    }
    
    /**
     * This method should be called after an inference run. Otherwise, the
     * results will not be correct.
     */
    protected abstract void updateBelief(boolean logSpace);

    protected void convertLogToProb(double[] message) {
        // Need to convert back from log space to probability space.
        // Use a trick described on page 360 in the PGM book in case 
        // probabilities are too small causing underflow
        // Find the largest values
        double max = message[0]; // Since we actually don't know what is the minimum value
        for (int i = 1; i < message.length; i++) {
            if (message[i] > max)
                max = message[i];
        }
        // Convert from log value to probability value
        for (int i = 0; i < message.length; i++) {
            message[i] = Math.exp(message[i] - max);
        }
    }
    
    protected void convertProbToLog(double[] message) {
        for (int i = 0; i < message.length; i++)
            message[i] = Math.log(message[i]);
    }
    
}
