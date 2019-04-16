/*
 * Created on Jun 5, 2014
 *
 */
package org.reactome.factorgraph;


/**
 * An edge between two FGNode objects. This edge has a direction.
 * @author gwu
 *
 */
public class Edge {
    private FGNode fromNode;
    private FGNode toNode;
    // Message along this edge
    // Use a double array should be faster than using a double ArrayList
    private double[] message;
    
    /**
     * Default constructor.
     */
    public Edge() {
    }
    
    public Edge(FGNode fromNode, FGNode toNode) {
        this.fromNode = fromNode;
        this.toNode = toNode;
    }

    public FGNode getFromNode() {
        return fromNode;
    }

    public void setFromNode(FGNode fromNode) {
        this.fromNode = fromNode;
    }

    public FGNode getToNode() {
        return toNode;
    }

    public void setToNode(FGNode toNode) {
        this.toNode = toNode;
    }

    public double[] getMessage() {
        return message;
    }

    public void setMessage(double[] message) {
        // Want to make a copy to avoid reuse
        if (this.message == null)
            this.message = new double[message.length];
        System.arraycopy(message, 0, this.message, 0, this.message.length);
    }
    
    /**
     * Initialize message to 1.0 based on variable.
     */
    protected void initializeMessage(double initialMessage, 
                                     boolean logSpace) {
        // Need to find the state
        int state = 0;
        if (fromNode instanceof Variable)
            state = ((Variable)fromNode).getStates();
        else
            state = ((Variable)toNode).getStates();
        message = new double[state];
        for (int i = 0; i < state; i++)
            message[i] = (logSpace ? Math.log(initialMessage) : initialMessage);
    }
}
