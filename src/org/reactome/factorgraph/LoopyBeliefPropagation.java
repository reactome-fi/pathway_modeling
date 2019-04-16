/*
 * Created on Jun 5, 2014
 *
 */
package org.reactome.factorgraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * This class implements the algorithm of loopy belief propagation.
 * @author gwu
 *
 */
public class LoopyBeliefPropagation extends AbstractInferencer {
    private static final Logger logger = Logger.getLogger(LoopyBeliefPropagation.class);
    // Type of inference: default is sum-product
    private InferenceType inferenceType = InferenceType.SUM_PRODUCT;
    // Flag if logspace should be used
    private boolean logSpace;
    // The initial message value: change this may get different results
    // because of the multiple local maxima of the objective function
    private double initialMessage = 1.0d; // Default value
    // The default update will iterate through variables. However,
    // iterating through factors may yield a different local maxima,
    // so it may need to try too.
    private boolean updateViaFactors;
    // Set dumping constant as described in the PGM book: page 408
    // As described at that page, the following is (1.0 - ramda)
    private double dumping = 0.0d; // Default there is no dumping
    // Enable convergence checking
    private boolean enableConvergenceCheck = true;
    
    /**
     * Default constructor.
     */
    public LoopyBeliefPropagation() {
    }
    
    public boolean isEnableConvergenceCheck() {
        return enableConvergenceCheck;
    }

    public void setEnableConvergenceCheck(boolean enableConvergenceCheck) {
        this.enableConvergenceCheck = enableConvergenceCheck;
    }

    public boolean isUpdateViaFactors() {
        return updateViaFactors;
    }
    
    public void setDumping(double dumping) {
        this.dumping = dumping;
    }
    
    public double getDumping() {
        return this.dumping;
    }
    
    /**
     * Set the message update via the factors. The default is via the variables.
     * Updating via factors may yield a different result because of multiple
     * local maximum of the joint energy functions.
     * @param updateViaFactors
     */
    public void setUpdateViaFactors(boolean updateViaFactors) {
        this.updateViaFactors = updateViaFactors;
    }
    
    public boolean getUpdateViaFactors() {
        return this.updateViaFactors;
    }

    /**
     * This should be the value in the probability space. So even if the log-space
     * is used, this value should not be logarithmized. The default value is 1.0.
     * Note: This method should be called before setUseLogSpace().
     * @param value
     */
    public void setInitialMessage(double value) {
        this.initialMessage = value;
    }
    
    public double getInitialMessage() {
        return this.initialMessage;
    }
    
    public void setUseLogSpace(boolean logSpace) {
        this.logSpace = logSpace;
    }
    
    public boolean getUseLogSpace() {
        return this.logSpace;
    }
    
    public void setInferenceType(InferenceType type) {
        this.inferenceType = type;
    }
    
    public InferenceType getInferenceType() {
        return this.inferenceType;
    }
    
    /**
     * The actual method to perform inference. This is a synchronized method so that only one
     * thread can run inference using this object.
     */
    public synchronized void runInference() throws InferenceCannotConvergeException {
        super.runInference(); // Do whatever the super class wants.
        truncateContinuousFactors();
        attachObservation();
        // Initialize messages
        initializeMessages(factorGraph);
        // Make a copy to avoid changing the original orders in the
        // passed factor graph object.
        List<Factor> factors = new ArrayList<Factor>(factorGraph.getFactors());
        List<Variable> variables = new ArrayList<Variable>(factorGraph.getVariables());
        // The max diff will be updated in method sendMessage. For the time being,
        // Set it as maximum so that the iteration can start
        maxDiff = Double.MAX_VALUE;
        iteration = 0;
        long time1 = System.currentTimeMillis();
        // Keep maxDiffList for checking if the inference cannot converge.
        List<Double> maxDiffList = new ArrayList<Double>();
        while (iteration <= maxIteration && maxDiff > tolerance) {
            maxDiff = 0.0d; // minimum. It should be reset in the sendMessage() method
            // The following message schedule is focusing on variables only.
            // First it update messages from factors to variables, and then
            // from variables to factors. This schedule proves a better convergence
            // for the FI network in the Ising model.
            if (updateViaFactors) {
                Collections.shuffle(factors);
                for (Factor factor : factors) {
                    for (Edge edge : factor.getOutEdges())
                        sendMessage(edge);
                    for (Edge edge : factor.getInEdges())
                        sendMessage(edge);
                }
            }
            else {
                Collections.shuffle(variables);
                for (Variable var : variables) {
                    for (Edge edge : var.getInEdges())
                        sendMessage(edge);
                    for (Edge edge : var.getOutEdges())
                        sendMessage(edge);
                }
            }
            iteration ++;
            if (debug)
                logger.info("Iteration: " + iteration + ", maxDiff: " + maxDiff);
            maxDiffList.add(maxDiff);
            if(!validateConverge(maxDiffList)) {
//                logger.error("Cannot converg...");
                // Don't forget remove observation. Otherwise, we will keep a modified FactorGraph.
                detachObservation();
                throw new InferenceCannotConvergeException("Inference for " + factorGraph + ": cannot converge.");
            }
        }
        long time2 = System.currentTimeMillis();
        if (debug)
            logger.info("Inference is done: " + iteration + ", maxDiff: " + maxDiff + ", using " + (time2 - time1) / 1000.0d + " seconds.");
        if (iteration > maxIteration) // No convergence. The client should be warned!
            logger.warn("Inferece for " + factorGraph + ": reach max iterations " + iteration + " with maxDiff " + maxDiff);
        calculateBeliefs(factorGraph);
        detachObservation();
        addBackContinuosFactors();
    }
    
    private boolean validateConverge(List<Double> maxDiffList) throws InferenceCannotConvergeException {
        if (!enableConvergenceCheck || maxDiffList.size() < 50)
            return true; // We will need at least 50 iterations
        int count = 0;
        for (int i = 0; i < maxDiffList.size() - 1; i++) {
            double diff1 = maxDiffList.get(i);
            double diff2 = maxDiffList.get(i + 1);
            if (diff2 > diff1)
                count ++;
        }
        return count < 10; // if there are 10 ups between two adjacent values. The cutoff is arbitrary.
    }
    
    
    /**
     * Calculate beliefs for variables in the FactorGraph object after LBP running.
     * @param fg
     */
    private void calculateBeliefs(FactorGraph fg) {
        // Calculate beliefs for variables
        for (Variable var : fg.getVariables()) {
            var.updateBelief(logSpace);
        }
        // Calculate beliefs for factors
        for (Factor factor : fg.getFactors()) {
            factor.updateBelief(logSpace);
        }
    }

    private void initializeMessages(FactorGraph fg) {
        // Refresh all edges in case they have been used previously
        for (Factor factor : fg.getFactors())
            factor.resetEdges();
        for (Variable var : fg.getVariables())
            var.resetEdges();
        // Generate in and out edges based on factors and variables
        for (Factor factor : fg.getFactors()) {
            initializeMessages(factor);
        }
    }

    private void initializeMessages(Factor factor) {
        for (Variable var : factor.getVariables()) {
            Edge factorToVarEdge = new Edge(factor, var);
            factorToVarEdge.initializeMessage(initialMessage, logSpace);
            Edge varToFactorEdge = new Edge(var, factor);
            varToFactorEdge.initializeMessage(initialMessage, logSpace);
            factor.addOutEdge(factorToVarEdge);
            factor.addInEdge(varToFactorEdge);
            var.addOutEdge(varToFactorEdge);
            var.addInEdge(factorToVarEdge);
        }
    }
    
    /**
     * Find the maximum assignment after a MAX_SUM inference is performed. The nature
     * of LBP determines that the returned result may not be an actual maximum, but 
     * should be so-called "Strong Local Maximum". See section 13.4 in the PGM book
     * for the detailed discussion, and also section 13.3.3.
     * @return
     */
    public Map<Variable, Integer> findMaximum() {
        if (inferenceType != InferenceType.MAX_PRODUCT)
            throw new IllegalArgumentException("findMaximum should be called after a MAX_PROD inference.");
        Map<Variable, Integer> varToAssign = new HashMap<Variable, Integer>();
        Set<Variable> variables = new HashSet<Variable>(factorGraph.getVariables());
        // First pass: If a Variable is unambiguous, we can choose its maximum state with
        // the largest value in its belief based on local optimality
        for (Variable var : variables) {
            double[] belief = var.getBelief();
            // If there is only one state, just pick it
            if (belief.length == 1) {
                varToAssign.put(var, belief.length - 1);
                continue;
            }
            int maxIndex = findMaxiumIndex(belief);
            if (maxIndex > -1)
                varToAssign.put(var, maxIndex);
        }
        variables.removeAll(varToAssign.keySet());
        if (variables.size() > 0 && debug) {
            logger.info("There are ambiguous variables: " + variables.size());
        }
        // Try to find max assignment for ambiguous variables by iterating factors
        for (Factor factor : factorGraph.getFactors()) {
            // Check if this factor should be used
            boolean isNeeded = false;
            for (Variable variable : factor.getVariables()) {
                if (variables.contains(variable)) {
                    isNeeded = true;
                    break;
                }
            }
            if (!isNeeded)
                continue;
            extractMaximumFromFactor(factor,
                                     varToAssign);
        }
        variables.removeAll(varToAssign.keySet());
        if (variables.size() > 0) {
            logger.error("There are variables that cannot be assigned: " + variables);
        }
        return varToAssign;
    }
    
    private void extractMaximumFromFactor(Factor factor,
                                          Map<Variable, Integer> varToAssign) {
        double[] belief = factor.getBelief();
        List<Integer> sortedIndices = sortIndices(belief);
        // We want to use the assignment with the maximum match to previously found
        // assignment via Variables only in case the factor is ambiguous
        int maxMatched = -1; // Use initialize value -1 so that we can get at least the first assignment 
                             // even if there is nothing known (aka no match at all).
        Map<Variable, Integer> maxFactorVarAssign = null;
        int matched = 0;
        for (int i = 0; i < sortedIndices.size(); i++) {
            if (i > 0) {
                // Check whether or not the same value is duplicated
                // If not, we should not check again
                int index0 = sortedIndices.get(i);
                int index1 = sortedIndices.get(i - 1);
                if (belief[index0] < belief[index1])
                    break;
            }
            int index = sortedIndices.get(i);
            Map<Variable, Integer> factorVarAssgn = factor.getAssignment(index);
            matched = 0;
            for (Variable var : factorVarAssgn.keySet()) {
                Integer known = varToAssign.get(var);
                if (known != null && known.equals(factorVarAssgn.get(var)))
                    matched ++;
            }
            if (matched > maxMatched) { 
                maxMatched = matched;
                maxFactorVarAssign = factorVarAssgn;
            }
        }
        // Theoretically assignments of variables here should be the same as in the previous
        // ones got from variables
        for (Variable var : maxFactorVarAssign.keySet()) {
            if (!varToAssign.containsKey(var))
                varToAssign.put(var, maxFactorVarAssign.get(var));
            // We can do this check if needed
            else if (varToAssign.get(var) != maxFactorVarAssign.get(var))
                logger.warn("Factor has a different assignment for Variable " + var + " in Factor " + factor);
        }
    }

    private List<Integer> sortIndices(final double[] belief) {
        // Need a sorted indices
        List<Integer> sortedIndices = new ArrayList<Integer>();
        for (int i = 0; i < belief.length; i++)
            sortedIndices.add(i);
        Collections.sort(sortedIndices, new Comparator<Integer>() {
            public int compare(Integer index1, Integer index2) {
                double value1 = belief[index1];
                double value2 = belief[index2];
                if (value1 > value2)
                    return -1;
                else if (value1 < value2)
                    return 1;
                return 0;
            }
        });
        return sortedIndices;
    }
    
    /**
     * Find the index of the maximum value in the belief. If the top two values
     * are the same, -1 is returned to indicate an ambiguous case.
     * @param belief
     * @return
     */
    private int findMaxiumIndex(double[] belief) {
        double max = Double.MIN_VALUE;
        int maxIndex = -1;
        double secondMax = Double.MIN_VALUE;
        for (int i = 0; i < belief.length; i++) {
            if (belief[i] > max) {
                max = belief[i];
                maxIndex = i;
            }
            else if (belief[i] > secondMax) {
                secondMax = belief[i];
            }
        }
        if (max == secondMax)
            return -1;
        return maxIndex;
    }
    
    private void sendMessage(Edge edge) {
        FGNode from = edge.getFromNode(); 
        FGNode to = edge.getToNode(); 
        double[] message = from.sendMessage(to, 
                                            inferenceType,
                                            logSpace);
//        System.out.println("Message before dumping: " + Arrays.toString(message));
        // Update messages
        double[] oldMessage = edge.getMessage();
        // Dumping is not supported for the time being. It doesn't help converge for some reason.
//        addDumping(message,
//                   oldMessage, 
//                   logSpace,
//                   to);
//        System.out.println("Message after dumping: " + Arrays.toString(message));
        double diff = 0.0d;
        for (int i = 0; i < message.length; i++) {
            if (Double.isNaN(message[i]))
                throw new IllegalStateException("A Message contains NaN: a possible numerical underflow occurs. Probably the log-space should be used for computation.");
//            if (Double.isInfinite(message[i])) // This should be allowed and hope to get back a meaningful value in the log-space.
//                throw new IllegalStateException("A Message contains Infinity: a possible numerical overflow occurs.");
//            if (logSpace) // Better to use the probability space to avoid infinity
//                diff = Math.abs(Math.exp(message[i]) - Math.exp(oldMessage[i]));
//            else
            // If the log space, the difference between two values when converging should
            // be exp(x) * dx, since exp(x)'s derivative is exp(x). Since x
            // is <= 0 (x is log of a probability), so exp(x) < 1, which makes
            // exp(x) * dx < dx. So if the log space is used and we use dx, instead of exp(x) * dx,
            // we may need several more steps. But we will avoid use exp(x), which is very slow. 
            // Overall we should be able to increase the performance.
            diff = Math.abs(message[i] - oldMessage[i]);
            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }
        edge.setMessage(message);
        //        logger.info("sendMessage: from " + from + " to " + to + " with message " + 
        //                message + " (old: " + oldMessage + "). Max: " + maxDiff);
    }
    
    private void addDumping(double[] newMessage,
                            double[] oldMessage,
                            boolean logSpace,
                            FGNode node) {
        if (logSpace) {
           node.convertLogToProb(newMessage);
           node.convertLogToProb(oldMessage);
        }
        for (int i = 0; i < newMessage.length; i++) {
            newMessage[i] = (1.0d - dumping) * newMessage[i] + dumping * oldMessage[i];
        }
        if (logSpace) {
            node.convertProbToLog(newMessage);
            node.convertProbToLog(oldMessage);
        }
    }

}
