/*
 * Created on Oct 15, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.reactome.factorgraph.Factor;
import org.reactome.factorgraph.Variable;

/**
 * Use this class to check how many components may be contained in a set of factors. If more than one
 * component, link them together.
 * @author gwu
 *
 */
public class ComponentHelper {
    private final static Logger logger = Logger.getLogger(ComponentHelper.class);
    
    /**
     * Default constructor.
     */
    public ComponentHelper() {
    }
    
    /**
     * The client should call this method to ensure factors are linked together.
     * @param factors
     */
    public void ensureOneComponent(Set<Factor> factors) {
        List<Set<Factor>> components = checkComponents(factors);
        logger.info("Total components: " + components.size());
        if (components.size() == 1)
            return;
        int i = 0;
        for (Set<Factor> comp : components) {
            logger.info("Component " + i + ": " + comp.size());
            i++;
        }
        linkComponents(components, 
                       factors);
//        useBiggestComponent(components, factors);
    }
    
    private void useBiggestComponent(List<Set<Factor>> components,
                                     Set<Factor> factors) {
        for (int i = 1; i < components.size(); i++) {
            Set<Factor> comp = components.get(i);
            factors.removeAll(comp);
        }
    }
    
    /**
     * Add fake Factors to link all variables together into a single FactorGraph.
     * @param components
     */
    private void linkComponents(List<Set<Factor>> components,
                                Set<Factor> factors) {
        // The first component should be the biggest one. All other factors
        // should be linked to the first components randomly
        Set<Factor> first = components.get(0);
        for (int i = 1; i < components.size(); i++) {
            Set<Factor> second = components.get(i);
            Variable var1 = pickUpVar(first);
            Variable var2 = pickUpVar(second);
            Factor factor = createEqualFactor(var1, var2);
            factors.add(factor);
        }
    }
    
    /**
     * Just to create a factor with equal probability for all combination.
     * @param var1
     * @param var2
     * @return
     */
    private Factor createEqualFactor(Variable var1,
                                     Variable var2) {
        List<Variable> varList = new ArrayList<Variable>();
        varList.add(var1);
        varList.add(var2);
        List<Double> values = new ArrayList<Double>();
        for (int i = 0; i < var1.getStates() * var2.getStates(); i++)
            values.add(1.0d);
        Factor factor = new Factor();
        factor.setVariables(varList);
        factor.setValues(values);
        factor.setName("LinkBetweenComponents");
        return factor;
    }
    
    private Variable pickUpVar(Set<Factor> factors) {
        Set<Variable> variables = new HashSet<Variable>();
        for (Factor factor : factors) {
            variables.addAll(factor.getVariables());
        }
        List<Variable> list = new ArrayList<Variable>(variables);
        int index = (int) (Math.random() * list.size());
        return list.get(index);
    }
    
    /**
     * This helper method is used to check connected graph components.
     * @param factors
     * @param keyToVar
     * @return
     */
    private List<Set<Factor>> checkComponents(Set<Factor> factors) {
        List<Set<Factor>> components = new ArrayList<Set<Factor>>();
        for (Factor factor : factors) {
            Set<Factor> comp = new HashSet<Factor>();
            comp.add(factor);
            components.add(comp);
        }
        
        List<Set<Factor>> toBeRemoved = new ArrayList<Set<Factor>>();
        while (components.size() > 1) {
            for (int i = 0; i < components.size() - 1; i++) {
                toBeRemoved.clear();
                Set<Factor> comp1 = components.get(i);
                for (int j = i + 1; j < components.size(); j++) {
                    Set<Factor> comp2 = components.get(j);
                    if (hasSharedVariables(comp1, 
                                           comp2)) {
                        comp1.addAll(comp2);
                        toBeRemoved.add(comp2);
                    }
                }
                if (toBeRemoved.size() > 0) {
                    components.removeAll(toBeRemoved);
                    break; // Restart from zero since the comparison should be changed now.
                }
            }
            if (toBeRemoved.size() == 0)
                break;
        }
        // Sort it based on size
        Collections.sort(components, new Comparator<Set<Factor>>() {
            public int compare(Set<Factor> comp1, Set<Factor> comp2) {
                return comp2.size() - comp1.size();
            }
        });
        return components;
    }
    
    private boolean hasSharedVariables(Set<Factor> comp1, 
                                       Set<Factor> comp2) {
        for (Factor f1 : comp1) {
            for (Factor f2 : comp2) {
                if(hasSharedVariable(f1, f2))
                    return true;
            }
        }
        return false;
    }
    
    private boolean hasSharedVariable(Factor factor1, 
                                      Factor factor2) {
        for (Variable var : factor1.getVariables()) {
            if (factor2.getVariables().contains(var))
                return true;
        }
        return false;
    }
    
}
