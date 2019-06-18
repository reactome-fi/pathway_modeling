/*
 * A BooleanNetwork should have a set of BooleanVariables and BooleanRelations. This is the class
 * should be used to model the whole pathway.
 * Created on Mar 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
@XmlAccessorType(XmlAccessType.FIELD)
public class BooleanNetwork {
	@XmlElement(name="variable")
	private Set<BooleanVariable> variables;
	@XmlElement(name="relation")
	private Set<BooleanRelation> relations;
	private String name;
	
	public BooleanNetwork() {
	}

	public Set<BooleanVariable> getVariables() {
		return variables;
	}

	public Set<BooleanRelation> getRelations() {
		return relations;
	}
	
	public Map<String, BooleanVariable> getNameToVar() {
	    Map<String, BooleanVariable> nameToVar = new HashMap<>();
	    for (BooleanVariable var : variables) {
	    		// The following exception may be thrown. The client code should never
	    		// use name to identify a BooleanVariable.
	        if (nameToVar.containsKey(var.getName()))
	            throw new IllegalStateException(var.getName() + " is duplicated in more than one variable!");
	        nameToVar.put(var.getName(), var);
	    }
	    return nameToVar;
	}
	
	/**
	 * Find a BooleanVariable for the passed varName
	 * @param varName
	 * @return
	 */
	public BooleanVariable find(String varName) {
		Optional<BooleanVariable> optional = variables.stream().filter(var -> var.getName().equals(varName)).findAny();
		return optional.get();
	}
	
	/**
	 * This method should be called after a network has been created and all relations 
	 * have been added to make sure variables are correct.
	 */
	public void validateVariables() {
	    // Need to extract all variables for this relation
	    if (variables == null)
	        variables = new HashSet<>();
	    if (relations == null)
	        relations = new HashSet<>();
	    for (BooleanRelation relation : relations) {
	        if (relation.getInputVariables() != null) {
	            for (BooleanVariable var : relation.getInputVariables())
	                variables.add(var);
	        }
	        if (relation.getOutputVariable() != null)
	            variables.add(relation.getOutputVariable());
	    }
	    resetIds();
	}
	
	private void resetIds() {
	    // Set ids
	    int index = 0;
	    if (variables != null) {
	        for (BooleanVariable var : variables)
	            var.setId(index ++);
	    }
	    // Reset relations: make sure relations and variables have different ids.
	    // Otherwise, JAXB unmarshalling will be wrong
	    if (relations != null) {
	        for (BooleanRelation relation : relations)
	            relation.setId(index ++);
	    }
	}

	public void setRelations(Set<BooleanRelation> relations) {
		this.relations = relations;
	}
	
	public void addRelation(BooleanRelation relation) {
		if (relations == null)
			relations = new HashSet<BooleanRelation>();
		relations.add(relation);
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("Name: " + name);
		if (variables != null) {
			builder.append("Variables: ").append(variables.size()).append("\n");
			for (BooleanVariable variable : variables)
				builder.append(variable).append("\n");
		}
		if (relations != null) {
			builder.append("Relations: ").append(relations.size()).append("\n");
			for (BooleanRelation relation : relations)
				builder.append(relation).append("\n");
		}
		return builder.toString();
	}

}
