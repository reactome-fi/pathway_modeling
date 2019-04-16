/*
 * Created on Apr 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.util.HashMap;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlID;

/**
 * An abstract superclass for relations and variables.
 * @author gwu
 *
 */
@XmlAccessorType(XmlAccessType.FIELD)
public abstract class BooleanObject {
    @XmlID
    protected String id;
    protected String name;
    // Customized information
    protected HashMap<String, String> properties;
    
    /**
     * Default constructor.
     */
    public BooleanObject() {
    }

    public void addProperty(String name, String value) {
        if (properties == null)
            properties = new HashMap<>();
        properties.put(name, value);
    }

    public String getProperty(String name) {
        if (properties == null)
            return null;
        return properties.get(name);
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
    
}
