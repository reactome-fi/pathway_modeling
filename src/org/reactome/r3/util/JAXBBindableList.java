/*
 * Created on Feb 11, 2015
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * This class wraps an ArrayList for easy JAXB binding.
 * @author gwu
 */
@XmlRootElement
public class JAXBBindableList<T> {
    private List<T> list;
    
    /**
     * Default constructor.
     */
    public JAXBBindableList() {
        list = new ArrayList<T>();
    }

    public List<T> getList() {
        return list;
    }

    public void setList(List<T> list) {
        this.list = list;
    }
    
    
}
