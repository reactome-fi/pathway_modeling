/*
 * Created on Mar 23, 2017
 *
 */
package org.reactome.booleannetwork;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * A class to read in a Boolean network in a variety of formats.
 * @author gwu
 *
 */
public class BooleanNetworkReader {
    
    /**
     * Default constructor.
     */
    public BooleanNetworkReader() {
    }
    
    @Test
    public void testReadToyModel() throws IOException {
        String fileName = "results/BooleanNetwork/ToyModel.txt";
        BooleanNetwork network = readToyModel(fileName);
        System.out.println(network);
    }
    
    @Test
    public void testJAXB() throws JAXBException, IOException {
        JAXBContext context = JAXBContext.newInstance(BooleanNetwork.class);
        Marshaller marshaller = context.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        String fileName = "results/BooleanNetwork/ToyModel.txt";
        BooleanNetwork network = readToyModel(fileName);
        BooleanVariable var = network.getVariables().iterator().next();
        var.addProperty("reactomeId", "12345");
        System.out.println("Total variables: " + network.getVariables().size());
        System.out.println("Total relatioins: " + network.getRelations().size());
        System.out.println("\nWriting:");
        StringWriter writer = new StringWriter();
        marshaller.marshal(network, writer);
        System.out.println(writer.toString());
        
        System.out.println("\nReading...");
        Unmarshaller unmarshaller = context.createUnmarshaller();
        network = (BooleanNetwork) unmarshaller.unmarshal(new StringReader(writer.toString()));
        System.out.println("Name: " + network.getName());
        System.out.println("Total variables: " + network.getVariables().size());
        for (BooleanVariable var1 : network.getVariables()) {
            String reactomeId = var1.getProperty("reactomeId");
            if (reactomeId != null)
                System.out.println(var1.getName() + " has Reactome id: " + reactomeId);
        }
        System.out.println("Total relatioins: " + network.getRelations().size());
    }
    
    /**
     * This method is used to load a BooleanNetwork output from a CellNOptR ToyModel.
     * @param fileName
     * @return
     * @throws IOException
     */
    public BooleanNetwork readToyModel(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index;
        Map<String, BooleanVariable> nameToVar = new HashMap<>();
        BooleanNetwork network = new BooleanNetwork();
        network.setName(fileName);
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("=");
            String varName1 = line.substring(0, index);
            boolean isNegated = false;
            if (varName1.startsWith("!")) {
                varName1 = varName1.substring(1);
                isNegated = true;
            }
            String varName2 = line.substring(index + 1);
            BooleanVariable var1 = getVariable(varName1, nameToVar);
            BooleanVariable var2 = getVariable(varName2, nameToVar);
            BooleanRelation relation = new BooleanRelation();
            relation.setOutputVariable(var2);
            relation.addInputVariable(var1, isNegated);
            network.addRelation(relation);
        }
        fu.close();
        network.validateVariables();
        return network;
    }
    
    private BooleanVariable getVariable(String name,
                                        Map<String, BooleanVariable> nameToVar) {
        BooleanVariable var = nameToVar.get(name);
        if (var == null) {
            var = new BooleanVariable();
            var.setName(name);
            nameToVar.put(name, var);
        }
        return var;
    }
}
