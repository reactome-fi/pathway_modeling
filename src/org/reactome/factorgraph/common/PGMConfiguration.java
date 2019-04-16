/*
 * Created on May 28, 2014
 *
 */
package org.reactome.factorgraph.common;

import java.io.FileInputStream;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.reactome.factorgraph.ExpectationMaximization;
import org.reactome.factorgraph.GibbsSampling;
import org.reactome.factorgraph.LoopyBeliefPropagation;

/**
 * This class is used to parse an XML based configuration file. This class should be used as 
 * singleton.
 * @author gwu
 */
@SuppressWarnings("unchecked")
public class PGMConfiguration {
    private static final Logger logger = Logger.getLogger(PGMConfiguration.class);
    // These names are used for central dogma nodes
    public static final String DNA = "DNA";
    public static final String mRNA = "mRNA";
    public static final String protein = "protein";
    // Member variables
    protected Map<DataType, double[]> typeToThreshold;
    protected Map<DataType, double[]> typeToFactorValues;
    protected Map<DataType, ObservationFactorHandler> typeToHandler;
    // In order to support gene expression in Gaussian distribution
    private boolean useGaussianForGeneExp;
    protected Map<DataType, String> typeToFileName;
    protected Map<String, String> properties;
    // Values configured for central dogma nodes
    protected double[] centralDogmaValues;
    protected int numberOfStates;
    // A flag to control if parameters should be learned
    private boolean learnParameters;
    private ExpectationMaximization em;
    // Set LoopyBeliefPropagation algorithm
    private LoopyBeliefPropagation lbp;
    // Gibbs Sampling
    private GibbsSampling gbs;
    
    protected PGMConfiguration() {
    }
    
    public boolean isUseGaussianForGeneExp() {
        return useGaussianForGeneExp;
    }

    public void setUseGaussianForGeneExp(boolean useGaussianForGeneExp) {
        this.useGaussianForGeneExp = useGaussianForGeneExp;
    }

    public Map<DataType, double[]> getTypeToThreshold() {
        return typeToThreshold;
    }
    
    public void setTypeToThreshold(DataType dataType, double[] thresholds) {
        if (typeToThreshold == null)
            typeToThreshold = new HashMap<DataType, double[]>();
        typeToThreshold.put(dataType, thresholds);
    }
    
    public Map<DataType, double[]> getTypeToFactorValues() {
        return typeToFactorValues;
    }
    
    public Map<DataType, String> getTypeToEvidenceFile() {
        return typeToFileName;
    }
    
    public Map<DataType, ObservationFactorHandler> getTypeToHandler() {
        return typeToHandler;
    }
    
    public Map<String, String> getProperties() {
        return this.properties;
    }
    
    public int getNumberOfStates() {
        return this.numberOfStates;
    }
    
    public void setNumberOfStates(int states) {
        this.numberOfStates = states;
    }
    
    /**
     * A flag indicating is parameters should be learned. The default is false.
     * @return
     */
    public boolean getLearnParameters() {
        return learnParameters;
    }
    
    public void setLearnParameters(boolean learn) {
        this.learnParameters = learn;
    }
    
    /**
     * Get a pre-configured ExpectationMaximization. The returned EM is cloned from a pre-configured one
     * to avoid thread issues.
     * @return
     */
    public ExpectationMaximization getEM() throws Exception {
        if (em == null)
            return null;
        // Make a copy
        ExpectationMaximization clone = (ExpectationMaximization) em.getClass().newInstance();
        clone.setDebug(em.getDebug());
        clone.setMaxIteration(em.getMaxIteration());
        clone.setTolerance(em.getTolerance());
        // Handle LBP
        clone.setInferenser(getLBP());
        return clone;
    }

    public void cloneLBP(LoopyBeliefPropagation target,
                         LoopyBeliefPropagation src) {
        target.setDebug(src.getDebug());
        target.setMaxIteration(src.getMaxIteration());
        target.setTolerance(src.getTolerance());
        target.setUseLogSpace(src.getUseLogSpace());
        target.setDumping(src.getDumping());
        target.setUpdateViaFactors(src.getUpdateViaFactors());
    }
    
    /**
     * A client should call this method to set up necessary configuration using a File.
     * @throws Exception
     */
    public void configure(String configFile) throws Exception {
        config(new FileInputStream(configFile));
    }
    
    /**
     * A client should call this method to set up necessary configuration using a generic InputStream object.
     * @param is
     * @throws Exception
     */
    public void config(InputStream is) throws Exception {
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(is);
        Element root = document.getRootElement();
        List<Element> children = root.getChildren();
        for (Element child : children) {
            String name = child.getName();
            if (name.equals("thresholds")) {
                parseThresholds(child);
            }
            else if (name.equals("dataFiles"))
                parseDataFiles(child);
            else if (name.equals("useGaussianForGeneExp"))
                parseUseGaussianForGeneExp(child);
            else if (name.equals("factorValues"))
                parseFactorValues(child);
            else if (name.equals("properties"))
                parseProperties(child);
            else if (name.equals("learnParameters"))
                parseLearnParameters(child);
            else if (name.equals("LoopyBeliefPropagation"))
                parseLBP(child);
            else if (name.equals("GibbsSampling"))
                parseGBS(child);
        }
    }
    
    private void parseUseGaussianForGeneExp(Element element) {
        String value = element.getText();
        if (value.equals("true"))
            useGaussianForGeneExp = true;
        else 
            useGaussianForGeneExp = false;
    }
    
    private void parseGBS(Element element) {
        gbs = new GibbsSampling();
        List<Element> list = element.getChildren();
        for (Element elm : list) {
            String name = elm.getName();
            String value = elm.getText();
            if (name.equals("debug"))
                gbs.setDebug(new Boolean(value));
            else if (name.equals("maxIteration"))
                gbs.setMaxIteration(new Integer(value));
            else if (name.equals("burnin"))
                gbs.setBurnin(new Integer(value));
            else if (name.equals("restart"))
                gbs.setRestart(new Integer(value));
        }
    }
    
    private void parseLearnParameters(Element element) throws Exception {
        String needToLearn = element.getAttributeValue("needTolearn");
        if (needToLearn.equalsIgnoreCase("true") || needToLearn.equals("1")) {
            learnParameters = true;
            // Get the defined class if any
            String clsName = element.getAttributeValue("class");
            if (clsName != null && clsName.length() > 0)
                em = (ExpectationMaximization) Class.forName(clsName).newInstance();
            else
                em = new ExpectationMaximization();
            List<Element> children = element.getChildren();
            for (Element child : children) {
                String name = child.getName();
                String value = child.getValue();
                if (name.equals("maxIteration"))
                    em.setMaxIteration(new Integer(value));
                else if (name.equals("tolerance"))
                    em.setTolerance(new Double(value));
                else if (name.equals("debug"))
                    em.setDebug(new Boolean(value));
            }
        }
        else
            learnParameters = false;
    }
    
    private void parseLBP(Element elm) {
        // Currently support LBP only
        LoopyBeliefPropagation lbp = new LoopyBeliefPropagation();
        parseLBPParameters(elm, lbp);
        this.lbp = lbp;
    }
    
    public LoopyBeliefPropagation getLBP() {
        if (lbp == null)
            return null;
        LoopyBeliefPropagation rtn = new LoopyBeliefPropagation();
        cloneLBP(rtn, lbp);
        return rtn;
    }
    
    public GibbsSampling getGibbsSampling() {
        if (gbs == null)
            return gbs;
        GibbsSampling rtn = new GibbsSampling();
        rtn.setDebug(gbs.getDebug());
        rtn.setMaxIteration(gbs.getMaxIteration());
        rtn.setRestart(gbs.getRestart());
        rtn.setBurnin(gbs.getBurnin());
        return rtn;
    }
    
    private void parseLBPParameters(Element elm, LoopyBeliefPropagation lbp) {
        List<Element> children = elm.getChildren();
        for (Element child : children) {
            String name = child.getName();
            String value = child.getValue();
            if (name.equals("debug"))
                lbp.setDebug(new Boolean(value));
            else if (name.equals("maxIteration"))
                lbp.setMaxIteration(new Integer(value));
            else if (name.equals("tolerance"))
                lbp.setTolerance(new Double(value));
            else if (name.equals("logSpace"))
                lbp.setUseLogSpace(new Boolean(value));
            else if (name.equals("updateViaFactors"))
                lbp.setUpdateViaFactors(new Boolean(value));
            else if (name.equals("dumping"))
                lbp.setDumping(new Double(value));
        }
    }
    
    private void parseThresholds(Element element) {
        List<Element> list = element.getChildren();
        typeToThreshold = new HashMap<DataType, double[]>();
        for (Element elm : list) {
            String type = elm.getAttributeValue("type");
            String value = elm.getAttributeValue("value");
            double[] thresholds = parseValues(value);
            typeToThreshold.put(DataType.valueOf(type),
                                thresholds);
        }
    }
    
    private void parseDataFiles(Element element) {
        typeToFileName = new HashMap<DataType, String>();
        List<Element> list = element.getChildren();
        for (Element elm : list) {
            String type = elm.getAttributeValue("type");
            String fileName = elm.getAttributeValue("value");
            typeToFileName.put(DataType.valueOf(type),
                               fileName);
            // Check if handler is specified
            Element handlerElm = elm.getChild("handler");
            if (handlerElm != null)
                parseDataHandler(type, handlerElm);
        }
    }
    
    private void parseDataHandler(String type,
                                  Element handlerElm) {
        String clsName = handlerElm.getAttributeValue("class");
        try {
            ObservationFactorHandler handler = (ObservationFactorHandler) Class.forName(clsName).newInstance();
            List<?> propElms = handlerElm.getChildren("property");
            if (propElms != null && propElms.size() > 0) {
                for (Object obj : propElms) {
                    Element propElm = (Element) obj;
                    String name = propElm.getAttributeValue("name");
                    String value = propElm.getAttributeValue("value");
                    String setMethodName = "set" + name.substring(0, 1).toUpperCase() + name.substring(1);
                    Method setMethod = handler.getClass().getMethod(setMethodName, String.class);
                    setMethod.invoke(handler, value);
                }
            }
            if (typeToHandler == null)
                typeToHandler = new HashMap<DataType, ObservationFactorHandler>();
            typeToHandler.put(DataType.valueOf(type),
                              handler);
        }
        catch(Exception e) {
            logger.error(e);
        }
    }
    
    private void parseFactorValues(Element element) {
        typeToFactorValues = new HashMap<DataType, double[]>();
        List<Element> list = element.getChildren();
        for (Element elm : list) {
            String type = elm.getAttributeValue("type");
            String value = elm.getAttributeValue("value");
            double[] factorValues = parseValues(value);
            typeToFactorValues.put(DataType.valueOf(type),
                                   factorValues);
        }
    }
    
    private double[] parseValues(String text) {
        String[] tokens = text.split(" ");
        double[] values = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++)
            values[i] = new Double(tokens[i]);
        return values;
    }
    
    private void parseProperties(Element element) {
        properties = new HashMap<String, String>();
        List<Element> list = element.getChildren();
        for (Element elm : list) {
            String name = elm.getName();
            String value = elm.getText();
            properties.put(name, value);
        }
    }
    
    /**
     * Different DataType has different thresholds for discretizing values.
     * This static method is the configuration for these threshold values.
     * @param type
     * @return
     */
    public double[] getThreshold(DataType type) {
        Map<DataType, double[]> typeToThreshold = getTypeToThreshold();
        if (typeToThreshold == null)
            throw new IllegalArgumentException("Uknown threshold valules for data type: " + type);
        double[] rtn = typeToThreshold.get(type);
        if (rtn == null)
            throw new IllegalArgumentException("Uknown threshold valules for data type: " + type);
        // There should be only one value configured
        return rtn;
    }
    
    /**
     * Get the pre-configured factor values for a passed DataType.
     * @param type
     * @return
     */
    public double[] getDataTypeValues(DataType type) {
        Map<DataType, double[]> typeToValues = getTypeToFactorValues();
        return typeToValues.get(type);
    }
    
    public double[] getCentralDogmaValues() {
        return this.centralDogmaValues;
    }
    
    public void setCentralDogmaValues(double[] values) {
        this.centralDogmaValues = values;
    }

}
