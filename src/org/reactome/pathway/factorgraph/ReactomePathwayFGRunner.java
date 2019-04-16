/*
 * Created on Jun 11, 2012
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;

import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.render.Renderable;
import org.gk.render.RenderablePathway;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.factorgraph.*;
import org.reactome.factorgraph.common.DataType;
import org.reactome.factorgraph.common.ObservationFileLoader;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.JAXBBindableList;
import org.reactome.r3.util.ReactomeDBBridge;

/**
 * This class is used to run PARADIGM for the converted Reactome FactorGraph. Most of the code here is ported
 * from PARADIGM c++ code.
 * @author gwu
 *
 */
public class ReactomePathwayFGRunner {
    private static final Logger logger = Logger.getLogger(ReactomePathwayFGRunner.class);
    // Used to load the pathways to be processed
    private PersistenceAdaptor adaptor;
    // Specify a list of pathways that should be run
    private List<Long> chosenPathwayIds;
    // Result directory
    private String resultDir = "tmp/";
    // Flag to indicate there is no learning and inferencing is needed. Just to check
    // converting pathway to factor graphs.
    private boolean convertToFGOnly = false;
    // Cached loaded data files
    private Map<DataType, ObservationFileLoader.ObservationData> typeToData;
    // Flag to show the inference is for perturbation analysis
    private boolean isForPerturbation;
    
    public static void main(String[] args) {
        ReactomePathwayFGRunner runner = null;
        if (args.length > 0) {
            if (args[0].equals("checkCPUNumber")) {
                System.out.println("Total CPU/Core: " + Runtime.getRuntime().availableProcessors());
                System.exit(0);
            }
            else if (args[0].matches("parallel=(\\d+)")) {
                int index = args[0].indexOf("=");
                int number = new Integer(args[0].substring(index + 1));
                runner = new ReactomePathwayFGRunner();
                FileUtility.initializeLogging();
                // There should be a second argument
                if (args.length < 2) {
                    logger.error("No remoteCommand has been provided!");
                    System.exit(1);
                }
                DRMAAJobScheduler scheduler = new DRMAAJobScheduler(args[1]);
                try {
                    String dirName = runner.resultDir;
                    if (args.length > 2) // Directory should be able to specify too.
                        dirName = args[2];
                    scheduler.runParallelParadigm(runner,
                                                  number,
                                                  dirName);
                }
                catch(Exception e) {
                    logger.error("Error in PARADIGMRunner.main()", e);
                }
                System.exit(0);
            }
            else {
                initializeLogging(args[0]);
                // This should be the directory pass
                runner = new ReactomePathwayFGRunner();
                runner.setResultDir(args[0]);
            }
        }
        try {
            if (runner == null) {
                runner = new ReactomePathwayFGRunner();
                FileUtility.initializeLogging();
                logger.info("Initialized PARADIGMRunner().");
            }
            runner.runAllPathways();
            logger.info("Finishing the main() method!");
            // Need to shutdown logging. For programmatically created FileAppender,
            // This is needed. Otherwise, the logging may be truncated.
            LogManager.shutdown();
        }
        catch (Exception e) {
            logger.error(e.getMessage(), e);
        }
    }
    
    public ReactomePathwayFGRunner() {
    }
    
    public void setResultDir(String dir) {
        logger.info("setResultDir: " + dir);
        this.resultDir = dir;
        // Make sure there is a file separator at the end of the resultDir
        if (!resultDir.endsWith(File.separator)) 
            resultDir = resultDir + File.separator;
        // Check if there is a Pathway list file
        File file = new File(dir, "Pathways.txt");
        if (file.exists()) {
            FileUtility fu = new FileUtility();
            try {
                fu.setInput(file.getAbsolutePath());
                String line = null;
                List<Long> dbIds = new ArrayList<Long>();
                while ((line = fu.readLine()) != null) {
                    String[] tokens = line.split("\t");
                    dbIds.add(new Long(tokens[0]));
                }
                fu.close();
                setChosenPathwayIds(dbIds);
            }
            catch(IOException e) {
                logger.error("Error in setResultDir(): ", e);
            }
        }
    }
    
    public void setChosenPathwayIds(List<Long> ids) {
        logger.info("setChoosenPathwayIds(): " + ids.size() + " in total.");
        this.chosenPathwayIds = ids;
    }
    
    public List<String> getNameForEscapeList() {
        // Add a list of escape small molecules
        String[] escapeNames = new String[] {
                "ATP",
                "ADP",
                "Pi",
                "H2O",
                "GTP",
                "GDP",
                "CO2",
                "H+"
        };
        List<String> escapeList = Arrays.asList(escapeNames);
        return escapeList;
    }
    

    /**
     * @param string
     * @param priorFactors
     * @param clamped
     * @param builder
     */
    private void calculateIPA(Map<Variable, double[]> varToPrior,
                              Map<Variable, double[]> varToPosterior,
                              Map<Variable, Double> varToIPA) {
        for (Variable var : varToPrior.keySet()) {
            double[] prior = varToPrior.get(var);
            double[] posterior = varToPosterior.get(var);
            double ipa = IPACalculator.calculateIPA(prior, 
                                                    posterior);
            varToIPA.put(var, ipa);
        }
    }
    
    private static void initializeLogging(String dir) {
        try {
            File file = new File(dir, "logging.txt");
            PatternLayout pattern = new PatternLayout("%d{ISO8601} [%t] %-5p %c %x - %m%n");
            FileAppender fa = new FileAppender(pattern, file.getAbsolutePath());
            fa.setName("FileLogger");
            fa.setAppend(true);
            fa.setBufferedIO(true);
            fa.setBufferSize(1024);
            fa.setThreshold(Level.INFO);
            fa.activateOptions();
            Logger.getRootLogger().addAppender(fa);
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
    @Test
    public void testGetPathways() throws Exception {
        List<GKInstance> pathways = getPathwayList();
        for (GKInstance pathway : pathways) {
            System.out.println(pathway);
        }
    }
    
    /**
     * Test to see how observations are generated.
     * @throws Exception
     */
    @Test
    public void testCheckObservations() throws Exception {
        FileUtility.initializeLogging();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "gk_current_ver50", 
                                            "root", 
                                            "macmysql01");
        Long dbId = 15869L; // Metabolism of nucleotides
        GKInstance pathway = dba.fetchInstance(dbId);
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        ConvertedFactorGraph cfg = convertPathway(pathway, converter);
        System.out.println("Total observations: " + cfg.observations.size());
        Set<Variable> geneExpVars = new HashSet<Variable>();
        Set<Variable> cnvVars = new HashSet<Variable>();
        for (Observation obs : cfg.observations) {
            Map<Variable, Integer> varToAssign = obs.getVariableToAssignment();
            for (Variable var : varToAssign.keySet()) {
                if (var.getName().endsWith(DataType.mRNA_EXP.toString()))
                    geneExpVars.add(var);
                else if (var.getName().endsWith(DataType.CNV.toString()))
                    cnvVars.add(var);
            }
        }
        System.out.println("Total exp vars: " + geneExpVars.size());
        System.out.println("Total cnv vars: " + cnvVars.size());
        String fileName = "tmp/" + pathway.getDisplayName() + "_GeneExp.txt";
        outputObservations(cfg, geneExpVars, fileName);
        fileName = "tmp/" + pathway.getDisplayName() + "_CNV.txt";
        outputObservations(cfg, cnvVars, fileName);
    }

    private void outputObservations(ConvertedFactorGraph cfg,
                                    Set<Variable> geneExpVars, String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (Variable var : geneExpVars) {
            builder.append("\t");
            builder.append(var.getName());
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (Observation obs : cfg.observations) {
            Map<Variable, Integer> varToAssign = obs.getVariableToAssignment();
            builder.append(obs.getName());
            for (Variable var : geneExpVars) {
                Integer state = varToAssign.get(var);
                builder.append("\t").append(state == null ? "na" : state);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkPathwayInference() throws Exception {
        FileUtility.initializeLogging();
        Long dbId = 15869L; // Metabolism of nucleotides: have different results in different runs
        
        List<GKInstance> pathways = getPathwayList();
        GKInstance pathway = null;
        for (GKInstance inst : pathways) {
            if (inst.getDBID().equals(dbId)) {
                pathway = inst;
                break;
            }
        }
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        ConvertedFactorGraph cfg = convertPathway(pathway, converter);
        List<Variable> pathwayVars = getPathwayVariables(pathway,
                                                         cfg.instToVar,
                                                         cfg.fg);
        // Use LBP as the default inferencer.
        // Be careful with the threading issue: we will use a new lbp object always
        // to avoid overwrite any member properties.
        LoopyBeliefPropagation lbp = PathwayPGMConfiguration.getConfig().getLBP();
        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.setFactorGraph(cfg.fg);
        // Just in case GibbsSampling is used. Since this is a very light initiailization,
        // we can afford to create an object that may not be used to make coding clean
        GibbsSampling gbs = PathwayPGMConfiguration.getConfig().getGibbsSampling();
        gbs.setFactorGraph(cfg.fg);
        // Save the output into a file
        String fileName = resultDir + InteractionUtilities.getFileNameFromInstance(pathway) + ".txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        int count = 0;
        exportHeaders(pathwayVars, cfg.instToVar, fu);
        for (int i = 0 ; i < 10; i++) {
            Observation observation = cfg.observations.get(0);
            logger.info("Performing inference on sample " + observation.getName() + " for " + pathway + "...");
            Map<Variable, double[]> varToPrior = new HashMap<Variable, double[]>();
            // We only need to calculate prior once
            logger.info("Performing inference on prior for " + pathway + "...");
            performInference(lbp,
                             gbs, 
                             pathwayVars,
                             varToPrior,
                             observation);
            Map<Variable, Double> varToValue = new HashMap<Variable, Double>();
            for (int j = 0; j < 3; j++) {
                for (Variable var : varToPrior.keySet()) {
                    double[] values = varToPrior.get(var);
                    varToValue.put(var, values[j]);
                }
                exportIPA(varToValue, 
                          fu, 
                          pathwayVars, 
                          "Prior_" + j);
            }
        }
        fu.close();
    }
    
    @Test
    public void testDumpAllFactorGraphs() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_51_plus_i", 
                                            "root", 
                                            "macmysql01");
        setAdaptor(dba);
        List<GKInstance> pathways = getPathwayList();
        System.out.println("Total pathways: " + pathways.size());
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        long time1 = System.currentTimeMillis();
        int count = 1;
        JAXBBindableList<FactorGraph> fgList = new JAXBBindableList<FactorGraph>();
        for (GKInstance pathway : pathways) {
            System.out.println(count +": " + pathway);
            FactorGraph fg = converter.convertPathway(pathway);
            fgList.getList().add(fg);
            count ++;
            if (count == 4)
                break;
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1) + " ms");
        System.out.println("Total memory: " + Runtime.getRuntime().totalMemory() / (1024 * 1024) + " mb");
        // Export into a file
        // Have to list both classes here.
        JAXBContext context = JAXBContext.newInstance(JAXBBindableList.class, FactorGraph.class);
        Marshaller marshaller = context.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
//        marshaller.marshal(fgList, new File("tmp/AllFGs.xml"));
        File file = new File("tmp/AllFGs.zip");
        FileOutputStream fos = new FileOutputStream(file);
        ZipOutputStream zos = new ZipOutputStream(fos);
        ZipEntry entry = new ZipEntry("AllFgs.xml");
        zos.putNextEntry(entry);
        marshaller.marshal(fgList, zos);
        zos.close();
        
        // To read
//        JAXBContext context = JAXBContext.newInstance(JAXBBindableList.class, FactorGraph.class);
//        File file = new File("tmp/AllFGs.zip");
//        Unmarshaller unmarshaller = context.createUnmarshaller();
//        FileInputStream fis = new FileInputStream(file);
//        ZipInputStream zis = new ZipInputStream(fis);
//        ZipEntry entry = zis.getNextEntry(); // Have to call this method
//        JAXBBindableList<FactorGraph> list = (JAXBBindableList<FactorGraph>) unmarshaller.unmarshal(zis);
//        System.out.println("Read: " + list.getList().size());
    }
    
    @Test
    public void testRunSinglePathway() throws Exception {
        isForPerturbation = true;
        FileUtility.initializeLogging();
        chosenPathwayIds = new ArrayList<Long>();
        
        // Used to work with a specific pathway
        Long dbId = 2032785L; // YAP1 and WWTR1-stimulated gene expression
        dbId = 1257604L; // PIP3 activates AKT signaling: a pathway has a known complex loop
        dbId = 5673001L; // RAF/MAP kinase cascade
////        Long dbId = 2028269L; // Signaling by Hippo
////        dbId = 400206L; // Regulation of Lipid Metabolism by... (not the same factor numbers)
////        dbId = -1L;
////        dbId = 1592389L; // Activation of Matrix Metalloproteinases
////        dbId = 186712L; // Regulation of beta-cell development: a gene regulatory network
////        dbId = 400206L; // Regulation of Lipid Metabolism by Peroxisome proliferator-activated receptor alpha (PPARalpha): Another gene regulatory network
////        dbId = 381340L; // Transcriptional Regulation of White Adipocyte Differentiation
//////        dbId = 381038L; // Activation of Chaperone Genes by XBP1(S)
//        dbId = 535734L; // Fatty acid, triacylglycerol, and ketone body metabolism: took over 12 hours for the TCGA BRCA data.
//      
//        dbId = 452723L; // Transcriptional Regulation of Pluripotent Stem Cells
//        dbId = 1980143L; // Signaling by NOTCH1
//        // A very slow process
//        dbId = 1483206L; // Glycerophospholipid biosynthesis
//        dbId = 1650814L; // Collagen biosynthesis and modifying enzymes
//        dbId = 73923L; 
        // [Pathway:1834949] Cytosolic sensors of pathogen-associated DNA: Got NaN using the log-space
//        dbId = 1834949L;
//        // [Pathway:111885] Opioid Signalling: cannot converge.
//        dbId = 111885L;
//        // [Pathway:400206] Regulation of lipid metabolism by Peroxis
//        dbId = 400206L;
//        // [Pathway:1059683] Interleukin-6 signaling
//        dbId = 1059683L;
//        // [Pathway:390918] Peroxisomal lipid metabolism: OutOfMemory: GC limit
//        dbId = 390918L;
//        // Passive transport by Aquaporins: Infinity IPA
//        dbId = 432047L;
//        // [Pathway:75153] Apoptotic execution  phase: null exception
//        dbId = 75153L;
//        
//        dbId = 15869L;
//        
//        dbId = 1483206L; // [Pathway:1483206] Glycerophospholipid biosynthesis\
//        dbId = 201451L; // [Pathway:201451] Signaling by BMP
                
//        chosenPathwayIds.add(dbId);
        
//        adaptor = new MySQLAdaptor("localhost",
//                                   "reactome_59_plus_i", 
//                                   "root", 
//                                   "macmysql01");
//        GKInstance pathway = adaptor.fetchInstance(dbId);
        
        // For a simple test pathway
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
//        fileAdaptor.setSource("results/PerturbationAnalysis/ChainPathway.rtpj");
        fileAdaptor.setSource("results/PerturbationAnalysis/ChainPathwayWithLoop.rtpj");
        Collection<GKInstance> pathways = fileAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
                                                                         ReactomeJavaConstants.name,
                                                                         "=",
                                                                         "Chain Pathway");
        GKInstance pathway = pathways.iterator().next();
        
        RunPathwayThread runner = new RunPathwayThread(pathway);
        // Force it to run in the same thread without generating another thread
        runner.run();
        
        
//        runAllPathways();
    }
    
    @Test
    public void testRunAllPathways() throws Exception {
//        convertToFGOnly = true;
        FileUtility.initializeLogging();
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_41_plus_i", 
//                                            "root", 
//                                            "macmysql01");
//        adaptor = dba;
        
        runAllPathways();
    }
    
    private void loadDataFiles() throws IOException {
        // Load data files
        Map<DataType, String> typeToFile = PathwayPGMConfiguration.getConfig().getTypeToEvidenceFile();
        if (typeToFile != null) {
            ObservationFileLoader observationLoader = new ObservationFileLoader();
            observationLoader.setPGMConfiguration(PathwayPGMConfiguration.getConfig());
            for (DataType type : typeToFile.keySet()) {
                String file = typeToFile.get(type);
                Map<String, Map<String, Float>> sampleToGeneToValue = observationLoader.loadObservationData(file, 
                                                                                                              type);
            }
            typeToData = observationLoader.getLoadedData();
        }
    }
    
    public void runAllPathways() throws Exception {
        List<GKInstance> pathwayList = getPathwayList();
        logger.info("Total pathways: " + pathwayList.size());
        if (chosenPathwayIds != null && chosenPathwayIds.size() > 0) {
            // Do a simple filtering
            for (Iterator<GKInstance> it = pathwayList.iterator(); it.hasNext();) {
                GKInstance inst = it.next();
                if (!chosenPathwayIds.contains(inst.getDBID()))
                    it.remove();
            }
        }
        
        // Just for test
//        Set<GKInstance> tmp = MathUtilities.randomSampling(pathwayList, 10);
//        pathwayList = new ArrayList<GKInstance>(tmp);
        
        logger.info("Subject to PARADIGM: " + pathwayList.size());
        // Pre-load data
        logger.info("Loading data files...");
        loadDataFiles();
        logger.info("Loading data is done.");
        // The following statements are used to run all pathways using multiple threads
        int numberOfThread = getNumberOfThreads(pathwayList);
        List<RunPathwayThread> threads = new ArrayList<RunPathwayThread>();
        while (pathwayList.size() > 0 || hasRunningThreads(threads)) {
//            logger.info("PathwayList size: " + pathwayList.size());
            // Define a thread for running pathway analysis
            for (int i = threads.size(); i < numberOfThread; i++) {    
                if (pathwayList.size() == 0)
                    break;
                // Popout the first pathway in the list
                GKInstance pathway = pathwayList.remove(0);
                RunPathwayThread t = createPathwayThread(pathway);
                t.start();
                threads.add(t);
            }
            // Wait for a little while
            try {
                Thread.sleep(1000); // 1 second
            }
            catch(InterruptedException e) {
                logger.error("runAllPathways: " + e.getMessage(), e);
            }
            removeFinishedThreads(threads);
        }
    }
    
    private boolean hasRunningThreads(List<RunPathwayThread> threads) {
        for (Iterator<RunPathwayThread> it = threads.iterator(); it.hasNext();) {
            RunPathwayThread t = it.next();
            if (!t.isDone()) {
                return true;
            }
        }
        return false;
    }
    
    private void removeFinishedThreads(List<RunPathwayThread> threads) {
        for (Iterator<RunPathwayThread> it = threads.iterator(); it.hasNext();) {
            RunPathwayThread t = it.next();
            if (!t.isDone()) {
                continue;
            }
            it.remove();
        }
    }
    
    private RunPathwayThread createPathwayThread(GKInstance pathway) {
        RunPathwayThread t = new RunPathwayThread(pathway);
        return t;
    }
    
    /**
     * A helper method to calculate how many threads can be created based on the total CPU in the machine.
     * @param pathwayList
     * @return
     */
    private int getNumberOfThreads(List<GKInstance> pathwayList) {
        int cpu = Runtime.getRuntime().availableProcessors();
        logger.info("Total CPU cores: " + cpu);
        // Want to use 80% cpu cores to spawning threads
        int threadNumber = (int) (cpu * 0.80d);        // Use a thread safe list
        if (threadNumber > 16)
            threadNumber = 16; // In order to control memory usage
        if (threadNumber > pathwayList.size())
            threadNumber = pathwayList.size();
        logger.info("Total planned number of threads: " + threadNumber);
        return threadNumber;
    }
    
    /**
     * There are many variables converted from instances contained by pathway participants (e.g.
     * members in an EntitySet). We don't want to calculate IPAs for these values.
     * @param pathway
     * @param instToVar
     * @return
     */
    List<Variable> getPathwayVariables(GKInstance pathway,
                                       Map<GKInstance, Variable> instToVar,
                                       FactorGraph fg) throws Exception {
        Set<GKInstance> containedInsts = InstanceUtilities.grepPathwayParticipants(pathway);
        List<Variable> vars = new ArrayList<Variable>();
        for (GKInstance inst : instToVar.keySet()) {
            if (containedInsts.contains(inst))
                vars.add(instToVar.get(inst));
        }
        // In case some variables may not be in the final factor graph
        vars.retainAll(fg.getVariables());
        // Sort based on names
        Collections.sort(vars, new Comparator<Variable>() {
            public int compare(Variable var1, Variable var2) {
                return var1.getName().compareTo(var2.getName());
            }
        });
        return vars;
    }
    
    private void runPathway(GKInstance pathway,
                            PathwayToFactorGraphConverter converter) throws Exception {
        long time1 = System.currentTimeMillis();
        ConvertedFactorGraph cfg = convertPathway(pathway, converter);
        FactorGraph fg = cfg.fg;
        Map<DataType, SharedEMFactors> typeToSharedFactors = cfg.typeToFactors;
        logger.info("FactorGraph for " + pathway + ": " + fg.getFactors().size() + " factors, " 
                + fg.getVariables().size() + " variables.");
        logger.info("FactorGraph is a tree: " + fg.isTree());
        if (convertToFGOnly)
            return;
        // The following two local variables may not be used if learning works fine.
        // Otherwise, we have to use the original assigned values to avoid a set of
        // half-learned parameter values, which are better than carefully tuned original
        // values.
        if (cfg.typeToFactors != null && cfg.typeToFactors.size() > 0) {
            List<EMFactor> emFactors = new ArrayList<EMFactor>(cfg.typeToFactors.values());
            Map<EMFactor, double[]> originalFactorToValues = recordOriginalValues(emFactors);
            try {
                performLearning(pathway,
                                converter,
                                cfg);
            }
            catch(InferenceCannotConvergeException e) {
                logger.error("Learning aborted for " + pathway + ": " + e.getMessage(), e);
                // Rollover
                recoverOriginalValues(emFactors,
                                      originalFactorToValues);
            }
        }
        // Though we may be able to learn the parameters, we should try to perform inference still.
        // Get variables we need to calculate IPAs
        performInference(pathway, 
                         converter,
                         cfg);
        long time2 = System.currentTimeMillis();
        logger.info("Time for " + pathway + ": " + (time2 - time1) / 1000.0d + " seconds.");
    }
    
    /**
     * The following method is synchronized so that a copy of ConvertedFactorGraph can be generated in a thread-safe
     * way. Attributes in the returned ConvertedFactorGraph object should be copies of individual values to avoid
     * overwriting from different threads.
     * @param pathway
     * @param converter
     * @return
     * @throws Exception
     */
    synchronized ConvertedFactorGraph convertPathway(GKInstance pathway, 
                                                     PathwayToFactorGraphConverter converter) throws Exception {
        FactorGraph fg = converter.convertPathway(pathway);
        Map<DataType, SharedEMFactors> typeToSharedFactors = converter.getTypeToSharedFactors();
        // Make a copy to avoid overwrite by other thread
        Map<DataType, SharedEMFactors> clone = null;
        if (typeToSharedFactors != null)
            clone = new HashMap<DataType, SharedEMFactors>(typeToSharedFactors);
        ConvertedFactorGraph rtn = new ConvertedFactorGraph();
        rtn.fg = fg;
        rtn.observations = converter.getObservationLoader().getObservations();
        rtn.instToVar = new HashMap<GKInstance, Variable>(converter.getInstToVarMap());
        rtn.typeToFactors = clone;
        return rtn;
    }
    
    /**
     * A helper method for performing parameter learning. The map of typeToSharedFactors
     * should NOT be got from the passed converter to avoid thread issues. The passed
     * typeToSharedFactors should be used.
     * @param pathway
     * @param converter
     * @param fg
     * @throws Exception
     */
    void performLearning(GKInstance pathway,
                         PathwayToFactorGraphConverter converter,
                         ConvertedFactorGraph cfg) throws Exception {
        if (!PathwayPGMConfiguration.getConfig().getLearnParameters())
            return; // There is no need to learn parameters
        long time1 = System.currentTimeMillis();
        logger.info("Starting learning parameters for pathway: " + pathway + "...");
        ExpectationMaximization em = PathwayPGMConfiguration.getConfig().getEM();
        em.setEvidences(Observation.convertListFromNumberToInteger(cfg.observations));
        Map<DataType, SharedEMFactors> typeToSharedFactors = cfg.typeToFactors;
        em.learn(cfg.fg, 
                 new ArrayList<EMFactor>(typeToSharedFactors.values()));
        for (DataType type : typeToSharedFactors.keySet()) {
            SharedEMFactors factors = typeToSharedFactors.get(type);
            logger.info("Learned parameters: " + type + ": " + Arrays.toString(factors.getValues()));
        }
        long time2 = System.currentTimeMillis();
        logger.info("Done learning for pathway: " + pathway + ": "+ (time2 - time1) / 1000.0d + " seconds.");
    }
    
    private void recoverOriginalValues(List<EMFactor> emFactors,
                                       Map<EMFactor, double[]> originalFactoToValues) {
        for (EMFactor factor : emFactors) {
            double[] values = originalFactoToValues.get(factor);
            if (values == null)
                continue;
            factor.setValues(values);
        }
    }
    
    private Map<EMFactor, double[]> recordOriginalValues(List<EMFactor> emFactors) {
        Map<EMFactor, double[]> factorToValues = new HashMap<EMFactor, double[]>();
        for (EMFactor factor : emFactors) {
            double[] values = factor.getValues();
            if (values == null)
                continue;
            double[] copy = Arrays.copyOf(values, values.length);
            factorToValues.put(factor, copy);
        }
        return factorToValues;
    }
    
    private Map<Variable, GKInstance> generateVarToInst(Map<GKInstance, Variable> instToVar) {
        Map<Variable, GKInstance> varToInst = new HashMap<Variable, GKInstance>();
        for (GKInstance inst : instToVar.keySet())
            varToInst.put(instToVar.get(inst), inst);
        return varToInst;
    }
    
    private List<Variable> getProteinVariables(FactorGraph fg) {
        List<Variable> list = new ArrayList<Variable>();
        Map<String, Variable> nameToVar = new HashMap<String, Variable>();
        for (Variable var : fg.getVariables())
            nameToVar.put(var.getName(), var);
        for (String name : nameToVar.keySet()) {
            if (name.endsWith("_protein")) {
                list.add(nameToVar.get(name));
            }
        }
        return list;
    }

    /**
     * A helper method for performing inference.
     * @param pathway
     * @param converter
     * @param fg
     * @throws Exception
     * @throws IOException
     */
    void performInference(GKInstance pathway,
                          PathwayToFactorGraphConverter converter,
                          ConvertedFactorGraph cfg) throws Exception, IOException {
        long time1 = System.currentTimeMillis();
        logger.info("Starting inference for pathway: " + pathway + "...");
        List<Variable> pathwayVars = getPathwayVariables(pathway,
                                                         cfg.instToVar,
                                                         cfg.fg);
//        List<Variable> pathwayVars = getProteinVariables(cfg.fg);
//        logger.info("Total variables to be output for " + pathway + ": " + pathwayVars.size());
        
        // Use LBP as the default inferencer.
        // Be careful with the threading issue: we will use a new lbp object always
        // to avoid overwrite any member properties.
        LoopyBeliefPropagation lbp = PathwayPGMConfiguration.getConfig().getLBP();
        // As of December 18, 2014, use MAX_PRODUCT for all pathway PGM inference.
//        lbp.setInferenceType(InferenceType.MAX_PRODUCT);
        lbp.setFactorGraph(cfg.fg);
        // Just in case GibbsSampling is used. Since this is a very light initiailization,
        // we can afford to create an object that may not be used to make coding clean
        GibbsSampling gbs = PathwayPGMConfiguration.getConfig().getGibbsSampling();
        gbs.setFactorGraph(cfg.fg);
        
        // Save the output into a file
        String fileName = resultDir + InteractionUtilities.getFileNameFromInstance(pathway) + ".txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        int count = 0;
        exportHeaders(pathwayVars, cfg.instToVar, fu);
        
        // We only need to calculate prior once
        Map<Variable, double[]> varToPrior = new HashMap<Variable, double[]>();
//        logger.info("Performing inference on prior for " + pathway + "...");
//        performInference(lbp,
//                         gbs, 
//                         pathwayVars,
//                         varToPrior,
//                         null);
//        logger.info("Done inference on prior for " + pathway);
//        if (varToPrior.size() == 0) {
//            fu.close(); // We will get an empty file
//            return; // Cannot calculate IPAs
//        }
//        System.out.println("Variable\tPrior");
//        for (Variable var : lbp.getFactorGraph().getVariables()) {
//            System.out.println(var.getName() + "\t" + var.getBelief()[0]);
//        }
//        exportInferenceResults(varToPrior, fu, pathwayVars, "Prior");
        
        // Just a test
        Observation<Number> observation = new Observation<Number>();
        for (Variable var : lbp.getFactorGraph().getVariables()) {
            if (var.getName().endsWith("_protein"))
                observation.addAssignment(var, 0);
        }
        logger.info("Performing inference on base for " + pathway + "...");
        performInference(lbp,
                         gbs, 
                         pathwayVars,
                         varToPrior,
                         observation);
        logger.info("Done inference on base for " + pathway);
        
        System.out.println("Variable\tBase");
        for (Variable var : lbp.getFactorGraph().getVariables()) {
            System.out.println(var.getName() + "\t" + var.getBelief()[0] + "\t" + var.getBelief()[1]);
        }
        
        // Start posterior inference
        Map<Variable, double[]> varToPosterior = new HashMap<Variable, double[]>();
        Map<Variable, Double> varToIPA = new HashMap<Variable, Double>();
        
//        String varName = "TBX5_protein";
//        varName = "WWTR1 [nucleoplasm]_1629775";
//        varName = "WWTR1_protein";
//        String varName = "PTEN [cytosol]_199420";
//        String varName = "Activator:PI3K [plasma membrane]_2316432";
//        varName = "TORC2 complex [cytosol]_198626";
        
        // Intrinsic Pathway for Apoptosis
//        String varName = "RAS GEFs [plasma membrane]_5672601";
        
        // For a test chain graph
        String varName = "Input11_-3";
        
        Variable var = lbp.getFactorGraph().getVariable(varName);
        observation.addAssignment(var, 1.0d);
        
//        for (Observation observation : cfg.observations) {
//            observation = cfg.observations.get(0);
            logger.info("Performing inference on sample " + observation.getName() + " for " + pathway + "...");
            // Prior inference
            varToPosterior.clear();
            performInference(lbp, 
                             gbs, 
                             pathwayVars, 
                             varToPosterior,
                             observation);
            
            System.out.println("Variable\tPosterior");
            for (Variable var1 : lbp.getFactorGraph().getVariables()) {
                System.out.println(var1.getName() + "\t" + var1.getBelief()[0] + "\t" + var1.getBelief()[1]);
            }
            System.out.println("\nVariable\tLogRatio");
            for (Variable var1 : pathwayVars) {
                double logRatio = IPACalculator.calculateLogRatio(varToPrior.get(var1)[1],
                                                                  var1.getBelief()[1]);
                System.out.println(var1.getName() + "\t" + logRatio);
            }
            System.out.println();
            exportInferenceResults(varToPosterior, fu, pathwayVars, "Posterior");
            varToPosterior.clear();
            
            if (varToPosterior.size() > 0) {
                calculateIPA(varToPrior,
                             varToPosterior,
                             varToIPA);
                exportIPA(varToIPA, 
                          fu, 
                          pathwayVars, 
                          observation.getName());
            }
            count ++;
            logger.info("Done inference for observation: " + observation.getName() + " for " + pathway);
//            if (count == 10)
//                break;
//        }
        fu.close();
        long time2 = System.currentTimeMillis();
        logger.info("Done inference for pathway: " + pathway + ": "+ (time2 - time1) / 1000.0d + " seconds.");
    }
    
    private void performInference(LoopyBeliefPropagation lbp,
                                  GibbsSampling gibbs,
                                  List<Variable> pathwayVars,
                                  Map<Variable, double[]> varToProbs,
                                  Observation observation) throws InferenceCannotConvergeException {
        try {
//            logger.info("Used Inferencer: " + lbp.getClass().getName());
            if (observation == null)
                lbp.clearObservation();
            else
                lbp.setObservation(observation);
            lbp.runInference();
            recordBeliefs(pathwayVars, varToProbs);
        }
        catch(InferenceCannotConvergeException e) {
            if (observation == null)
                logger.warn("LBP cannot converge for prior!");
            else
                logger.warn("LBP cannot converge for " + observation.getName());
            // Check if MAX_PRODUCT is used. If true, don't use Gibbs
            if (lbp.getInferenceType() == InferenceType.MAX_PRODUCT)
                return; // Gibbs sampling cannot generate marginal max probabilities
            if (observation == null)
                logger.warn("Switch to Gibbs for prior!");
            else
                logger.warn("Switch to Gibbs for " + observation.getName());
            // Try to remedy once
            if (observation == null)
                gibbs.clearObservation();
            else
                gibbs.setObservation(observation.getVariableToAssignment());
            gibbs.runInference();
            recordBeliefs(pathwayVars, varToProbs);
        }
    }

    private void exportHeaders(List<Variable> vars,
                               Map<GKInstance, Variable> instToVar,
                               FileUtility fu) throws IOException {
        Map<Variable, GKInstance> varToInst = generateVarToInst(instToVar);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (Variable var : vars) {
            GKInstance inst = varToInst.get(var);
            if (inst == null)
                builder.append("\t").append(var.getName());
            else
                builder.append("\t").append(inst.getDisplayName());
        }
        fu.printLine(builder.toString());
    }
    
    private void exportIPA(Map<Variable, Double> varToIPA,
                           FileUtility fu,
                           List<Variable> orderedVars,
                           String sample) throws IOException {
        StringBuilder builder = new StringBuilder();
        builder.append(sample);
        for (Variable var : orderedVars) {
            builder.append("\t").append(varToIPA.get(var));
        }
        fu.printLine(builder.toString());
    }
    
    private void exportInferenceResults(Map<Variable, double[]> varToValues,
                                        FileUtility fu,
                                        List<Variable> orderedVars,
                                        String sample) throws IOException {
        StringBuilder builder = new StringBuilder();
        builder.append(sample);
        for (Variable var : orderedVars) {
            double[] values = varToValues.get(var);
            builder.append("\t").append(Arrays.toString(values));
        }
        fu.printLine(builder.toString());
    }
    
    /**
     * Record inferred beliefs to avoid overwriting by further inference.
     * @param vars
     * @return
     */
    private void recordBeliefs(List<Variable> vars,
                               Map<Variable, double[]> varToBelief) {
        varToBelief.clear();
        for (Variable var : vars) {
            double[] belief = var.getBelief();
            double[] copy = new double[belief.length];
            System.arraycopy(belief,
                             0, 
                             copy, 
                             0,
                             belief.length);
            varToBelief.put(var, copy);
        }
    }
    
    @Test
    public void testGetPathwayList() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_47_plus_i", 
                                            "root", 
                                            "macmysql01");
        adaptor = dba;
        List<GKInstance> pathways = getPathwayList();
        System.out.println("Total pathways: " + pathways.size());
    }

    public List<GKInstance> getPathwayList() throws Exception {
        List<GKInstance> pathways = queryPDPathway();
        return pathways;
    }
    
    public void setAdaptor(PersistenceAdaptor adaptor) {
        this.adaptor = adaptor;
    }
    
    private void initReactomeAdaptor() throws Exception {
        String reactomeFileName = PathwayPGMConfiguration.getConfig().getProperties().get("ReactomeFile");
        if (reactomeFileName == null || reactomeFileName.length() == 0)
            reactomeFileName = ReactomeDBBridge.REACTOME_PATHWAY_FILE;
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(reactomeFileName);
        adaptor = fileAdaptor;
    }

    /**
     * Get a list of pathways that have been drawn in pathway diagram. Only human pathways with
     * ReacctionlikeEvent instances drawn are used.
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    private List<GKInstance> queryPDPathway() throws Exception {
        // Try once
        if (adaptor == null) {
            initReactomeAdaptor();
        }
        if (adaptor == null) 
            throw new IllegalStateException("The data source for pathways has not been assigned!");
        SchemaClass cls = adaptor.getSchema().getClassByName(ReactomeJavaConstants.PathwayDiagram);
        Collection<GKInstance> pds = adaptor.fetchInstancesByClass(cls);
        DiagramGKBReader diagramReader = new DiagramGKBReader();
        List<GKInstance> rtn = new ArrayList<GKInstance>();
        for (GKInstance pd : pds) {
            List<GKInstance> pathways = pd.getAttributeValuesList(ReactomeJavaConstants.representedPathway);
            if (pathways == null || pathways.size() == 0)
                continue;
            GKInstance species = (GKInstance) pathways.get(0).getAttributeValue(ReactomeJavaConstants.species);
            if (!species.getDBID().equals(48887L)) { // For human pathways only
                continue;
            }
            // Check if this is a super pathway
            RenderablePathway diagram = diagramReader.openDiagram(pd);
            if (diagram.getComponents() == null || diagram.getComponents().size() == 0)
                continue;
            boolean isSuperPathway = true;
            for (Object o : diagram.getComponents()) {
                Renderable r = (Renderable) o;
                if (r.getReactomeId() == null)
                    continue;
                GKInstance inst = adaptor.fetchInstance(r.getReactomeId());
                if (inst != null && inst.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                    isSuperPathway = false;
                    break;
                }
            }
            if (!isSuperPathway) {
                // No disease pathway for the time being
                for (GKInstance pathway : pathways) {
                    if (pathway.isShell())
                        continue; // Don't use shell pathway.
                    GKInstance disease = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.disease);
                    if (disease != null)
                        continue;
                    rtn.add(pathway);
                }
            }
        }
        InstanceUtilities.sortInstances(rtn);
        return rtn;
    }
    
    /**
     * A simple data structure to be used in a multiple threading environment.
     * @author gwu
     *
     */
    class ConvertedFactorGraph {
        
        FactorGraph fg;
        Map<DataType, SharedEMFactors> typeToFactors;
        List<Observation<Number>> observations;
        Map<GKInstance, Variable> instToVar;
        
    }
    
    /**
     * A simple thread for call runPathway() using multiple threads.
     */
    private class RunPathwayThread extends Thread {
        private GKInstance pathway;
        private boolean isDone;
        
        public RunPathwayThread(GKInstance pathway) {
            this.pathway = pathway;
        }
        
        public void run() {
            try {
                // Each thread needs its own converter since concerter
                // has some thread-not-safe variables.
                PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
                converter.setForPerturbation(isForPerturbation);
                converter.setNamesForEscape(getNameForEscapeList());
                if (typeToData != null)
                    converter.getObservationLoader().setLoadedData(typeToData);
                runPathway(pathway, converter);
                isDone = true;
                // Call gc and hope to get some memory back to avoid
                // GC over limit error.
                System.gc();
            }
            catch(Exception e) { // want to catch any error
                isDone = true; // Have to set this flag. Otherwise, the main thread will not stop.
                logger.error("runPathway: " + e.getMessage() + " for " + pathway, e);
            }
            catch(Error e) { // Handle some errors like OutOfMemory to avoid a thread is hanging there
                isDone = true;
                logger.error("runPathway: " + e.getMessage() + " for " + pathway, e);
            }
        }
        
        public boolean isDone() {
            return this.isDone;
        }
    }
    
}
