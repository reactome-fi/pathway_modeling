/*
 * Created on Dec 7, 2012
 *
 */
package org.reactome.r3.util;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.database.EventCheckOutHandler;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.PersistenceManager;
import org.gk.persistence.XMLFileAdaptor;
import org.junit.Test;

/**
 * This class is used to link to Reactome pathways for the PARADIGM implementation.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactomeDBBridge {
    private static final Logger logger = Logger.getLogger(ReactomeDBBridge.class);
//    public static final String DUMP_DIR_NAME = "results/paradigm/ReactomePathways/";
    public static final String DUMP_DIR_NAME = "results/ReactomePathways/";
//    public static final String REACTOME_PATHWAY_FILE = DUMP_DIR_NAME + "ReactomePathways_Release42.rtpj";
//    public static final String REACTOME_PATHWAY_FILE = DUMP_DIR_NAME + "ReactomePathways_Release50.rtpj";
//    public static final String REACTOME_PATHWAY_FILE = DUMP_DIR_NAME + "ReactomePathways_Release59.rtpj";
    public static final String REACTOME_PATHWAY_FILE = DUMP_DIR_NAME + "ReactomePathways_Release63.rtpj";
    
    public ReactomeDBBridge() {
    }
    
    /**
     * Get the released pathways from a released database.
     * @return
     * @throws Exception
     */
    private Set<Long> getReleasedPathwayIds(MySQLAdaptor dba) throws Exception {
        // Get a list of top-level human pathways from the FrontPage instance
        Collection<GKInstance> c = dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        GKInstance frontPage = c.iterator().next();
        List<GKInstance> frontPageItem = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        Set<GKInstance> humanPathways = new HashSet<GKInstance>();
        Set<Long> dbIds = new HashSet<Long>();
        for (GKInstance inst : frontPageItem) {
            GKInstance species = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.species);
            if (species.getDisplayName().equals("Homo sapiens")) {
                dbIds.add(inst.getDBID());
            }
        }
        // Don't want to include these two human pathways
        // Disease
        dbIds.remove(1643685L);
        dbIds.remove(879392L); // Mycobacterium tuberculois biological processes
        return dbIds;
    }
    
    /**
     * Dump pathways in the Reactome db to a curator tool project
     * file so that they can be used in a cluster.
     * @throws Exception
     */
    @Test
    public void dumpPathwaysToFile() throws Exception {
        FileUtility.initializeLogging();
        logger.info("Starting dumping...");
        // Note: Have to use a slice database. Otherwise, there is too many instances to be checked out,
        // which takes a huge time and much large memory without any use!!!
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_slice_63", // reactome_59_plus_i can be used too to generate the same thing.
                                            "root", 
                                            "macmysql01");
        // Get a list of top-level human pathways from the FrontPage instance
        Set<Long> pathwayIds = getReleasedPathwayIds(dba);
        Set<GKInstance> humanPathways = new HashSet<GKInstance>();
        for (Long dbId : pathwayIds) {
            GKInstance inst = dba.fetchInstance(dbId);
            humanPathways.add(inst);
        }
        
        // For test
//        Long dbId = pathwayIds.iterator().next();
//        pathwayIds.clear();
//        pathwayIds.add(dbId);
        
        EventCheckOutHandler checkOutHelper = new EventCheckOutHandler();
//            logger.info("Handling " + inst + "...");
//            // This is for test
//            if (!inst.getDisplayName().equals("Apoptosis"))
//                continue;
//            logger.info("Checking out " + inst + "...");
            // This is needed by checkOutHelper
            XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
            PersistenceManager.getManager().setActiveFileAdaptor(fileAdaptor);
            PersistenceManager.getManager().setActiveMySQLAdaptor(dba);
            logger.info("Checking out Events...");
            checkOutHelper.checkOutEvents(humanPathways, fileAdaptor);
            // The following two classes should be taken care of in checkOtuEvent().
//            logger.info("Checking out CatalystActivities...");
//            checkOutAll(fileAdaptor, 
//                        dba, 
//                        ReactomeJavaConstants.CatalystActivity);
//            logger.info("Checking out Regulations...");
//            checkOutAll(fileAdaptor,
//                        dba,
//                        ReactomeJavaConstants.Regulation);
            // Make sure all needed instances are fully checked out
            logger.info("Checking out PhysicalEntities...");
            checkOutAll(fileAdaptor, 
                        dba, 
                        ReactomeJavaConstants.PhysicalEntity);
            logger.info("Checking out ReferenceSequences...");
            checkOutAll(fileAdaptor,
                        dba, 
                        ReactomeJavaConstants.ReferenceSequence);
            logger.info("Saving to file...");
            fileAdaptor.save(REACTOME_PATHWAY_FILE);
//            fileAdaptor.save(DUMP_DIR_NAME + inst.getDisplayName() + ".rtpj");
            logger.info("Saved!");
//        }
    }
    
    private void checkOutAll(XMLFileAdaptor fileAdaptor,
                             MySQLAdaptor dba,
                             String className) throws Exception {
        PersistenceManager manager = PersistenceManager.getManager();
        while (true) {
            Collection<GKInstance> c = fileAdaptor.fetchInstancesByClass(className);
            boolean hasShell = false;
            for (GKInstance inst : c) {
                if (inst.isShell()) {
                    GKInstance dbCopy = dba.fetchInstance(inst.getDBID());
//                    dba.fastLoadInstanceAttributeValues(dbCopy);
                    manager.updateLocalFromDB(inst, 
                                              dbCopy);
                    hasShell = true;
                }
            }
            if (!hasShell)
                break;
        }
    }
    
    public Set<GKInstance> getAllOutputs(GKInstance pathway) throws Exception {
        Set<GKInstance> rtn = new HashSet<GKInstance>();
        Set<GKInstance> events = InstanceUtilities.grepPathwayEventComponents(pathway);
        for (GKInstance event : events) {
            if (!event.getSchemClass().isValidAttribute(ReactomeJavaConstants.output))
                continue;
            List<GKInstance> output = event.getAttributeValuesList(ReactomeJavaConstants.output);
            rtn.addAll(output);
        }
        return rtn;
    }
    
    public Set<Long> getAllOutputIds(GKInstance pathway) throws Exception {
        Set<Long> ids = new HashSet<Long>();
        Set<GKInstance> outputs = getAllOutputs(pathway);
        for (GKInstance output : outputs)
            ids.add(output.getDBID());
        return ids;
    }
    
    @Test
    public void testGetAllOutputIds() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(REACTOME_PATHWAY_FILE);
        Long dbId = 452723L; 
        GKInstance pathway = fileAdaptor.fetchInstance(dbId);
        Set<Long> outputs = getAllOutputIds(pathway);
        System.out.println("Outputs for " + pathway + ": " + outputs.size());
        for (Long output : outputs) {
            GKInstance inst = fileAdaptor.fetchInstance(output);
            System.out.println(inst);
        }
    }
    
    @Test
    public void checkGeneRegulationInPathways() throws Exception {
    	MySQLAdaptor dba = new MySQLAdaptor("localhost",
    									    "gk_central_010914",
    									    "root", 
    										"macmysql01");
    	Collection<GKInstance> pds = dba.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram);
    	for (GKInstance pd : pds) {
    		List<GKInstance> pathways = pd.getAttributeValuesList(ReactomeJavaConstants.representedPathway);
    		if (pathways.size() != 1)
    			continue;
    		GKInstance pathway = pathways.get(0);
    		Set<GKInstance> events = InstanceUtilities.grepPathwayEventComponents(pathway);
    		int count = 0;
    		for (GKInstance event : events) {
    			if (!event.getSchemClass().isa(ReactomeJavaConstants.BlackBoxEvent)) {
    				continue;
    			}
    			// Want to check event's regulationss
    			Collection<GKInstance> regulators = InstanceUtilities.getRegulations(event); 
    			if (regulators == null || regulators.size() == 0)
    				continue;
    			count ++;
    		}
    		if (count > 0)
    			System.out.println(pathway + "\t" + events.size() + "\t" + count + "\t" + (double)count / events.size());
    	}
    }
    
    @Test
    public void testLoadFile() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(REACTOME_PATHWAY_FILE);
        System.out.println("Finish loading!");
        try {
            Thread.sleep(1000 * 30);
        }
        catch(InterruptedException e) {
            e.printStackTrace();
        }
    }
    
}
