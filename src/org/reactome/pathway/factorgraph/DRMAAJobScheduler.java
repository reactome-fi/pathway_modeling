/*
 * Created on Jan 29, 2014
 *
 */
package org.reactome.pathway.factorgraph;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.JobInfo;
import org.ggf.drmaa.JobTemplate;
import org.ggf.drmaa.Session;
import org.ggf.drmaa.SessionFactory;
import org.gk.model.GKInstance;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used as a client to sun the OCIR Sun grid based on DRMAA Java API.
 * For details, see http://gridscheduler.sourceforge.net/howto/drmaa_java.html.
 * @author gwu
 *
 */
public class DRMAAJobScheduler {
    private static final Logger logger = Logger.getLogger(DRMAAJobScheduler.class);
    private String remoteCommand;
    
    /**
     * The default constructor.
     */
    public DRMAAJobScheduler(String remoteCommand) {
        this.remoteCommand = remoteCommand;
    }
    
    public void runParallelParadigm(ReactomePathwayFGRunner runner,
                                    int totalNodes,
                                    String dir) throws Exception {
        if (remoteCommand == null)
            throw new IllegalStateException("remoteCommand has not been provided!");
        long start = System.currentTimeMillis();
        logger.warn("Starting runParallelParadigm at " + new Date());
        logger.warn("Specified total nodes: " + totalNodes);
        // Get the whole list of pathways.
        List<GKInstance> pathways = runner.getPathwayList();
        splitPathways(pathways, 
                      totalNodes,
                      dir);
//        if (true) {
//            System.out.println("Total nodes: " + totalNodes);
//            System.out.println("Dir: " + dir);
//            return;
//        }
        boolean isOk = runJobs(totalNodes, dir);
        if (!isOk) {
            logger.error("Jobs cannot be finished all! Some jobs may failed. Check the above job information!");
            return;
        }
        long end = System.currentTimeMillis();
        logger.warn("Total time used: " + (end - start) / (1000 * 60) + " minutes");
        logger.warn("Ending runParallelParadigm at " + new Date());
    }
    
    private boolean runJobs(int totalNodes, String dir) throws DrmaaException, IOException {
        SessionFactory factory = SessionFactory.getFactory();
        final Session session = factory.getSession();
        session.init(null);
        Thread t = new Thread() {
            public void run() {
                try {
                    session.exit();
                }
                catch(DrmaaException e) {
                    logger.error("Error in exit of session: " + e.getMessage(), e);
                }
            }
        };
        Runtime.getRuntime().addShutdownHook(t);
        JobTemplate jt = session.createJobTemplate();
        List<String> jobIds = new ArrayList<String>();
        for (int i = 0; i < totalNodes; i++) {
            File nodeDir = getNodeDir(dir, i);
            jt.setRemoteCommand(remoteCommand);
            jt.setArgs(new String[] {nodeDir.getAbsolutePath()});
            // Need to override the default setting of DRMAA in order to send jobs
            // Don't forget the following parameters used in qsub!
            jt.setNativeSpecification("-V -w n -l h_vmem=16G");
            String jobId = session.runJob(jt);
            logger.info(jobId + " has been submitted.");
            jobIds.add(jobId);
        }
        session.deleteJobTemplate(jt);
        session.synchronize(jobIds,
                            Session.TIMEOUT_WAIT_FOREVER,
                            false);
        boolean finsihedAll = true;
        for (String jobId : jobIds) {
            JobInfo info = session.wait(jobId,
                                        Session.TIMEOUT_WAIT_FOREVER);
            outputJobInfo(info);
            if (!info.hasExited()) {
                finsihedAll = false;
            }
        }
        return finsihedAll;
    }
    
    private void outputJobInfo(JobInfo info) throws DrmaaException {
        if (info.wasAborted())
            logger.warn("Job " + info.getJobId() + " never ran");
        else if (info.hasExited())
            logger.warn("Job " + info.getJobId() + " finished regularly with exit status " + info.getExitStatus());
        else if (info.hasSignaled())
            logger.warn("Job " + info.getJobId() + " finished due to signal " + info.getTerminatingSignal());
        else
            logger.warn("Job " + info.getJobId() + " finished with unclear conditions");
        logger.warn("Job Usage for " + info.getJobId() + ": ");
        @SuppressWarnings("rawtypes")
        Map rmap = info.getResourceUsage();
        for (Object key : rmap.keySet()) {
            Object value = rmap.get(key);
            logger.warn(" " + key + "=" + value);
        }
    }
    
    /**
     * Get a dedicated file directory for a specified node in the cluster.
     * @param dir
     * @param nodeIndex
     * @return
     * @throws IOException
     */
    private File getNodeDir(String dir,
                            int nodeIndex) throws IOException {
        File nodeDir = new File(dir, "Node" + nodeIndex);
        if (!nodeDir.exists())
            nodeDir.mkdir();
        return nodeDir;
    }

    private void splitPathways(List<GKInstance> pathways,
                               int totalNodes,
                               String dir) throws Exception {
        logger.info("Total pathways: " + pathways.size());
        // Split pathways into pieces
        Collections.shuffle(pathways);
        // Use the ceil method to make sure all pathways have been assigned
        int sizeInNode = (int) Math.ceil((double)pathways.size() / totalNodes);
        int start = 0;
        int end = start + sizeInNode;
        for (int i = 0; i < totalNodes; i++) {
            if (end > pathways.size())
                end = pathways.size();
            if (start >= end)
                break;
            List<GKInstance> subList = pathways.subList(start, end);
            // Output this list
            File nodeDir = getNodeDir(dir, i);
            File file = getPathwayFile(nodeDir);
            outputPathways(subList, file);
            start = end;
            end = start + sizeInNode;
        }
    }

    private File getPathwayFile(File nodeDir) {
        File file = new File(nodeDir, "Pathways.txt");
        return file;
    }
    
    private void outputPathways(List<GKInstance> pathways, File file) throws IOException {
        logger.info("Pathways in " + file.getAbsolutePath() + ": " + pathways.size());
        FileUtility fu = new FileUtility();
        fu.setOutput(file.getAbsolutePath());
        for (GKInstance pathway : pathways) {
            fu.printLine(pathway.getDBID() + "\t" + pathway.getDisplayName());
        }
        fu.close();
    }
    
}
