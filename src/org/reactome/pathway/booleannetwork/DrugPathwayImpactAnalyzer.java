/*
 * Created on Sep 20, 2018
 *
 */
package org.reactome.pathway.booleannetwork;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.HibernateUtil;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.ProcessRunner;

import edu.ohsu.bcb.druggability.dataModel.Drug;
import edu.ohsu.bcb.druggability.dataModel.ExpEvidence;
import edu.ohsu.bcb.druggability.dataModel.Interaction;
import smile.validation.AdjustedRandIndex;

/**
 * Perform some type of systematic analysis for drugs based ob Boolean impact analysis results.
 * @author wug
 *
 */
public class DrugPathwayImpactAnalyzer {
//    private final Double IMPACT_CUTOFF = 1.0E-3;
    private final Double IMPACT_CUTOFF = 0.01d;
    private final String DIR = "results/BooleanNetwork/drugs/";
    private final FileUtility fu = new FileUtility();
    
    public DrugPathwayImpactAnalyzer() {
    }
    
    private Session getSession() throws Exception {
        String configFileName = "lib/drug/drugHibernate.cfg.xml";
        File configFile = new File(configFileName);
        SessionFactory sf = HibernateUtil.getSessionFactory(configFile);
        Session session = sf.openSession();
        return session;
    }
    
    /**
     * Use the following command to perform MCL clustering:
     * mcl {in_file} --abc -I 2 --force-connected=y -o {out_file}
     * Use --force_connected=y to make sure nodes in clusters are 
     * connected since the network is a bipartite graph.
     * Note: Use 2 for inflation. 3.5 is too big and generate too
     * many clusters.
     * @throws Exception
     */
    @Test
    public void convertMCLResultsToAttFile() throws Exception {
        String fileName = DIR + "Targetome_Drug_Pathway_Interaction_01_100218_mcl_I_2_out_102418.txt";
        String outFileName = DIR + "Targetome_Drug_Pathway_Interaction_01_100218_mcl_I_2_out_102418.att";
        fu.setInput(fileName);
        fu.setOutput(outFileName);
        fu.printLine("Node\tmodule");
        String line = null;
        int c = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            c ++;
            for (String token : tokens)
                fu.printLine(token + "\t" + c);
        }
        fu.close();
    }
    
    private Map<String, Integer> loadMCLDrugClusters(String fileName, Set<String> drugs) throws IOException {
        fu.setInput(fileName);
        int cluster = 0;
        String line = null;
        Map<String, Integer> drugToCluster = new HashMap<>();
        while ((line = fu.readLine()) != null) {
//            System.out.println(cluster + ": " + line);
            String[] tokens = line.split("\t");
            cluster ++;
            for (String token : tokens) {
                if (drugs.contains(token))
                    drugToCluster.put(token, cluster);
            }
        }
        fu.close();
        return drugToCluster;
    }
    
    @Test
    public void testLoadMCLDrugClusters() throws IOException {
        String fileName = DIR + "Targetome_Drug_Pathway_Interaction_01_100218_no_weight_mcl_I_2_out_102518.txt";
        Map<String, Integer> drugToCluster = loadMCLDrugClusters(fileName, new HashSet<>());
        
    }
    
    @Test
    public void compareTwoNetworkFiles() throws IOException {
        String fileName = "Targetome_Drug_Pathway_Interaction_01_100218_011819.txt";
        Set<String> interactions1 = Files.lines(Paths.get(DIR, fileName))
                .map(line -> line.replaceAll("\"", ""))
                .map(line -> line.split("\t"))
                .map(tokens -> tokens[0] + "\t" + tokens[1])
                .collect(Collectors.toSet());
        System.out.println("File1: " + fileName);
        System.out.println("Interactions1: " + interactions1.size());
        fileName = "Targetome_Drug_Pathway_Interaction_100218_011819.txt";
        Set<String> interactions2 = Files.lines(Paths.get(DIR, fileName))
                .map(line -> line.replaceAll("\"", ""))
                .map(line -> line.split("\t"))
                .map(tokens -> tokens[0] + "\t" + tokens[1])
                .collect(Collectors.toSet());
        System.out.println("File2: " + fileName);
        System.out.println("Interactions2: " + interactions2.size());
        Set<String> shared = InteractionUtilities.getShared(interactions1, interactions2);
        System.out.println("Shared: " + shared.size());
        interactions1.removeAll(shared);
        System.out.println("Unique in interactions1: " + interactions1.size());
        if (interactions1.size() < 10)
            interactions1.forEach(System.out::println);
        interactions2.removeAll(shared);
        System.out.println("Unique in interactions2: " + interactions2.size());
    }
    
    @Test
    public void permutationTestCompareCancerRxDrugsInClusters() throws Exception {
        // Get all drugs
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        Set<String> drugNames = drugs.stream().map(drug -> drug.getDrugName()).collect(Collectors.toSet());
        session.close();

        Map<String, String> cancerrxDrugToTargetomeDrug = loadCancerRxToTargetomeMap();

        String mclFile = "Targetome_Drug_Pathway_Interaction_100218_011819_mcl_I_2_out.txt";
        Map<String, Integer> targetomeDrugInMCLClusters = loadMCLDrugClusters(DIR + mclFile,
                                                                              drugNames);
        int cancerrxClusterNumber = 5;
        String cancerrxFileName = DIR + "cancerrx/cancerrx.k=" + cancerrxClusterNumber + ".consensusClass.csv";
        Map<String, Integer> cancerrxDrugToCluster = loadDrugToConsensusCluster(cancerrxFileName);
        Map<String, Integer> targetomeDrugInCancerrxClusters = new HashMap<>();
        Map<String, Integer> targetomeDrugToCancerrxId = new HashMap<>();
        cancerrxDrugToTargetomeDrug.forEach((key, value) -> {
            Integer cluster = cancerrxDrugToCluster.get(key);
            int drugId = new Integer(key.split("_")[1]);
            if (targetomeDrugToCancerrxId.containsKey(value)) {
                if (drugId > targetomeDrugToCancerrxId.get(value))
                    return;
            }
            targetomeDrugInCancerrxClusters.put(value, cluster);
            targetomeDrugToCancerrxId.put(value, drugId);
        });
        System.out.println("MCL file: " + mclFile);
        System.out.println("Cancerrx file: " + cancerrxFileName);
        double adjustRandIndex = calculateAdjustRandIndex(targetomeDrugInMCLClusters, 
                                                          targetomeDrugInCancerrxClusters);
        System.out.println(adjustRandIndex);
        // Perform permutation test
        int permutation = 10000;
        List<Double> randomARIs = new ArrayList<>();
        for (int i = 0; i < permutation; i++) {
            Map<String, Integer> randomMCLDrugToCluster = MathUtilities.permutate(targetomeDrugInMCLClusters);
            double randomARI = calculateAdjustRandIndex(randomMCLDrugToCluster, 
                                                        targetomeDrugInCancerrxClusters);
            randomARIs.add(randomARI);
        }
        System.out.println("\nResults from " + permutation + " permutations:");
        randomARIs.stream().sorted(Comparator.reverseOrder()).forEach(System.out::println);
    }
    
    @Test
    public void checkDrugClusterOverlaps() throws Exception {
        // Some file variables
        String[] mclFiles = {
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_mcl_out.txt",
//                "Targetome_Drug_Target_102518_030719_mcl_out.txt",
//                "Targetome_Drug_Target_Weight_102518_030719_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032019_mcl_out.txt"
                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Target_102518_L01XE_032119_032119_mcl_out.txt",
                "Targetome_Drug_Target_Weight_102518_L01XE_032119_032119_mcl_out.txt"
        };

        // Get all drug names for MCL clustering results
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        Map<String, String> drugNameToATCCode = drugs.stream().collect(Collectors.toMap(drug -> drug.getDrugName(),
                                                                                        drug -> drug.getAtcClassID()));
        Set<String> drugNames = drugs.stream().map(drug -> drug.getDrugName()).collect(Collectors.toSet());
        session.close();


        for (String mclFile : mclFiles) {
            System.out.println("MCL file: " + mclFile);
//            String cancerrxTitle = "cancerrx_hc_030819_1/cancerrx_hc_030819_1.k=";
//            int cancerrxClusterNumber = 5;
            
            String cancerrxTitle = "cancerrx_hc_kinase_inhibitors_032119/cancerrx_hc_kinase_inhibitors_032119.k=";
            int cancerrxClusterNumber = 3;
            
            String cancerrxFileName = cancerrxTitle + cancerrxClusterNumber + ".consensusClass.csv";

             // Load two clustering results
            Map<String, Integer> targetomeDrugToMCLCluster = loadMCLDrugClusters(DIR + mclFile,
                                                                                 drugNames);
            Map<String, Integer> cancerrxDrugToCluster = loadDrugToConsensusCluster(DIR + cancerrxFileName);
            Map<String, Integer> targetomeDrugToCancerrxCluster = mapCancerRxClustersToTargetomeClusters(cancerrxDrugToCluster);

            // Check overlapping
            Set<String> sharedDrugs = InteractionUtilities.getShared(targetomeDrugToCancerrxCluster.keySet(),
                                                                     targetomeDrugToMCLCluster.keySet());
            System.out.println("Total MCL drugs: " + targetomeDrugToMCLCluster.size());
            System.out.println("Total cancerx drugs: " + cancerrxDrugToCluster.size());
            System.out.println("\tAfter mapping to targetome drugs: " + targetomeDrugToCancerrxCluster.size());
            System.out.println("Total shared drugs: " + sharedDrugs.size());
//            sharedDrugs.forEach(drug -> System.out.println(drug + "\t" + drugNameToATCCode.get(drug)));

            targetomeDrugToMCLCluster.keySet().retainAll(sharedDrugs);
            targetomeDrugToCancerrxCluster.keySet().retainAll(sharedDrugs);

            double ari = calculateAdjustRandIndex(targetomeDrugToCancerrxCluster,
                                                  targetomeDrugToMCLCluster);
            System.out.println("ARI: " + ari);

            Map<Integer, Set<String>> cancerrxClusterToDrugs = new HashMap<>();
            targetomeDrugToCancerrxCluster.forEach((drug, cluster) -> {
                cancerrxClusterToDrugs.compute(cluster, (key, set) -> {
                    if (set == null)
                        set = new HashSet<>();
                    set.add(drug);
                    return set;
                });
            });
            Map<Integer, Set<String>> mclClusterToDrugs = new HashMap<>();
            targetomeDrugToMCLCluster.forEach((drug, cluster) -> {
                mclClusterToDrugs.compute(cluster, (key, set) -> {
                    if (set == null)
                        set = new HashSet<>();
                    set.add(drug);
                    return set;
                });
            });
            StringBuilder builder = new StringBuilder();
            builder.append("cancerrx_cluster");
            List<Integer> mclClusterList = mclClusterToDrugs.keySet().stream().sorted().collect(Collectors.toList());
            mclClusterList.forEach(cluster -> builder.append("\t").append(cluster));
            builder.append("\n");
            for (Integer cancerrxCluster : cancerrxClusterToDrugs.keySet()) {
                builder.append(cancerrxCluster);
                Set<String> rxDrugs = cancerrxClusterToDrugs.get(cancerrxCluster);
                for (Integer mclCluster : mclClusterList) {
                    Set<String> mclDrugs = mclClusterToDrugs.get(mclCluster);
                    Set<String> shared = InteractionUtilities.getShared(rxDrugs, mclDrugs);
                    builder.append("\t").append(shared.size());
                }
                builder.append("\n");
            }
            System.out.println(builder.toString());
        }
    }

    private Map<String, Integer> mapCancerRxClustersToTargetomeClusters(Map<String, Integer> cancerrxDrugToCluster) throws IOException {
        // Now we can compare
        Map<String, Integer> targetomeDrugInCancerrxToCluster = new HashMap<>();
        Map<String, Integer> targetomeDrugToCancerrxId = new HashMap<>();
        Map<String, String> cancerrxDrugToTargetomeDrug = loadCancerRxToTargetomeMap();
        cancerrxDrugToTargetomeDrug.forEach((key, value) -> {
            Integer cluster = cancerrxDrugToCluster.get(key);
            if (cluster == null)
                return; // In case there is nothing there
            int drugId = new Integer(key.split("_")[1]);
            //            System.out.println(key + "\t" + value + "\t" + cluster);
            if (targetomeDrugToCancerrxId.containsKey(value)) {
                // We will choose the drug having smaller id
                if (drugId > targetomeDrugToCancerrxId.get(value))
                    return;
            }
            targetomeDrugInCancerrxToCluster.put(value, cluster);
            targetomeDrugToCancerrxId.put(value, drugId);
        });
        return targetomeDrugInCancerrxToCluster;
    }
    
    @Test
    public void compareCancerRxDrugsInClusters() throws Exception {
        // Get all drugs
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        Set<String> drugNames = drugs.stream().map(drug -> drug.getDrugName()).collect(Collectors.toSet());
        session.close();
        
        Map<String, String> cancerrxDrugToTargetomeDrug = loadCancerRxToTargetomeMap();
        
        StringBuilder builder = new StringBuilder();
        builder.append("TargetomeMCLFile\tCancerRxFile\tRxConsensusCluster\tDrugsFromCancerRx\tDrugsInClusters\tAdjustedRandIndex\n");
        String[] mclFiles = getMCLFiles();
//        for (int j = 2; j <= 16; j++) {
//////            String mclFile = DIR + "ConsensusClustering/Targetome_100218/Targetome_100218.k=" + j + ".consensusClass.csv";
//            String mclFile = DIR + "targetome_target_111518/targetome_target_111518.k=" + j + ".consensusClass.csv";
////            String mclFile = DIR + "targetome_111518/targetome_111518.k=" + j + ".consensusClass.csv";
//////            System.out.println(mclFile);
//            Map<String, Integer> targetomeDrugInMCLClusters = loadDrugToConsensusCluster(mclFile);
        
//        String cancerrxTitle = "cancerrx_hc_031119_2/cancerrx_hc_031119_2.k=";
        String cancerrxTitle = "cancerrx_hc_kinase_inhibitors_032119/cancerrx_hc_kinase_inhibitors_032119.k=";
        
        for (String mclFile : mclFiles) {
            Map<String, Integer> targetomeDrugInMCLClusters = loadMCLDrugClusters(DIR + mclFile,
                                                                                  drugNames);
            for (int i = 2; i <= 16; i++) {
                //                String fileName = DIR + "cancerrx_km_111518/cancerrx_km_111518.k=" + i + ".consensusClass.csv";
                //                String fileName = DIR + "cancerrx_pam_111518/cancerrx_pam_111518.k=" + i + ".consensusClass.csv";
                //                
                //                String fileName = DIR + "cancerrx/cancerrx.k=" + i + ".consensusClass.csv";
                // Another run
                //                String fileName = DIR + "cancerrx_102618/cancerrx_102618.k=" + i + ".consensusClass.csv";
                // Another run
                //String fileName = DIR + "cancerrx_102618_2nd/cancerrx_102618_2nd.k=" + i + ".consensusClass.csv";
                // Repeated cancerx consensus clustering with reps = 10,000 (1,000) before to hope to
                // get consistent results
                String fileName = DIR + cancerrxTitle + i + ".consensusClass.csv";

                _compareCancerRxDrugsInClusters(cancerrxDrugToTargetomeDrug,
                                                targetomeDrugInMCLClusters,
                                                mclFile,
                                                i,
                                                fileName,
                                                builder);
            }
        }
        
        if (true) {
            System.out.println(builder.toString());
            return;
        }
        
        String[] iGraphFiles = {
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_infomap_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_infomap_undirected_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_infomap_unweighted_undirected_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_infomap_unweighted_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_louvain_undirected_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_louvain_unweighted_undirected_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_betweeness_unweighted_undirected_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_walktrap_031219.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_walktrap_unweighted_031219.txt"
        };
        for (String iGraphFile : iGraphFiles) {
            Map<String, Integer> drugsInIGraphClusters = loadiGraphClusters(DIR + iGraphFile,
                                                                            drugNames);
            for (int i = 2; i <= 16; i++) {
                //                String fileName = DIR + "cancerrx_km_111518/cancerrx_km_111518.k=" + i + ".consensusClass.csv";
                //                String fileName = DIR + "cancerrx_pam_111518/cancerrx_pam_111518.k=" + i + ".consensusClass.csv";
                //                
                //                String fileName = DIR + "cancerrx/cancerrx.k=" + i + ".consensusClass.csv";
                // Another run
                //                String fileName = DIR + "cancerrx_102618/cancerrx_102618.k=" + i + ".consensusClass.csv";
                // Another run
                //String fileName = DIR + "cancerrx_102618_2nd/cancerrx_102618_2nd.k=" + i + ".consensusClass.csv";
                // Repeated cancerx consensus clustering with reps = 10,000 (1,000) before to hope to
                // get consistent results
                String fileName = DIR + cancerrxTitle + i + ".consensusClass.csv";

                _compareCancerRxDrugsInClusters(cancerrxDrugToTargetomeDrug,
                                                drugsInIGraphClusters,
                                                iGraphFile,
                                                i,
                                                fileName,
                                                builder);
            }
        }
        
        System.out.println(builder.toString());
    }

    private void _compareCancerRxDrugsInClusters(Map<String, String> cancerrxDrugToTargetomeDrug,
                                                 Map<String, Integer> networkDrugClusters,
                                                 String networkClusterFileName,
                                                 int clusterNumber,
                                                 String cancerrxFileName,
                                                 StringBuilder builder) throws IOException {
        Map<String, Integer> cancerrxDrugToCluster = loadDrugToConsensusCluster(cancerrxFileName);
        //            cancerrxDrugToCluster.forEach((key, value) -> System.out.println(key + "\t" + value));

//                String fileName = DIR + "CancerRxSampleDrugNetwork_111518_mcl_I_8_out.txt";
//                Map<String, Integer> cancerrxDrugToCluster = loadMCLDrugClusters(fileName, cancerrxDrugToTargetomeDrug.keySet());
        
        Map<String, Integer> targetomeDrugInCancerrxClusters = mapCancerRxClustersToTargetomeClusters(cancerrxDrugToCluster);

        double adjustRandIndex = calculateAdjustRandIndex(networkDrugClusters, targetomeDrugInCancerrxClusters);
        int index = cancerrxFileName.lastIndexOf("/");
        cancerrxFileName = cancerrxFileName.substring(index + 1);
        index = networkClusterFileName.lastIndexOf("/");
        if (index >= 0)
            networkClusterFileName = networkClusterFileName.substring(index + 1);
        builder.append(networkClusterFileName + "\t" + 
                       cancerrxFileName + "\t" + 
                       clusterNumber + "\t" + 
                       targetomeDrugInCancerrxClusters.size() + "\t" + 
                       networkDrugClusters.size() + "\t" + 
                       adjustRandIndex + "\n");
    }
    
    /**
     * Load the mapping file generated by method mapTargetomeDrugs in class CancerRxDataAnalyzer
     * from the Druggability project. The output has been edited a little bit for mapping drugs
     * ended with "rescreen".
     * @return
     * @throws IOException
     */
    private Map<String, String> loadCancerRxToTargetomeMap() throws IOException {
        String fileName = "CancerRxDrugToTargetomeMap_102518.txt";
        try (Stream<String> stream = Files.lines(Paths.get(DIR, fileName))) {
            Map<String, String> cancerrxToDrug = new HashMap<>();
            stream.skip(1)
                  .forEach(line -> {
                      String[] tokens = line.split("\t");
                      if (tokens.length == 4)
                          cancerrxToDrug.put("Drug_" + tokens[0] + "_IC50", tokens[3]);
                  });
            return cancerrxToDrug;
        }
    }
    
    @Test
    public void calculateAdjustedRandIndexForATCAndClusters() throws Exception {
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        System.out.println("Total drugs in targetome: " + drugs.size());
        // Since a drug may have multiple ATC codes, targetome picked up the first code currently.
        // Collect drugs with codes starting with L01 for cancer drugs directly to minimize this issue.
        drugs = drugs.stream().filter(drug -> drug.getAtcClassID().startsWith("L01")).collect(Collectors.toList());
        System.out.println("Filter to L01 drugs: " + drugs.size());
        
        Map<String, String> drugToATCClass = new HashMap<>();
        drugs.forEach(drug -> drugToATCClass.put(drug.getDrugName(), drug.getAtcClassName()));
        List<String> atcClassList = drugToATCClass.values().stream().distinct().sorted().collect(Collectors.toList());
        Map<String, Integer> drugToATCClassIndex = new HashMap<>();
        for (String drug : drugToATCClass.keySet()) {
            String atcClass = drugToATCClass.get(drug);
            int index = atcClassList.indexOf(atcClass);
            drugToATCClassIndex.put(drug, index);
        }
        session.close();
        
        if (true) {
            Map<String, String> atcCodeToName = new HashMap<>();
            drugs.stream().forEach(drug -> atcCodeToName.put(drug.getAtcClassID(), drug.getAtcClassName()));
            atcCodeToName.keySet().stream().sorted().forEach(code -> System.out.println(code + "\t" + atcCodeToName.get(code)));
            return;
        }
        
        System.out.println("\nAdjusted Rand Index between Drug Clusters and ATC classes:");
        StringBuilder builder = new StringBuilder();
        
        // Consensus clustering
        for (int i = 2; i < 16; i++) {
//            String fileName = DIR + "/ConsensusClustering/Targetome_100218/Targetome_100218.k=" + i + ".consensusClass.csv";
//            String fileName = DIR + "targetome_111518/targetome_111518.k=" + i + ".consensusClass.csv";
            String fileName = DIR + "targetome_011919_030419/targetome_011919_030419.k=" + i + ".consensusClass.csv";
            Map<String, Integer> drugToCluster = loadDrugToConsensusCluster(fileName);
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        // Another consensus clustering
        for (int i = 2; i < 16; i++) {
            String fileName = DIR + "targetome_011919_1_030419/targetome_011919_1_030419.k=" + i + ".consensusClass.csv";
            Map<String, Integer> drugToCluster = loadDrugToConsensusCluster(fileName);
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        // Yet another consensus clustering
        for (int i = 2; i < 16; i++) {
            String fileName = DIR + "targetome_012119_030419/targetome_012119_030419.k=" + i + ".consensusClass.csv";
            Map<String, Integer> drugToCluster = loadDrugToConsensusCluster(fileName);
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        for (int i = 2; i < 16; i++) {
//          String fileName = DIR + "/ConsensusClustering/Targetome_100218/Targetome_100218.k=" + i + ".consensusClass.csv";
          String fileName = DIR + "targetome_target_111518/targetome_target_111518.k=" + i + ".consensusClass.csv";
          Map<String, Integer> drugToCluster = loadDrugToConsensusCluster(fileName);
          calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        Map<String, String> cancerrxDrugToTargetomeDrug = loadCancerRxToTargetomeMap();
        for (int i = 2; i <= 16; i++) {
//            String fileName = DIR + "cancerrx_km_111518/cancerrx_km_111518.k=" + i + ".consensusClass.csv";
//            String fileName = DIR + "cancerrx_102618_2nd/cancerrx_102618_2nd.k=" + i + ".consensusClass.csv";
//            // Another run
            String fileName = DIR + "cancerrx_102618/cancerrx_102618.k=" + i + ".consensusClass.csv";
            Map<String, Integer> cancerrxDrugToCluster = loadDrugToConsensusCluster(fileName);
            Map<String, Integer> drugToCluster = new HashMap<>();
            cancerrxDrugToCluster.forEach((rxDrug, cluster) -> {
                String ttDrug = cancerrxDrugToTargetomeDrug.get(rxDrug);
                if (ttDrug == null)
                    return;
                drugToCluster.put(ttDrug, cluster);
            });
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        String[] fileNames = new String[] {
                "CancerRxSampleDrugNetwork_111518_mcl_I_2_out.txt",
                "CancerRxSampleDrugNetwork_111518_mcl_I_4_out.txt",
                "CancerRxSampleDrugNetwork_111518_mcl_I_6_out.txt",
                "CancerRxSampleDrugNetwork_111518_mcl_I_8_out.txt",
        };
        for (String fileName : fileNames) {
            Map<String, Integer> cancerrxDrugToCluster = loadMCLDrugClusters(DIR + fileName, cancerrxDrugToTargetomeDrug.keySet());
            Map<String, Integer> drugToCluster = new HashMap<>();
            cancerrxDrugToCluster.forEach((rxDrug, cluster) -> {
                String ttDrug = cancerrxDrugToTargetomeDrug.get(rxDrug);
                if (ttDrug == null)
                    return;
                drugToCluster.put(ttDrug, cluster);
            });
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        // Two runs of consensus clustering
        String[] bipartiteFiles = new String[] {
                "Targetome_Drug_Pathway_Bipartite_Clusters_101218.txt",
                "Targetome_Drug_Pathway_Bipartite_Clusters_2nd_Run_101218.txt"
        };
        for (String fileName : bipartiteFiles) {
            fileName = DIR + fileName;
            Map<String, Integer> drugToCluster = loadBipartiteClusters(fileName, drugToATCClass.keySet());
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        String[] mclFiles = getMCLFiles();
        for (String mclFile : mclFiles) {
            String fileName = DIR + mclFile;
            Map<String, Integer> drugToCluster = loadMCLDrugClusters(fileName, drugToATCClass.keySet());
            calculateAdjustedRandIndexForATCAndClusters(drugToATCClassIndex, drugToCluster, fileName, builder);
        }
        
        System.out.println("FileName\tDrugsInClusters\tAdjustedRandIndex");
        System.out.println(builder.toString());
    }

    private String[] getMCLFiles() {
        // mcl based clustering
        String[] mclFiles = new String[] {
                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032119_mcl_out.txt",
                "Targetome_Drug_Target_102518_L01XE_032119_032119_mcl_out.txt",
                "Targetome_Drug_Target_Weight_102518_L01XE_032119_032119_mcl_out.txt"
                
                // Repeated results after the code improvement regarding cycle handling
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_mcl_out.txt",
//                // Yet another impact analysis
//                "Targetome_Drug_Pathway_Interaction_012119_012119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_012119_no_weight_012119_mcl_out.txt",
////                "Targetome_Drug_Pathway_Interaction_01_012119_012119_mcl_out.txt",
//                // Yet another impact analysis with force_connected = false
//                "Targetome_Drug_Pathway_Interaction_012119_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_012119_no_weight_012119_no_force_mcl_out.txt",
////                "Targetome_Drug_Pathway_Interaction_01_012119_012119_no_force_mcl_out.txt",
//                // Another impact analysis
//                "Targetome_Drug_Pathway_Interaction_011919_1_012119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_1_no_weight_012119_mcl_out.txt",
////                "Targetome_Drug_Pathway_Interaction_01_011919_1_012119_mcl_out.txt",
//                // No force_connected = true
//                "Targetome_Drug_Pathway_Interaction_011919_1_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_1_no_weight_012119_no_force_mcl_out.txt",
////                "Targetome_Drug_Pathway_Interaction_01_011919_1_012119_no_force_mcl_out.txt",
//                // Used a new drug pathway impact file to repeat experiments
//                "Targetome_Drug_Pathway_Interaction_011919_011919_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_no_weight_011919_mcl_out.txt",
////                "Targetome_Drug_Pathway_Interaction_01_011919_011919_mcl_out.txt",
//                // No force_connected=y
//                "Targetome_Drug_Pathway_Interaction_011919_011919_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_no_weight_011919_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_011919_no_force_mcl_out.txt",
                // Repeat mcl clustering for drug/pathways impact on January 18, 2019
//                "Targetome_Drug_Pathway_Interaction_100218_011819_mcl_I_2_out.txt",
//                "Targetome_Drug_Pathway_Interaction_100218_no_weight_011819_mcl_I_2_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_100218_011819_mcl_I_2_out.txt",
//                // Repeat mcl clustering for drug/pathways impact on January 18, 2019
//                // But without forcing connected
//                "Targetome_Drug_Pathway_Interaction_100218_011819_mcl_I_2_no_force_out.txt",
//                "Targetome_Drug_Pathway_Interaction_100218_no_weight_011819_mcl_I_2_no_force_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_100218_011819_mcl_I_2_no_force_out.txt",
//                // Drug pathway impact network
//                "Targetome_Drug_Pathway_Interaction_01_100218_mcl_I_2_out_102418.txt",
//                // As above but without threshold = 0.01
////                "Targetome_Drug_Pathway_Interaction_100118_mcl_I_2_out_102518.txt",
//                "Targetome_Drug_Pathway_Interaction_100218_011519_mcl_I_2_out_011519.txt",
//                // As the first one, but without weights
//                "Targetome_Drug_Pathway_Interaction_01_100218_no_weight_mcl_I_2_out_102518.txt",
//                // As the second one, but without weights
////                "Targetome_Drug_Pathway_Interaction_100118_no_weight_mcl_I_2_out_102518.txt",
//                "Targetome_Drug_Pathway_Interaction_100218_011519_no_weight_mcl_I_2_out_011519.txt",
                // Drug target network without edges
//                "Targetome_Drug_Target_102518_mcl_I_2_out.txt",
                // Drug target network with edges
//                "Targetome_Drug_Target_Weight_102518_mcl_I_2_out.txt"
//                "Targetome_Drug_Target_102518_030719_mcl_out.txt",
//                "Targetome_Drug_Target_Weight_102518_030719_mcl_out.txt",
//                "Targetome_Drug_Target_102518_030719_no_force_mcl_out.txt",
//                "Targetome_Drug_Target_Weight_102518_030719_no_force_mcl_out.txt"
        };
        return mclFiles;
    }

    private void calculateAdjustedRandIndexForATCAndClusters(Map<String, Integer> drugToATCClassIndex,
            Map<String, Integer> drugToCluster, String fileName, StringBuilder builder) {
        double value = calculateAdjustRandIndex(drugToCluster, drugToATCClassIndex);
        int index = fileName.lastIndexOf("/");
        builder.append(fileName.substring(index + 1) + "\t" + 
                       drugToCluster.size() + "\t" + 
                       value + "\n");
    }
    
    /**
     * It is found that two MCL clustering results having diffrent drugs. This method 
     * is designed to see why.
     * @throws Exception
     */
    @Test
    public void checkDrugsInMCLClusters() throws Exception {
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        System.out.println("Total drugs in targetome: " + drugs.size());
        Map<String, String> drugToATCClass = new HashMap<>();
        drugs.forEach(drug -> drugToATCClass.put(drug.getDrugName(), drug.getAtcClassName()));
        List<String> atcClassList = drugToATCClass.values().stream().distinct().sorted().collect(Collectors.toList());
        Map<String, Integer> drugToATCClassIndex = new HashMap<>();
        for (String drug : drugToATCClass.keySet()) {
            String atcClass = drugToATCClass.get(drug);
            int index = atcClassList.indexOf(atcClass);
            drugToATCClassIndex.put(drug, index);
        }
        session.close();
        
        String[] mclFiles = getMCLFiles();
        int counter = 0;
        Set<String> drugs1 = null;
        Set<String> drugs2 = null;
        for (String mclFile : mclFiles) {
            String fileName = DIR + mclFile;
            Map<String, Integer> drugToCluster = loadMCLDrugClusters(fileName, drugToATCClass.keySet());
            System.out.println(mclFile + ": " + drugToCluster.size());
            counter ++;
            if (counter == 1)
                drugs1 = drugToCluster.keySet();
            else if (counter == 2)
                drugs2 = drugToCluster.keySet();
//            if (counter == 2)
//                break;
        }
        Set<String> shared = InteractionUtilities.getShared(drugs1, drugs2);
        System.out.println("Shared: " + shared.size());
        drugs1.removeAll(shared);
        System.out.println("Unique in first: " + drugs1);
        drugs2.removeAll(shared);
        System.out.println("Unique in second: " + drugs2);
    }
    
    @Test
    public void checkDrugClusters() throws Exception {
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        System.out.println("Total drugs: " + drugs.size());
        Map<String, Drug> nameToDrug = new HashMap<>();
        drugs.forEach(drug -> nameToDrug.put(drug.getDrugName(), drug));
        
//        String fileName = DIR + "/ConsensusClustering/Targetome_100218/Targetome_100218.k=7.consensusClass.csv";
//        Map<String, Integer> drugToCluster = loadDrugToConsensusCluster(fileName);
        
//        String fileName = DIR + "Targetome_Drug_Pathway_Bipartite_Clusters_101218.txt";
//        Map<String, Integer> drugToCluster = loadBipartiteClusters(fileName, nameToDrug.keySet());
        
        String fileName = DIR + "Targetome_Drug_Pathway_Interaction_01_100218_mcl_I_2_out_102418.txt";
        fileName = DIR + "Targetome_Drug_Pathway_Interaction_100118_mcl_I_2_out_102518.txt";
        
        // Clustering based on drug/target network
//        fileName = DIR + "Targetome_Drug_Target_102518_mcl_I_2_out.txt";
//        fileName = DIR + "Targetome_Drug_Target_Weight_102518_mcl_I_2_out.txt";
        
        Map<String, Integer> drugToCluster = loadMCLDrugClusters(fileName, nameToDrug.keySet());
        
        System.out.println("Total drug to cluster: " + drugToCluster.size());
        
        System.out.println("Drug\tCluster\tATCClass\tEPCClass");
        Map<String, String> drugToATCClass = new HashMap<>();
        drugToCluster.forEach((key, cluster) -> {
            Drug drug = nameToDrug.get(key);
            System.out.println(key + "\t" + 
                               cluster + "\t" + 
                               drug.getAtcClassName() + "\t" + 
                               drug.getEpcClassName());
            drugToATCClass.put(key, drug.getAtcClassName());
        });
        session.close();
        List<String> atcClassList = drugToATCClass.values().stream().distinct().sorted().collect(Collectors.toList());
        Map<String, Integer> drugToATCClassIndex = new HashMap<>();
        for (String drug : drugToATCClass.keySet()) {
            String atcClass = drugToATCClass.get(drug);
            int index = atcClassList.indexOf(atcClass);
            drugToATCClassIndex.put(drug, index);
        }
        
        double value = calculateAdjustRandIndex(drugToCluster, drugToATCClassIndex);
        System.out.println("\nAdjust Rand Index between drug MCL clusters and ATC Classes: " + value);
    }
    
    private double calculateAdjustRandIndex(Map<String, Integer> drugToCluster1,
                                            Map<String, Integer> drugToCluster2) {
        Set<String> shared = InteractionUtilities.getShared(drugToCluster1.keySet(), 
                                                            drugToCluster2.keySet());
        int[] cluster1 = new int[shared.size()];
        int[] cluster2 = new int[shared.size()];
        int index = 0;
        for (String drug : shared) {
            cluster1[index] = drugToCluster1.get(drug);
            cluster2[index] = drugToCluster2.get(drug);
            index ++;
        }
        
//        StringBuilder builder = new StringBuilder();
//        for (int i : cluster1)
//            builder.append(i).append(", ");
//        System.out.println(builder.toString());
//        builder.setLength(0);
//        for (int i : cluster2)
//            builder.append(i).append(", ");
//        System.out.println(builder.toString());
        
        AdjustedRandIndex adjustRandIndex = new AdjustedRandIndex();
        double value = adjustRandIndex.measure(cluster1, cluster2);
        return value;
    }
    
    /**
     * Load the clusters for drugs based on clustering analysis results in R.
     * @param fileName
     * @return
     * @throws IOException
     */
    private Map<String, Integer> loadDrugToConsensusCluster(String fileName) throws IOException {
        Map<String, Integer> drugToCluster = new HashMap<>();
        try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
            stream.forEach(line -> {
                String[] tokens = line.split(",");
                String drug = tokens[0].substring(1, tokens[0].length() - 1);
                drugToCluster.put(drug, new Integer(tokens[1]));
            });
        }
        return drugToCluster;
    }
    
    private Map<String, Integer> loadiGraphClusters(String fileName,
                                                    Set<String> drugs) throws IOException {
        Map<String, Integer> nodeToCluster = new HashMap<>();
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            lines.map(line -> line.split("\t"))
                 .forEach(tokens -> {
                     nodeToCluster.put(tokens[0],
                                       Integer.parseInt(tokens[1]));
                 });
        }
        nodeToCluster.keySet().retainAll(drugs);
        return nodeToCluster;
    }
    
    private Map<String, Integer> loadBipartiteClusters(String fileName, Set<String> drugs) throws IOException {
        Map<String, Integer> nodeToCluster = new HashMap<>();
        fu.setInput(fileName);
        String line = null;
        int cluster = -1;
        while ((line = fu.readLine()) != null) {
            line = line.trim();
            if (line.length() == 0)
                continue;
            if (line.startsWith("Rownames") || line.startsWith("Colnames") ||
                line.startsWith("Depth") || line.startsWith("__"))
                continue;
            if (line.startsWith("Nr of module:")) {
                cluster = Integer.parseInt(line.split(":")[1].trim());
                continue;
            }
            // Come to the node part
            if (drugs != null && !drugs.contains(line))
                continue;
            nodeToCluster.put(line, cluster);
        }
        fu.close();
        return nodeToCluster;
    }
    
    @Test
    public void convertBipartiteClusterToAttFile() throws IOException {
        String fileName = DIR + "Targetome_Drug_Pathway_Bipartite_Clusters_101218.txt";
        String outFileName = DIR + "Targetome_Drug_Pathway_Bipartite_Clusters_101218.att";
        Map<String, Integer> nodeToModule = loadBipartiteClusters(fileName, null);
        nodeToModule.forEach((node, module) -> System.out.println(node + "\t" + module));
        fu.setOutput(outFileName);
        fu.printLine("Node\tBipartiteModule");
        for (String node : nodeToModule.keySet())
            fu.printLine(node + "\t" + nodeToModule.get(node));
        fu.close();
    }
    
    /**
     * Load drugs from the database.
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    public List<Drug> loadDrugs(Session session) throws Exception {
        Query query = session.createQuery("FROM Drug as drug");
        List<Drug> drugs = query.list();
        return drugs;
    }
    
    /**
     * Filter drugs that don't have any impact or impact less than a pre-defined cutoff value.
     * @param drugToResults
     */
    private void filterImpactResults(Map<String, List<DrugImpactResult>> drugToResults) {
        Set<String> drugs = drugToResults.keySet();
        for (Iterator<String> it = drugs.iterator(); it.hasNext();) {
            String drug = it.next();
            List<DrugImpactResult> results = drugToResults.get(drug);
            for (Iterator<DrugImpactResult> resultIt = results.iterator(); resultIt.hasNext();) {
                DrugImpactResult result = resultIt.next();
                if (result.value < IMPACT_CUTOFF)
                    resultIt.remove();
            }
            if (results.isEmpty())
                it.remove();
        }
    }
    
    @Test
    public void generateDrugPathwayMatrix() throws IOException {
        // For targetome
//        String fileName = "Targetome_Impact_091918.txt";
        String fileName = "Targetome_Impact_100218.txt";
////        String outFileName = "Targetome_Impact_Matrix_091918.txt";
//        String outFileName = "Targetome_Impact_Matrix_Filtered_091918.txt";
//        String outFileName = "Targetome_Impact_Matrix_Filtered_01_091918.txt";
//        String outFileName = "Targetome_Impact_Matrix_Filtered_05_091918.txt";
        String outFileName = "Targetome_Impact_Matrix_Filtered_01_100218.txt";
        boolean needFilter = true;
        
        // Repeated impact analysis with bug fix
        fileName = "Targetome_Impact_011919.txt";
        outFileName = "Targetome_Impact_011919_Matrix_030419.txt";
        needFilter = false;
        
        // Repeat
        fileName = "Targetome_Impact_011919_1.txt";
        outFileName = "Targetome_Impact_011919_1_Matrix_030419.txt";
        needFilter = false;
        
        // Another repeat
        fileName = "Targetome_Impact_012119.txt";
        outFileName = "Targetome_Impact_012119_Matrix_030419.txt";
        needFilter = false;
        
        // For DrugCentral
//        fileName = "DrugCentral_Impact_091918.txt";
//        String outFileName = "DrugCentral_Impact_Matrix_091918.txt";
//        String outFileName = "DrugCentral_Impact_Matrix_Filtered_091918.txt";
//        String outFileName = "DrugCentral_Impact_Matrix_Filtered_01_091918.txt";
//        outFileName = "DrugCentral_Impact_Matrix_Filtered_05_091918.txt";
        Map<String, List<DrugImpactResult>> drugToResults = loadDrugToResults(DIR + fileName);
        generateOutput(drugToResults, needFilter, DIR + outFileName);
    }
    
    @Test
    public void generateDrugTargetMatrix() throws IOException {
        String fileName = "Targetome_Drug_Target_Weight_102518.txt";
        String outFileName = "Targetome_Drug_Target_Weight_102518_Matrix_111518.txt";
        Map<String, List<DrugImpactResult>> drugToResults = loadDrugToTarget(DIR + fileName);
        generateOutput(drugToResults, false, DIR + outFileName);
    }
    
    @Test
    public void generateDrugTargetNetwork() throws Exception {
        Session session = getSession();
        Query query = session.createQuery("FROM " + Interaction.class.getName());
        List<Interaction> interactions = query.list();
        System.out.println("Total interactions: " + interactions.size());

        //        String fileName = DIR + "Targetome_Drug_Target_102518.txt";
        //        fu.setOutput(fileName);
        //        for (Interaction i : interactions) {
        //            fu.printLine(i.getIntDrug().getDrugName() + "\t" + i.getIntTarget().getTargetName());
        //        }
        //        fu.close();

        String fileName = DIR + "Targetome_Drug_Target_Weight_102518.txt";
        fu.setOutput(fileName);
        DefaultAffinityToModificationMap mapper = new DefaultAffinityToModificationMap();
        for (Interaction i : interactions) {
            Double minValue = getMinValue(i);
            if (minValue == null)
                continue;
            double weight = mapper.getModificationStrenth(minValue);
            // The affinity is too small to be considered
            if (weight <= 0.0d)
                continue;
            fu.printLine(i.getIntDrug().getDrugName() + "\t" + i.getIntTarget().getTargetName() + "\t" + weight);
        }
        fu.close();

        session.close();
    }
    
    // This method is copied from the caBigR3Web project. There is a bug in the used code ExpEvidence
    // that has not considered median value is "". However, the code in the Druggabilty project has 
    // fixed this problem.
    /**
     * Use this method to pick up whatever type of minimum value
     * @param interaction
     * @return
     */
    private Double getMinValue(Interaction interaction) {
        Double rtn = null;
        if (interaction.getExpEvidenceSet() == null)
            return rtn;
        for (ExpEvidence evidence : interaction.getExpEvidenceSet()) {
            if (evidence.getAssayValueMedian() == null ||
                evidence.getAssayValueMedian().trim().length() == 0 ||
                evidence.getAssayType() == null)
                continue;
            Number current = evidence.getAssayValue();
            if (current == null)
                continue;
            if (rtn == null || current.doubleValue() < rtn.doubleValue())
                rtn = current.doubleValue();
        }
        return rtn;
    }
    
    @Test
    public void generateCancerRxDrugCellLinesNetwork() throws IOException {
        String fileName = "datasets/CancerRxGene/IC50_v17.csv";
        String outFileName = DIR + "CancerRxSampleDrugNetwork_111518.txt";
        fu.setInput(fileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        String[] headers = line.split(",");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            for (int i = 1 ; i < tokens.length; i++) {
                if (tokens[i].length() == 0)
                    continue;
                Double value = new Double(tokens[i]); 
                value = Math.pow(2.0d, -value);
                fu.printLine(tokens[0] + "\t" + headers[i] + "\t" + value);
            }
        }
        fu.close();
    }
    
    @Test
    public void performMCLClustering() throws Exception {
        String[] inFileNames = new String[] {
                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032119.txt",
                "Targetome_Drug_Target_102518_L01XE_032119.txt",
                "Targetome_Drug_Target_Weight_102518_L01XE_032119.txt"
                
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032019.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032019.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032019.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032019.txt"
//                "Targetome_Drug_Target_Weight_102518.txt",
//                "Targetome_Drug_Target_102518.txt"
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119.txt"
//                "Targetome_Drug_Pathway_Interaction_011919_011919.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_no_weight_011919.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_011919.txt"
//                "Targetome_Drug_Pathway_Interaction_011919_1_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_1_no_weight_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_1_012119.txt"
//                "Targetome_Drug_Pathway_Interaction_012119_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_012119_no_weight_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_01_012119_012119.txt"
        };
        String[] outFileNames = new String[] {
              "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032119_mcl_out.txt",
              "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032119_mcl_out.txt",
              "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032119_mcl_out.txt",
              "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032119_mcl_out.txt",
              "Targetome_Drug_Target_102518_L01XE_032119_032119_mcl_out.txt",
              "Targetome_Drug_Target_Weight_102518_L01XE_032119_032119_mcl_out.txt"
              
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_L01XE_032019_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_L01XE_032019_mcl_out.txt"
//                "Targetome_Drug_Target_Weight_102518_030819_mcl_out.txt",
//                "Targetome_Drug_Target_102518_030819_mcl_out.txt"
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_no_force_mcl_out.txt"
//                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_1_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_2_0307119_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_030619_3_0307119_mcl_out.txt"
//                "Targetome_Drug_Pathway_Interaction_011919_011919_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_no_weight_011919_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_011919_no_force_mcl_out.txt"
//                "Targetome_Drug_Pathway_Interaction_011919_1_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_1_no_weight_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_1_012119_no_force_mcl_out.txt"
//                "Targetome_Drug_Pathway_Interaction_012119_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_012119_no_weight_012119_no_force_mcl_out.txt",
//                "Targetome_Drug_Pathway_Interaction_01_012119_012119_no_force_mcl_out.txt"
        };
        String script = "resources/mcl_script.sh";
//        String script = "resources/mcl_script_no_force.sh";
        for (int i = 0; i < inFileNames.length; i++) {
            System.out.println(inFileNames[i]);
            ProcessRunner runner = new ProcessRunner();
            List<String> parameters = new ArrayList<String>();
            // Have to use a bash script. Make sure bash is in the server-side
            parameters.add("bash");
            parameters.add(script);
            // The following parameters are passed into the shell script.
            parameters.add(DIR + inFileNames[i]);
            parameters.add("2.0");
            parameters.add(DIR + outFileNames[i]);
            String[] output = runner.runScript(parameters.toArray(new String[0]));
            System.out.println("System.out: " + output[0]);
            System.out.println("System.err: " + output[1]);
        }
//        checkDrugClusterOverlaps();
    }
    
    @Test
    public void filterCancerrxToKinaseInhibitors() throws Exception {
//        Set<String> inhibitors = loadKinaseInhibitors();
        Set<String> inhibitors = loadSharedL01XEDrugs();
        System.out.println("Total kinase inhibitors: " + inhibitors.size());
        
        Map<String, String> cancerrxIdToDrug = loadCancerRxToTargetomeMap();
        // Perform a filter
        for (Iterator<String> it = cancerrxIdToDrug.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            String drug = cancerrxIdToDrug.get(id);
            if (!inhibitors.contains(drug))
                it.remove();
        }
        System.out.println("Total selected cancerrx ids: " + cancerrxIdToDrug.size());
        
        String inFile = "datasets/CancerRxGene/IC50_v17.csv";
        String outFile = DIR + "cancerrx_ic50_v17_kinase_inhibitors.csv";
        
        fu.setInput(inFile);
        fu.setOutput(outFile);
        
        String line = fu.readLine();
        // Get a list of drugs we want
        Set<Integer> neededCols = new HashSet<>();
        String[] headers = line.split(",");
        StringBuilder builder = new StringBuilder();
        builder.append(headers[0]);
        for (int i = 1; i < headers.length; i++) {
            if (cancerrxIdToDrug.keySet().contains(headers[i])) {
                neededCols.add(i);
                builder.append(",").append(headers[i]);
            }
        }
        fu.printLine(builder.toString());
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            builder.setLength(0);
            builder.append(tokens[0]);
            for (int i = 1; i < tokens.length; i++) {
                if (!neededCols.contains(i))
                    continue;
                builder.append(",").append(tokens[i]);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * Load drugs annotated with ATC code L01XE from the cancer targetome.
     * @return
     * @throws Exception
     */
    private Set<String> loadKinaseInhibitors() throws Exception {
        // Get drugs with ATC code = L01XE
        Session session = getSession();
        List<Drug> drugs = loadDrugs(session);
        String atcCode = "L01XE";
        Set<String> drugNames = drugs.stream()
                .filter(drug -> drug.getAtcClassID().equals(atcCode))
                .map(drug -> drug.getDrugName())
                .collect(Collectors.toSet());
        session.close();
        System.out.println("Total L01XE drugs: " + drugNames.size());
        return drugNames;
    }
    
    /**
     * Perform this method to filter drug target network to drugs with ATC Code = L01XE,
     * which means Protein Kinase Inhibitor.
     * @throws IOException
     */
    @Test
    public void filterDrugPathwayNetworkToKinaseInhibitors() throws Exception {
//        Set<String> drugNames = loadKinaseInhibitors();
//        String filePostfix = "_L01XE_032019";
        Set<String> drugNames = loadSharedL01XEDrugs();
        String filePostfix = "_L01XE_032119";
        String[] inFiles = {
                "Targetome_Drug_Pathway_Interaction_030619_3_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_2_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119.txt",
                "Targetome_Drug_Target_102518.txt",
                "Targetome_Drug_Target_Weight_102518.txt"
        };
        
        String line = null;
        for (String inFile : inFiles) {
            String outFile = inFile.split("\\.")[0] + filePostfix + ".txt";
            fu.setInput(DIR + inFile);
            fu.setOutput(DIR + outFile);
            while ((line = fu.readLine()) != null) {
                if (drugNames.contains(line.split("\t")[0]))
                    fu.printLine(line);
            }
            fu.close();
        }
    }
    
    @Test
    public void generateDrugPathwayNetworkInBatch() throws IOException {
        String[] inFiles = {
                "Targetome_Impact_030619_3.txt",
                "Targetome_Impact_030619_2.txt",
                "Targetome_Impact_030619_1.txt",
                "Targetome_Hit_030619_1.txt"
        };
        String[] outFiles = {
                "Targetome_Drug_Pathway_Interaction_030619_3_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_2_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_0307119.txt",
                "Targetome_Drug_Pathway_Interaction_030619_1_Hit_0307119.txt"
        };
        boolean[] needWeights = {true, true, true, false};
        
        for (int i = 0; i < inFiles.length; i++) {
            Map<String, List<DrugImpactResult>> drugToResults = loadDrugToResults(DIR + inFiles[i]);
            fu.setOutput(DIR + outFiles[i]);
            for (String drug : drugToResults.keySet()) {
                List<DrugImpactResult> results = drugToResults.get(drug);
                for (DrugImpactResult result : results) {
                    if (!needWeights[i]) { // Some of pathways have targets mapped but not impacted 
                        // (e.g. inhibitors to one of activators working together on a reaction)
                        fu.printLine(drug + "\t" + result.pathway.trim());
                    }
                    else {
                        fu.printLine(drug + "\t" + result.pathway.trim() + "\t" + result.value);
                    }
                }
            }
            fu.close();
        }
    }
    
    @Test
    public void generateDrugPathwayNetwork() throws IOException {
        // For targetome
        String fileName = "Targetome_Impact_100218.txt";
        // There is a bug in the code generating this file
        fileName = "Targetome_Impact_011819.txt";
        // After the bug is fixed
        fileName = "Targetome_Impact_011919.txt";
        // Another run
        fileName = "Targetome_Impact_011919_1.txt";
        // yet another run
        fileName = "Targetome_Impact_012119.txt";

        boolean[] needWeights = new boolean[] {false, true, true};
        Double[] minValues = new Double[] {0.0d, 0.0d, 0.01d};
//        String[] outFileNames = new String[] {
//                "Targetome_Drug_Pathway_Interaction_011919_1_no_weight_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_011919_1_012119.txt",
//                "Targetome_Drug_Pathway_Interaction_01_011919_1_012119.txt"
//        };
        String[] outFileNames = new String[] {
                "Targetome_Drug_Pathway_Interaction_012119_no_weight_012119.txt",
                "Targetome_Drug_Pathway_Interaction_012119_012119.txt",
                "Targetome_Drug_Pathway_Interaction_01_012119_012119.txt"
        };
        for (int i = 0; i < needWeights.length; i++) {
            boolean needWeight = needWeights[i];
            Double minValue = minValues[i];
            String outFileName = outFileNames[i];
            
            //        String outFileName = "Targetome_Drug_Pathway_Interaction_100118.txt";
            //        String nodeAttFileName = "Targetome_Drug_Pathway_Interaction_100118.att";
            //        boolean needWeight = true;
            //        Double minValue = 0.01d;
            //        String outFileName = "Targetome_Drug_Pathway_Interaction_01_011919_011919.txt";
            //        String outFileName = "Targetome_Drug_Pathway_Interaction_100218_no_weight_011819_1.txt";
            //        String nodeAttFileName = "Targetome_Drug_Pathway_Interaction_100218_011519.att";

            //        Double minValue = 0.01d;
            //        String outFileName = "Targetome_Drug_Pathway_Interaction_01_100218.txt";
            //        String nodeAttFileName = "Targetome_Drug_Pathway_Interaction_01_100218.att";

            Map<String, String> nodeToType = new HashMap<>();
            // For DrugCentral
            //        fileName = "DrugCentral_Impact_091918.txt";
            //        String outFileName = "DrugCentral_Impact_Matrix_091918.txt";
            //        String outFileName = "DrugCentral_Impact_Matrix_Filtered_091918.txt";
            //        String outFileName = "DrugCentral_Impact_Matrix_Filtered_01_091918.txt";
            //        outFileName = "DrugCentral_Impact_Matrix_Filtered_05_091918.txt";
            Map<String, List<DrugImpactResult>> drugToResults = loadDrugToResults(DIR + fileName);
            fu.setOutput(DIR + outFileName);
            //        fu.printLine("Drug\tPathway\tImpact");
            for (String drug : drugToResults.keySet()) {
                List<DrugImpactResult> results = drugToResults.get(drug);
                for (DrugImpactResult result : results) {
                    if (!needWeight) { // Some of pathways have targets mapped but not impacted 
                        // (e.g. inhibitors to one of activators working together on a reaction)
                        fu.printLine(drug + "\t" + result.pathway.trim());
                    }
                    else if (result.value > minValue) {
                        fu.printLine(drug + "\t" + result.pathway.trim() + "\t" + result.value);
                    }
                    nodeToType.put(drug, "Drug");
                    nodeToType.put(result.pathway, "Pathway");
                }
            }
            fu.close();
        }
        
//        if (true)
//            return;
//        
//        fu.setOutput(DIR + nodeAttFileName);
//        fu.printLine("Node\tType");
//        for (String node : nodeToType.keySet())
//            fu.printLine(node + "\t" + nodeToType.get(node));
//        fu.close();
    }
    
    private Map<String, List<DrugImpactResult>> loadDrugToTarget(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, List<DrugImpactResult>> drugToResults = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            DrugImpactResult result = new DrugImpactResult();
            result.drug = tokens[0];
            result.pathway = tokens[1];
            result.value = new Double(tokens[2]);
            drugToResults.compute(result.drug, (key, list) -> {
                if (list == null)
                    list = new ArrayList<>();
                list.add(result);
                return list;
            });
        }
        fu.close();
        return drugToResults;
    }
    
    /**
     * Drugs in this shared list were manually collected by running method checkDrugClusterOverlaps().
     * @return
     * @throws IOException
     */
    private Set<String> loadSharedL01XEDrugs() throws IOException {
        String fileName = "CommonL01XEDrugs_032119.txt";
        try (Stream<String> stream = Files.lines(Paths.get(DIR, fileName))) {
            return stream.map(line -> line.split("\t")[0]).collect(Collectors.toSet());
        }
    }
    
    private Map<String, List<DrugImpactResult>> loadDrugToResults(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        int count = 0;
        Map<String, List<DrugImpactResult>> drugToResults = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 3) {
                System.out.println(line);
                count ++;
                continue;
            }
            DrugImpactResult result = new DrugImpactResult();
            result.drug = tokens[0];
            result.pathway = tokens[2];
            if (tokens.length > 4) // We may see drug pathway hit results
                result.value = new Double(tokens[4]);
            drugToResults.compute(result.drug, (key, list) -> {
                if (list == null)
                    list = new ArrayList<>();
                list.add(result);
                return list;
            });
        }
        fu.close();
        return drugToResults;
    }

    private void generateOutput(Map<String, List<DrugImpactResult>> drugToResults,
                                boolean needFilter,
                                String outFileName) throws IOException {
        if (needFilter)
            filterImpactResults(drugToResults);
        System.out.println("Total drugs having impact values: " + drugToResults.size());
        Set<String> pathways = drugToResults.values()
                                            .stream()
                                            .flatMap(results -> results.stream())
                                            .map(result -> result.pathway)
                                            .collect(Collectors.toSet());
        System.out.println("Total pathways: " + pathways.size());
        List<String> pathwayList = pathways.stream().sorted().collect(Collectors.toList());
        List<String> drugList = drugToResults.keySet().stream().sorted().collect(Collectors.toList());
        StringBuilder builder = new StringBuilder();
        builder.append("Drug");
        pathwayList.forEach(pathway -> builder.append("\t").append(pathway));
        fu.setOutput(outFileName);
        fu.printLine(builder.toString());
        builder.setLength(0);
        Map<String, Double> pathwayToValue = new HashMap<>();
        for (String drug : drugList) {
            builder.append(drug);
            pathwayToValue.clear();
            List<DrugImpactResult> results = drugToResults.get(drug);
            results.forEach(result -> pathwayToValue.put(result.pathway, result.value));
            pathwayList.forEach(pathway -> {
                Double value = pathwayToValue.get(pathway);
                builder.append("\t");
                if (value == null)
                    builder.append("0");
                else
                    builder.append(value);
            });
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private class DrugImpactResult {
        String drug;
        String pathway;
        Double value;
    }

}
