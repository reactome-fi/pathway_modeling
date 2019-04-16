/*
 * Created on Mar 28, 2007
 *
 */
package org.reactome.r3.util;


/**
 * This class holding constants used in analyze
 * @author guanming
 *
 */
public class R3Constants {
    // THe local database that is used as the data source
    public static final String REACTOME_SOURCE_DB_NAME = "reactome_40_plus_i";
    public static final String DB_USER = "root";
    public static final String DB_PWD = "macmysql01";
    //public static final String WORKING_DIR = "/Users/wgm/Documents/EclipseWorkspace/caBigR3/";
    public static final String WORKING_DIR = "/Users/wgm/Documents/EclipseWorkspace/caBigR3/";
    //public static final String RESULT_DIR = "results/v2/";
//    public static final String RESULT_DIR = "results/v3/";
//    public static final String RESULT_DIR = "results/v4/";
//    public static final String RESULT_DIR = "/Users/gwu/Documents/EclipseWorkspace/FINetworkBuild/results/2012/";
    public static final String RESULT_DIR = "/Users/wug/Documents/eclipse_workspace/FINetworkBuild/results/2016/";
//    public static final String RESULT_DIR = "results/FI_2012/";
    // These two files are for Gene expression data
//    public static final String LEE_GENE_EXP_FILE = "results/v3/PavlidisCoExp_Norm.txt";
    public static final String LEE_GENE_EXP_FILE = RESULT_DIR + "LeeGeneExp.txt";
//    public static final String PRIETO_GENE_EXP_FILE = "results/v3/CarlosCoExp_Norm.txt";
    public static final String PRIETO_GENE_EXP_FILE = RESULT_DIR + "PrietoGeneExp.txt";
    // FIs extracted from Reactome and converted pathway DBs
    public static final String REACTOME_FI_FILE = RESULT_DIR + "FIs_Reactome.txt";
    public static final String KEGG_FI_FILE = RESULT_DIR + "FIs_Kegg.txt";
    public static final String NCI_PID_FI_FILE = RESULT_DIR + "FIs_Pathway Interaction Database.txt";
    public static final String NCI_PID_BIOCARTA_FI_FILE = RESULT_DIR + "FIs_BioCarta - Imported by PID.txt";
    public static final String PANTHER_FI_FILE = RESULT_DIR + "FIs_pantherdb.txt";
    public static final String TRED_FI_FILE = RESULT_DIR + "FIs_TRED.txt";
    
    public static final String PREDICTED_FI_FILE = RESULT_DIR + "PredictedFIs_031512.txt";
    //public static final String PROTEIN_ID_TO_TOPIC = RESULT_DIR + "ProteinIdToTopic051109.txt";
//    public static final String PROTEIN_ID_TO_TOPIC = RESULT_DIR + "ProteinIdToTopic080410.txt";
    public static final String PROTEIN_ID_TO_TOPIC = RESULT_DIR + "ProteinIdToTopic101110.txt";
    public static final String GENE_TO_TOPIC = RESULT_DIR + "ProteinNameToTopics022717.txt";
    public static final String MCL_RESULT_DIR = RESULT_DIR + "mcl/";
    public static final String CLUSTER_RESULT_FILE = MCL_RESULT_DIR + "FI73_042108.class_I24_c26";
    //public static final String INTERACTION_FILE_NAME = MCL_RESULT_DIR  + "FIInteractions67.txt";
    //public static final String INTERACTION_FILE_NAME = RESULT_DIR  + "FIInteractions57.txt";
    public static final String INTERACTION_FILE_NAME = RESULT_DIR  + "FIs_043009.txt";
    public static final String INTERACTION_BIG_COMP_FILE_NAME = RESULT_DIR + "FIs_042109_BigComp.txt";
//    public static final String GENE_FI_BIG_COMP_FILE_NAME = RESULT_DIR + "FIsInGene_041709_BigComp.txt";
//    public static final String GENE_FI_BIG_COMP_FILE_NAME = RESULT_DIR + "FIsInGene_071012_BigComp.txt";
//    public static final String GENE_FI_BIG_COMP_FILE_NAME = RESULT_DIR + "FIsInGene_121013_BigComp.txt";
    public static final String GENE_FI_BIG_COMP_FILE_NAME = RESULT_DIR + "FIsInGene_022717_BigComp.txt";
//    public static final String GENE_FI_FILE_NAME = RESULT_DIR + "FIsInGene_041709.txt";
    public static final String GENE_FI_FILE_NAME = RESULT_DIR + "FIsInGene_071012.txt";
//    public static final String GENE_FI_FILE_NAME = RESULT_DIR + "FIsInGene_031612.txt";
    public static final String GENE_FI_PATHWAY_FILE_NAME = RESULT_DIR + "FIsInGene_Pathway_031612.txt";
    public static final String GENE_FI_PREDICTED_FILE_NAME = RESULT_DIR + "FIsInGene_Predicted_031612.txt";
//    public static final String GENE_FI_FILE_NAME = RESULT_DIR + "FIsInGene_No_ZNF_042810.txt";
//    public static final String GENE_FI_BIG_COMP_FILE_NAME = RESULT_DIR + "FIsInGene_No_ZNF_042810_BigComp.txt";
    
    //public static final String INTERACTION_FILE_NAME = RESULT_DIR  + "FI73_Pathway_042108.txt";
    //public static final String INTERACTION_FILE_NAME = RESULT_DIR  + "FI73_Predicated_042108.txt";
    // Negative dataset from random pairs from SwissProt
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "ReactomeInteractions020708_Random.model";
    // Negative dataset from different compartments
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "ReactomeInteractions020708_DiffCompart.model";
    // Positive and negative from Reactome only
    public static final String NBC_MODEL_NAME = RESULT_DIR + "ReactomeInteractions031908.model";
    // Learned model name
    // Six DBs but with random pairs
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions013108_SameIDsSizeRandom.model";
    // Four dbs (No INOH and Panther)
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "FourDBInteractions013108_SameIDs10Random_01Pos.model";
    // repeat using random pairs from all SwissProt identifiers
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions121307_2.model";
    // Repeat negative datasets is from pathways
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions121307_1.model";
    // Negative/positive ratio is based on pathways negative dataset. However, this negative
    // dataset was picked from nucleus and extracelluar region proteins
    // public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions121307.model";
    // Negative/positive ratio is around 1.0. Negative is picked from the same pathways
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions121207_1.model";
    // Negative dataset from distant cellular compartment (NoInteractions090506.txt)
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions121207.model";
    // Interactions remove HIVs
    // public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions061507.model";
    // Model from Xin using refined data set.
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "Xin092507.model";
//  Using random pairs as the negative data set
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions082807.model";
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions041807.model";
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractions020907.model";
    // A test using random pairs from all uniprot ids
    //public static final String NBC_MODEL_NAME = RESULT_DIR + "SixDBInteractionsWithRandomPairs.model";
    // Total line number of the above file (***57.txt)
    public static final int TOTAL_FI_INTERACTION = 194546;
    public static final int TOTAL_PROTEIN_IDS = 14346;
    public static final int TOTAL_HPRD_IDS = 25202;
    // The following number is from http://en.wikipedia.org/wiki/Human_genome
    public static final int TOTAL_HUMAN_GENES = 23000;
    // 25205 is the number of total proteins in HPRD. This number is
    // used as the total of human proteins
    public static final double COVERAGE = TOTAL_PROTEIN_IDS / (double) 25205;
    // Cutoff value used for NBC
    //public static final double CUT_OFF_VALUE = 0.73d;
    public static final double CUT_OFF_VALUE = 0.50d; // April, 2009.
    // The following are the list of interaction files 
    public static final String HUMAN_INTERACTION_FILE_NAME = RESULT_DIR + "HumanInteractions031808.txt";
    public static final String ORTHO_INTERAACTION_FILE_NAME = RESULT_DIR + "OrthoInteractions031908.txt";
    public static final String YEAST_INTERACTION_FILE_NAME = RESULT_DIR + "YeastInteractions031908.txt";
    public static final String GENE_EXP_FILE_NAME = RESULT_DIR + "GeneExpWith3FromPavlidis_041108.txt";
//    public static final String HUMAN_INTERACTION_FILE_NAME = RESULT_DIR + "HumanInteractions020507.txt";
//    public static final String ORTHO_INTERAACTION_FILE_NAME = "results/interaction/OrthoInteractions.txt";
//    public static final String YEAST_INTERACTION_FILE_NAME = "results/interaction/YeastInteractions.txt";
//    public static final String GENE_EXP_FILE_NAME = "results/microarray/GeneExpWith3FromPavlidis.txt";
    // the data sets directory
    public static final String DATA_SET_DIR = "datasets/";
    // Used to control the UniProt version
//    public static final String UNIPROT_DIR = DATA_SET_DIR + "UniProt" + File.separator + "release_14_9" + File.separator;
    public static final String UNIPROT_DIR = DATA_SET_DIR + "UniProt/release_2015_11/";
    
    // Directory for pFam
    public static final String PFAM_DIR_NAME = DATA_SET_DIR + "Pfam/26.0/";
    // Directory for RH file
    public static final String RH_DIR_NAME = DATA_SET_DIR + "GenomeResearch_RH_Map/fully_combined_RH_network/";
    // These parameters are related to KEGG pathways
    public static final String KEGG_DIR = DATA_SET_DIR + "KEGG/011112/";
    // Unzipped human KGML files should be in this directory
    public static final String KEGG_HSA_KGML_DIR = KEGG_DIR + "KGML/hsa/";
    // The converted KEGG pathways is saved in this project file
    public static final String KEGG_CONVERTED_FILE = KEGG_DIR + "030112.rtpj";
    // This file is used to map KEGG gene ids to UniProt ids
    public static final String KEGG_ID_TO_UNIPROT_MAP_FILE = KEGG_DIR + "hsa_genes_uniprot.list";
//    // KEGG parameters for 2009 version
//    public static final String KEGG_DIR = DATA_SET_DIR + "KEGG/031209/";
//    // Unzipped human KGML files should be in this directory
//    public static final String KEGG_HSA_KGML_DIR = KEGG_DIR + "hsa/";
//    // The converted KEGG pathways is saved in this project file
//    public static final String KEGG_CONVERTED_FILE = KEGG_DIR + "011612.rtpj";
//    // This file is used to map KEGG gene ids to UniProt ids
//    public static final String KEGG_ID_TO_UNIPROT_MAP_FILE = KEGG_DIR + "hsa_uniprot.list";
    
    // Used for the INOH directory
    public static final String INOH_DIR = DATA_SET_DIR + "INOH/AllDiagram_071128/";
    public static final String INOH_MOLECLE_ROLE_ONTOLOGY_FILE = DATA_SET_DIR + "INOH/MoleculeRoleOntology_222/MoleculeRoleOntology_222.obo";
    // Used for the panther database files
//    public static final String PANTHER_DIR = DATA_SET_DIR + "Panther/Version2.5/";
    public static final String PANTHER_DIR = DATA_SET_DIR + "Panther/Version3.0.1/"; // Download on Jan 18, 2011
    public static final String PANTHER_FILES_DIR = PANTHER_DIR + "SBML/";
    public static final String PANTHER_MAPPING_FILE = PANTHER_DIR + "SequenceAssociationPathway3.01.txt";
    public static final String PANTHER_CONVERTED_FILE = PANTHER_DIR + "Panther_3_0_1.rtpj";
    
    // Used for the Nature-PID database files
    public static final String NATURE_PID_DIR = DATA_SET_DIR + "NCI-Pathways/011612/";
    public static final String NATURE_PID_CURATED = NATURE_PID_DIR + "NCI-Nature_Curated.bp2.owl";
    public static final String NATURE_PID_CURATED_CONVERTED = NATURE_PID_DIR + "NCI-Nature_Curated.bp2.rtpj";
    public static final String NATURE_PID_BIOCARTA = NATURE_PID_DIR + "BioCarta.bp2.owl";
    public static final String NATURE_PID_BIOCARTA_CONVERTED = NATURE_PID_DIR + "BioCarta.bp2.rtpj";
//    public static final String NATURE_PID_DIR = DATA_SET_DIR + "NCI-Pathways/031709/";
    
    // This file is used to map Entrez id to UniProt accession number
    public static final String ENTREZ_TO_UNIPROT_MAP_FILE_NAME = DATA_SET_DIR + "iproclass/011612/EntrezToUniProt.txt";
    
    // For IntAct data set
//    public static final String INTACT_DIR = DATA_SET_DIR + "IntAct/032009/";
    public static final String INTACT_DIR = DATA_SET_DIR + "IntAct/022412/";
    public static final String INTACT_HUMAN_DIR = INTACT_DIR + "human/";
    // For BioGrid data set
    public static final String BIOGRID_DIR = DATA_SET_DIR + "BioGrid/2.0.50/";
    // For HPRD
    public static final String HPRD_DIR = DATA_SET_DIR + "HPRD/090107/";
    // For NCI_Nature
    public static final String NCI_NATURE_DIR = DATA_SET_DIR + "NCI-Pathways/031709/";
    // For OrthMCL
    public static final String ORTHO_MCL_DIR = DATA_SET_DIR + "OrthoMCL/v2/";
    // For geneways related data
    public static final String GENEWAYS_DIR = DATA_SET_DIR + "GeneWays/";
    // For ensembl related files
//    public static final String ENSEMBL_DIR = DATA_SET_DIR + "Ensembl/release_53/";
    public static final String ENSEMBL_DIR = DATA_SET_DIR + "Ensembl/release_62/";
    public static final String ENSEMBL_COMPARA_DATABASE = "ensembl_compara_62";
    public static final String ENSEMBL_PROTEIN_FAMILIES = ENSEMBL_DIR + "ProteinFamilies.txt";
    // For GO related files
//    public static final String GO_DIR = DATA_SET_DIR + "GO/032509/";
    public static final String GO_DIR = DATA_SET_DIR + "GO/030712/";
    // For TRED files
    public static final String TRED_DIR = DATA_SET_DIR + "TRED/";
    public static final String TRED_CONVERTED_FILE = TRED_DIR + "TRED_030112.rtpj";
    // For UCSC
//    public static final String UCSC_DIR = DATA_SET_DIR + "ucsc/072909/";
    public static final String UCSC_DIR = DATA_SET_DIR + "ucsc/080113/";
    
    public static final String NCBI_DIR = DATA_SET_DIR + "ncbi/080409/";
    
    public static final String BREAST_DIR = DATA_SET_DIR + "BreastCancer/";
    public static final String GBM_DIR = DATA_SET_DIR + "TCGA/GBM/";
    public static final String OVARIAN_DIR_NAME = DATA_SET_DIR + "TCGA/OvarianCancer/";
    
    public static final String mclScript = "/Users/gwu/Documents/EclipseWorkspace/caBigR3WebApp/WebContent/WEB-INF/mcl_script.sh";        
    public static final String survivalScript = "/Users/gwu/Documents/EclipseWorkspace/caBigR3WebApp/WebContent/WEB-INF/CGISurvivalAnalysis.R";
    public static final String R_SRC_DIR = "RSource/";
    public static final String TEMP_DIR = "tmp/";
    // Constants for iRefIndex data files
    public static final String IREFINDEX_DIR = DATA_SET_DIR + "iRefIndex/9.0/";
    public static final String IREFINDEX_HUMAN_FILE = IREFINDEX_DIR + "9606.mitab.10182011.txt";
    public static final String IREFINDEX_HUMAN_PPI_FILE = IREFINDEX_DIR + "HumanPPIsInUniProt022712.txt";
    public static final String IREFINDEX_YEAST_FILE = IREFINDEX_DIR  + "559292.mitab.10182011.txt";
    public static final String IREFINDEX_YEAST_PPI_FILE = IREFINDEX_DIR + "YeastPPIsInUniProt022812.txt";
    public static final String IREFINDEX_YEAST_TO_HUMAN_PPI_FILE = IREFINDEX_DIR + "HumanPPIsFromYeastInUniProt030112.txt";
    public static final String IREFINDEX_FLY_FILE = IREFINDEX_DIR + "7227.mitab.10182011.txt";
    public static final String IREFINDEX_FLY_PPI_FILE = IREFINDEX_DIR + "FLyPPIsInUniProt022812.txt";
    public static final String IREFINDEX_FLY_TO_HUMAN_PPI_FILE = IREFINDEX_DIR + "HumanPPIsFromFLyInUniProt030112.txt";
    public static final String IREFINDEX_WORM_FILE = IREFINDEX_DIR + "6239.mitab.10182011.txt";
    public static final String IREFINDEX_WORM_PPI_FILE = IREFINDEX_DIR + "WormPPIsInUniProt022812.txt";
    public static final String IREFINDEX_WORM_TO_HUMAN_PPI_FILE = IREFINDEX_DIR + "HumanPPIsFromWormInUniProt030112.txt";
    public static final String IREFINDEX_MOUSE_FILE = IREFINDEX_DIR + "10090.mitab.10182011.txt";
    public static final String IREFINDEX_MOUSE_PPI_FILE = IREFINDEX_DIR + "MousePPIsInUniProt031412.txt";
    public static final String IREFINDEX_MOUSE_TO_HUMAN_PPI_FILE = IREFINDEX_DIR + "HumanPPIsFromMouseInUniProt031412.txt";
    // For storing normalized PPIs mapped from non-human species
    public static final String HUMAN_PPIS_FROM_YEAST_FILE = RESULT_DIR + "HumanPPIsFromYeast030112.txt";
    public static final String HUMAN_PPIS_FROM_WORM_FILE = RESULT_DIR + "HumanPPIsFromWorm030112.txt";
    public static final String HUMAN_PPIS_FROM_FLY_FILE = RESULT_DIR + "HumanPPIsFromFly030112.txt";
    public static final String HUMAN_PPIS_FROM_MOUSE_FILE = RESULT_DIR + "HumanPPIsFromMouse031412.txt";
}
