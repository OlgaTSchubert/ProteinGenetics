# Libraries --------------------------------------------------------------------

library(tidyverse)
library(artMS)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Functions
source("../BE00_gRNALibraryDesign/code/BE00-00_ImportExtData.R")


# artMS results directory
resdir <- "results/artMS/"
dir.create(resdir)


# artMS QC results directory
qcdir <- "results/artMS/QC/"
dir.create(qcdir)


# MaxQuant files and results directory
mqdir <- "results/MaxQuant/"
dir.create(mqdir)
# Now manually copy evidence.txt and summary.txt from MQ into this folder

evidencefile_all <- "results/MaxQuant/evidence.txt"
evidencefile     <- "results/MaxQuant/evidence_uniqueOnly.txt"
summaryfile      <- "results/MaxQuant/summary.txt"


# Import MaxQuant output and filter for unique peptides only
evid_all <- read.delim2(evidencefile_all)
evid     <- evid_all %>% filter(!str_detect(Proteins, ";"))
write_delim(evid, evidencefile, delim = "\t")


# artMS files
configfile    <- "code/BE22-01_artMS_config.yaml"  # To get config template: artmsWriteConfigYamlFile()
keysfile      <- "code/BE22-01_artMS_keys.txt"
contrastsfile <- "code/BE22-01_artMS_contrasts.txt"


# Gene annotations from Uniprot and SGD
sgd <- importExtData(dataset = "SGD_features", localfile = T) %>%
        select(Feature_name, Standard_gene_name) %>%
        dplyr::rename(geneSys = Feature_name, 
                      gene    = Standard_gene_name) %>%
        mutate(gene = ifelse(is.na(gene), geneSys, gene)) %>% print()


# Proteins and mutations lists
muts  <- read_delim("annotations/MutationLookup.tsv", delim = "\t") %>% print()




# QC ---------------------------------------------------------------------------

#?artmsQualityControlEvidenceBasic
artmsQualityControlEvidenceBasic(
        evidence_file = evidencefile, 
        keys_file = keysfile, 
        prot_exp = "AB",
        fractions = 0,
        output_dir = qcdir,
        output_name = "qcBasic",
        isSILAC = FALSE,
        plotINTDIST = TRUE,
        plotREPRO = TRUE,
        plotCORMAT = TRUE,
        plotINTMISC = TRUE,
        plotPTMSTATS = TRUE,
        printPDF = TRUE,
        verbose = TRUE)


#?artmsQualityControlEvidenceExtended
artmsQualityControlEvidenceExtended(
        evidence_file = evidencefile, 
        keys_file = keysfile, 
        output_dir = qcdir,
        output_name = "qcExtended",
        isSILAC = FALSE,
        plotPSM = TRUE,
        plotIONS = TRUE,
        plotTYPE = TRUE,
        plotPEPTIDES = TRUE,
        plotPEPTOVERLAP = TRUE,
        plotPROTEINS = TRUE,
        plotPROTOVERLAP = TRUE,
        plotPIO = TRUE,
        plotCS = TRUE,
        plotME = TRUE,
        plotMOCD = TRUE,
        plotPEPICV = TRUE,
        plotPEPDETECT = TRUE,
        plotPROTICV = TRUE,
        plotPROTDETECT = TRUE,
        plotIDoverlap = TRUE,
        plotPCA = FALSE,       # Takes very long
        plotSP = TRUE,
        printPDF = TRUE,
        verbose = TRUE)


#?artmsQualityControlSummaryExtended
artmsQualityControlSummaryExtended(
        summary_file = summaryfile, 
        keys_file = keysfile, 
        output_dir = qcdir,
        output_name = "qcExtended_summary",
        isFractions = FALSE,
        plotMS1SCANS = TRUE,
        plotMS2 = TRUE,
        plotMSMS = TRUE,
        plotISOTOPE = TRUE,
        printPDF = TRUE,
        verbose = TRUE)




# Analysis ---------------------------------------------------------------------

# Relative quantification
#?artmsQuantification
artmsQuantification(yaml_config_file = configfile,
                    display_msstats = TRUE)


# Further analyses
#?artmsAnalysisQuantifications
artmsAnalysisQuantifications(
        log2fc_file  = paste0(resdir, "artMS_results.txt"),
        modelqc_file = paste0(resdir, "artMS_results_ModelQC.txt"),
        output_dir   = "results",
        species = "YEAST",  # This is for enrichment, which is currently only supported for HUMAN and MOUSE
        enrich = FALSE,     # This is automatically disabled for YEAST
        l2fc_thres = 1,
        isBackground = "nobackground",
        outliers = "keep",           # c("keep", "iqr", "std")
        choosePvalue = "adjpvalue",  # c("adjpvalue", "pvalue")
        isPtm = "global",            # c("global", "ptmsites")
        mnbr = 2,
        pathogen = "nopathogen",
        plotPvaluesLog2fcDist = TRUE,
        plotAbundanceStats = TRUE,
        plotReproAbundance = TRUE,
        plotCorrConditions = TRUE,
        plotCorrQuant = TRUE,
        plotPCAabundance = TRUE,
        plotFinalDistributions = TRUE,
        plotPropImputation = TRUE,
        plotHeatmapsChanges = TRUE,
        plotTotalQuant = TRUE,
        plotClusteringAnalysis = TRUE,
        data_object = FALSE,
        verbose = TRUE)




# Reformatting artMS output ----------------------------------------------------

res <- read_delim(paste0(resdir, "results_adjpvalue/artMS_results-log2fc-long.txt"), delim = "\t") %>% 
        as_tibble() %>%
        rename(geneSys  = "Protein") %>%
        left_join(sgd, by = "geneSys") %>% 
        mutate(mutant = str_sub(Comparison, 1, -4)) %>% 
        mutate(adj.pvalue = ifelse(is.na(pvalue), NA, adj.pvalue)) %>%
        rename(log2FC   = "log2FC",
               log2FC.i = "iLog2FC",
               adj.p    = "adj.pvalue",
               adj.p.i  = "iPvalue") %>%
        select(geneSys, gene, mutant, log2FC, adj.p, log2FC.i, adj.p.i) %>%
        mutate_at(vars(log2FC, log2FC.i, adj.p, adj.p.i), list(as.numeric)) %>%
        filter(mutant != "SAP155") %>%
        print()
write_delim(res, paste0(resdir, "artMS_results_OS.tsv"), delim = "\t")

res.wide <- res %>%
        mutate(mutant = factor(mutant, levels = muts$alt)) %>%
        arrange(mutant) %>%
        pivot_wider(names_from = mutant,
                    names_sep = "_",
                    values_from = c(log2FC, adj.p, log2FC.i, adj.p.i),
                    values_fill = NA) %>% print()
write_delim(res.wide, paste0(resdir, "artMS_results_OS_wide.tsv"), delim = "\t")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-01_SessionInfo_NEW.txt")

