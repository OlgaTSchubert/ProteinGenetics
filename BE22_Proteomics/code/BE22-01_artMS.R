# Libraries --------------------------------------------------------------------

library(tidyverse)
library(artMS)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Results directory
resdir <- "results/artMS/"
dir.create(resdir)


# QC results directory
qcdir <- "results/02_artMS_v4_QC/"
dir.create(qcdir)


# artMS input files
read.delim2("results/01_MaxQuant_v4/evidence.txt") %>% 
        filter(!str_detect(Proteins, ";")) %>%
        write_delim("results/01_MaxQuant_v4/evidence_uniqueOnly.txt", delim = "\t")
evidencefile  <- "results/01_MaxQuant_v4/evidence_uniqueOnly.txt"
summaryfile   <- "results/01_MaxQuant_v4/summary.txt"
configfile    <- "results/02_artMS_shared/config_v4_all_quantile.yaml"  # To get config template: artmsWriteConfigYamlFile()
keysfile      <- "results/02_artMS_shared/keys.txt"
contrastsfile <- "results/02_artMS_shared/contrasts_all.txt"




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
        verbose = TRUE
)




# Analysis ---------------------------------------------------------------------

# Relative quantification
#?artmsQuantification
artmsQuantification(yaml_config_file = configfile,
                    display_msstats = TRUE)


# Further analyses
# artMS asks to install org.Sc.sgd.db: install.packages('org.Sc.sgd.db')
# This gave an error (related to R version), but installation from Bioconductor worked.
# However, when loading the library (library(org.Sc.sgd.db)) and running the command,
# there is an error regarding a missing column "SYMBOL". Indeed there is no such
# column in the yeast db file: # columns(org.Sc.sgd.db)
# Had to remove.packages("org.Sc.sgd.db") and restart R before I could run it again.

?artmsAnalysisQuantifications
artmsAnalysisQuantifications(
        log2fc_file  = paste0(resdir, "results.txt"),
        modelqc_file = paste0(resdir, "results_ModelQC.txt"),
        output_dir   = "results",
        
        # Enrichment - only for HUMAN or MOUSE
        species = "YEAST",  # This is for enrichment, which currently only supported for HUMAN and MOUSE
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




# Reformatting artMS output to wide table --------------------------------------
# The pre-printed results-summary.xlsx contains 0 instead of NA for 
# log2FC and p-values of missing comparisons - misleading!)

res <- read.delim2(paste0(resdir, "results_adjpvalue/results-log2fc-long.txt")) %>% 
        as_tibble() %>% 
        mutate(mutant = str_sub(Comparison, 1, -4)) %>%
        left_join(sgd, by = c("Protein" = "geneSys")) %>%
        mutate(gfp = ifelse(gene %in% prots, T, F),
               adj.pvalue = ifelse(is.na(pvalue), NA, adj.pvalue)) %>%
        select(Protein, gene, gfp, Description, mutant, log2FC, adj.pvalue, 
               iLog2FC, iPvalue) %>%
        rename(geneSys  = "Protein", 
               descr    = "Description",
               log2FC   = "log2FC",
               log2FC.i = "iLog2FC",
               adj.p    = "adj.pvalue",
               adj.p.i  = "iPvalue") %>%
        mutate_at(vars(log2FC, log2FC.i, adj.p, adj.p.i), list(as.numeric)) %>%
        print()

res.wide <- res %>%
        mutate(mutant = factor(mutant, levels = muts)) %>%
        arrange(mutant) %>%
        pivot_wider(names_from = mutant,
                    names_sep = "_",
                    values_from = c(log2FC, adj.p, log2FC.i, adj.p.i),
                    values_fill = NA) %>% print()
write_delim(res.wide, paste0(resdir, "artMS_results_OS.tsv"), delim = "\t")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-01_SessionInfo.txt")

