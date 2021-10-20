# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE21_FitnessScreen")
getwd()


# Data directory
datadir <- "data/"
dir.create(datadir)


# Reads directory
readsdir <- "~/BigDataFiles/2021_ProteinGenetics/BE21/ReadsRDS/"


# Processed reads files
reads.files <- list.files(readsdir, pattern = ".*_Reads.RDS", full.names = T); reads.files


# List of experiment names
experiments <- reads.files %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_Reads.RDS)") %>%
        str_sort() %>%
        print()


# Guides table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set != "Can1Ade2") %>%
        select(guide, set, geneSys, gene) %>%
        print()




# Read counts per guide --------------------------------------------------------

# Initiate results table (based on guides table)
res <- guides


# Add read counts to results table
for(exp in seq_along(experiments)) {
        
        print(experiments[exp])
        
        samplesDF <- readRDS(paste0(readsdir, experiments[exp], "_Reads.RDS"))
        reads     <- samplesDF$results
        
        for (i in seq_along(reads)) {
                
                print(names(reads)[i])
                
                readsX <- reads[[i]] %>%
                        as_tibble() %>%
                        group_by(guide) %>%
                        summarize(n_reads = n()) 
                
                res <- left_join(res, readsX, by = "guide") %>% 
                        rename_at(vars(starts_with("n_reads")), ~names(reads)[i])
        }
}; res


# Save results tables
saveRDS(res, paste0(datadir, "readCounts.RDS"))

res_wES  <- res %>% select(-c(T00h_woES_A:Tpre_woES_B))
saveRDS(res_wES, paste0(datadir, "readCounts_wES.RDS"))

res_woES <- res %>% select(-c(T00h_wES_A:Tpre_wES_B))
saveRDS(res_woES, paste0(datadir, "readCounts_woES.RDS"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-02_SessionInfo.txt")

