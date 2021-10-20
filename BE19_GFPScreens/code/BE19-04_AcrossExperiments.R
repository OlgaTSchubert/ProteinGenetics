# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/"


# Processing directory
procdir <- "results/processing/"


# List of all results tables (BCstats.csv)
BCstats.files <- list.files(path = procdir, pattern = "BCstats.csv", 
                            recursive = T, full.names = T); BCstats.files


# List of experiment names (proteins)
experiments <- BCstats.files %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_processing/BCstats.csv)") %>%
        str_sort() %>% print()


# Import all results tables and store in a list
comblist <- map(BCstats.files, ~ read_csv(.))
names(comblist) <- experiments
saveRDS(comblist, paste0(procdir, "comblist.RDS"))


# Guide table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set %in% c("eProvs", "neProvs", "neStops")) %>%
        select(c("guide", "PAMstrand", "chr", "chr.start", "chr.end",
                 "mut1", "mut2", "mut3", "provean", "maxAbsProvean", 
                 "essential", "set", "geneSys", "gene")) %>%
        print()


# Gene annotations from SGD
sgd <- importExtData(dataset = "SGD_features", localfile = T) %>%
        select(Feature_name, Description) %>%
        rename(geneSys   = Feature_name, 
               geneDescr = Description) %>% print()


# Add gene descriptions to guides table
guides <- guides %>% left_join(sgd, by = "geneSys") %>% print()




# Combine GUIDE results across experiments -------------------------------------

# Combine results into a table
combdf.gd <- map(comblist, ~ select(., guide, log2fc, q)) %>%
        purrr::reduce(full_join, by = "guide") %>% print()


# Add protein names to column headers
nam1  <- rep(experiments, each = 2)
nam2  <- rep(c("log2fc", "q"), times = length(experiments))
names(combdf.gd)[2:(1+length(experiments)*2)] <- paste0(nam1, "_", nam2)


# Sort columns by type (log2fc, q) instead of experiment
combdf.gd <- relocate(combdf.gd, ends_with("_log2fc"), .after = "guide") %>% print()


# Add count of experiments for which the guide is a significant regulator
combdf.gd <- combdf.gd %>%
        mutate(FDR0.05_count = select(., ends_with("_q")) %>% 
                       apply(1, function(x) sum(x < 0.05, na.rm = T)))


# Add more annotations for each guide
combdf.gd <- right_join(guides, combdf.gd, by = "guide") %>% print()


# Save table
saveRDS(combdf.gd, paste0(procdir, "combdf_gd.RDS"))
write_excel_csv(combdf.gd, paste0(resdir, "combined_perGuide.csv"))




# Combining GENE results across experiments ------------------------------------

# Combine results into a table
combdf.gn <- comblist %>%
        map(~ .[!duplicated(.$geneSys), ]) %>%
        map(~ mutate(., essential = ifelse(set %in% c("eStops", "eProvs"), T, F))) %>%
        map(~ select(., geneSys, gene, essential, log2fc.of.q.min, q.comb)) %>%
        purrr::reduce(left_join, by = c("geneSys", "gene", "essential")) %>% 
        print()


# Add protein names to column headers
nam1  <- rep(experiments, each = 2)
nam2  <- rep(c("log2fc", "q"), times = length(experiments))
names(combdf.gn)[4:(3+length(experiments)*2)] <- paste0(nam1, "_", nam2)


# Sort columns by type (instead of experiment)
combdf.gn <- relocate(combdf.gn, ends_with("_log2fc"), .after = "essential") %>% print()


# Add count of experiments for which the gene is a significant regulator
combdf.gn <- combdf.gn %>%
        mutate(FDR0.05_count = select(., ends_with("_q")) %>% 
                       apply(1, function(x) sum(x < 0.05, na.rm = T)))


# Add more annotations for each gene
combdf.gn <- right_join(sgd, combdf.gn, by = "geneSys") %>% 
        relocate(gene, .after = geneSys) %>%print()


# Save table
saveRDS(combdf.gn, paste0(procdir, "combdf_gn.RDS"))
write_excel_csv(combdf.gn, paste0(resdir, "combined_perGene.csv"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-04_SessionInfo.txt")

