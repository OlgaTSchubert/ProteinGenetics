# Libraries --------------------------------------------------------------------

library(tidyverse)
library(STRINGdb)   # To see all available functions: STRINGdb$methods()




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/broadRegulators_Enrichments/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()




# Prepare to run STRING --------------------------------------------------------

# Initialize STRING DB
string_db <- STRINGdb$new(version = "11", 
                          species = 4932,
                          score_threshold = 400, 
                          input_directory = "../ExternalData/STRING/")


# Add STRING identifiers to results (doesn't work on tibble)
#STRINGdb$help("map")
res.m <- string_db$map(my_data_frame = as.data.frame(combdf.gn), 
                       my_data_frame_id_col_names = "geneSys", 
                       takeFirst = TRUE, 
                       removeUnmappedRows = TRUE) %>%
    as_tibble() %>% 
    print()


# Define background (all proteins quantified across all samples)
bg <- res.m$STRING_id
string_db$set_background(background_vector = bg)
length(bg) # 4511




# Generate STRING input lists --------------------------------------------------

broadregs <- combdf.gn %>%
    filter(FDR0.05_count > 7) %>%
    pull(geneSys) %>%
    print()

specregs <- combdf.gn %>%
    filter(FDR0.05_count %in% c(1,2)) %>%
    pull(geneSys) %>%
    print()




# Enrichment analysis ----------------------------------------------------------

broadregs.enr <- string_db$get_enrichment(broadregs) %>%
    as_tibble() %>% 
    mutate(term = str_replace(term, "\\.", ":")) %>%
    mutate(nall = length(bg)) %>%
    mutate(nset = length(broadregs)) %>%
    rename(GOID = term,
           GOterm = description,
           nhits  = number_of_genes,
           nrefs  = number_of_genes_in_background,
           hits   = inputGenes,
           hits2  = preferredNames,
           p      = p_value,
           fdr    = fdr) %>%
    select(category, GOID, GOterm, nall, nset, nhits, nrefs, p, fdr, hits, hits2) %>%
    print()
write_tsv(broadregs.enr, file = paste0(resdir, "broadregs_enr.tsv"))

specregs.enr <- string_db$get_enrichment(specregs) %>%
    as_tibble() %>% 
    mutate(term = str_replace(term, "\\.", ":")) %>%
    mutate(nall = length(bg)) %>%
    mutate(nset = length(specregs)) %>%
    rename(GOID = term,
           GOterm = description,
           nhits  = number_of_genes,
           nrefs  = number_of_genes_in_background,
           hits   = inputGenes,
           hits2  = preferredNames,
           p      = p_value,
           fdr    = fdr) %>%
    select(category, GOID, GOterm, nall, nset, nhits, nrefs, p, fdr, hits, hits2) %>%
    print()
write_tsv(specregs.enr, file = paste0(resdir, "specregs_enr.tsv"))


# Combine results into a big table
combined <- full_join(broadregs.enr, specregs.enr, 
                      by = c("category", "GOID", "GOterm", "nall", "nrefs")) %>%
    relocate(nrefs, .after = nall) %>%
    select(-starts_with("p.")) %>%
    print()

names(combined)[6:ncol(combined)] <- paste0(rep(c("broadregs", "specregs"), each = 5), "_", 
                                            names(combined)[6:ncol(combined)]) %>%
    str_replace_all("\\..*", "")
combined


# Filter categories and write to file
combined.GO <- combined %>% filter(category == "Process") %>% print()
write_tsv(combined.GO, file = paste0(resdir, "combined_enr.tsv"))




# Manual selection of GO terms for plotting ------------------------------------

selectedGOs <- rbind(broadregs.enr %>% filter(GOID == "GO:0042254" & GOterm == "ribosome biogenesis"),
                     broadregs.enr %>% filter(GOID == "GO:0006399" & GOterm == "tRNA metabolic process"),
                     specregs.enr  %>% filter(GOID == "GO:0006351" & GOterm == "transcription, DNA-templated"),
                     specregs.enr  %>% filter(GOID == "GO:0006325" & GOterm == "chromatin organization")) %>%
    mutate(group = c("broad", "broad", "spec", "spec")) %>%
    relocate(group) %>%
    print()
write_delim(selectedGOs, paste0(resdir, "selectedGOs.tsv"), delim = "\t")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-14b_SessionInfo.txt")

