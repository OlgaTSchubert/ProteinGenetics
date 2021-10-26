# Libraries --------------------------------------------------------------------

library(tidyverse)
library(STRINGdb)   # To see all available functions: STRINGdb$methods()




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Results directory
resdir <- "results/STRING-enrichments/"
dir.create(resdir)


# Import annotations
muts  <- read_delim("annotations/MutationLookup.tsv", delim = "\t") %>% 
        pull(alt) %>% print()


# Import artMS results
res <- read_delim("results/artMS/artMS_results_OS_wide.tsv", delim = "\t") %>% print()




# Prepare to run STRING --------------------------------------------------------

# Initialize STRING DB
string_db <- STRINGdb$new(version = "11", 
                          species = 4932,
                          score_threshold = 400, 
                          input_directory = "../ExternalData/STRING/")


# Add STRING identifiers to results (doesn't work on tibble)
#STRINGdb$help("map")
res.m <- string_db$map(my_data_frame = as.data.frame(res), 
                       my_data_frame_id_col_names = "geneSys", 
                       takeFirst = TRUE, 
                       removeUnmappedRows = TRUE) %>%
        as_tibble() %>% 
        print()


# Define background (just take all proteins quantified across all samples)
bg <- res.m$STRING_id
string_db$set_background(background_vector = bg)




# Generate STRING input lists --------------------------------------------------

# List of vectors of sig. up- or down-regulated proteins in each sample
hits        <- vector(mode = "list", length = length(muts)*2)
names(hits) <- paste0(rep(muts, each = 2), "_", rep(c("up", "dn"), length(muts)))
str(hits)

for (i in seq_along(names(hits))) {
        
        if (str_detect(names(hits)[i], "up")) {
                hits[[i]] <- res %>% filter(get(str_c("log2FC.i_", str_sub(names(hits)[i], 1, -4))) > 1,
                                            get(str_c("adj.p.i_", str_sub(names(hits)[i], 1, -4))) < 0.05) %>%
                        pull(geneSys)
        } else {
                hits[[i]] <- res %>% filter(get(str_c("log2FC.i_", str_sub(names(hits)[i], 1, -4))) < -1,
                                            get(str_c("adj.p.i_", str_sub(names(hits)[i], 1, -4))) < 0.05) %>%
                        pull(geneSys)
        }
}
str(hits)




# Enrichment analysis ----------------------------------------------------------

enr        <- vector(mode = "list", length = length(muts)*2)
names(enr) <- paste0(rep(muts, each = 2), "_", rep(c("up", "dn"), length(muts)))
str(enr)

for (i in seq_along(names(enr))) {

        print(i)
        
        # Enrichment analysis
        #STRINGdb$help("get_enrichment")
        enr[[i]] <- string_db$get_enrichment(hits[[i]]) %>%
                as_tibble() %>%
                mutate_at(vars(category, term, inputGenes, preferredNames, description), list(as.character)) %>%
                mutate_at(vars(number_of_genes, number_of_genes_in_background, ncbiTaxonId), list(as.integer)) %>%
                mutate_at(vars(p_value, fdr), list(as.double))
}


# Write individual files with enrichments
names(enr) %>% map(~ write_tsv(enr[[.]], file = paste0(resdir, "enriched_", ., ".tsv")))


# Combine all results into a big table
enr.tb <- map(enr, function(x) dplyr::select(x, "category", "term", "description",
                                             "number_of_genes_in_background",
                                             "number_of_genes", "fdr", 
                                             "inputGenes", "preferredNames")) %>%
        purrr::reduce(full_join, by = c("category", "term", "description", 
                                        "number_of_genes_in_background")) %>%
        as_tibble()

names(enr.tb)[5:ncol(enr.tb)] <- paste0(rep(names(enr), each = 4), "_", 
                                        names(enr.tb)[5:ncol(enr.tb)]) %>%
        str_replace_all("\\..*", "")


# Add counts for significant enrichments across samples
enr.tb <- enr.tb %>%
        mutate(up_sig_005 = rowSums(dplyr::select(., ends_with("_up_fdr")) < 0.05, na.rm = T),
               dn_sig_005 = rowSums(dplyr::select(., ends_with("_dn_fdr")) < 0.05, na.rm = T),
               up_sig_001 = rowSums(dplyr::select(., ends_with("_up_fdr")) < 0.01, na.rm = T),
               dn_sig_001 = rowSums(dplyr::select(., ends_with("_dn_fdr")) < 0.01, na.rm = T)) %>% 
        print()


# Keep only useful categories
enr.tb <- enr.tb %>%
        filter(category %in% c("NetworkNeighborAL", "KEGG", 
                               "Process", "Component"))


# Write file
write_tsv(enr.tb, file = paste0(resdir, "enriched_all.tsv"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-02_SessionInfo.txt")

