# Libraries --------------------------------------------------------------------

library(tidyverse)
library(STRINGdb)   # To see all available functions: STRINGdb$methods()




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Results directory
resdir <- "results/enrichments/"
dir.create(resdir)


# Import annotations
muts <- c("SSY5", "POP1", "SIT4") %>% print()


# Import artMS results
res <- read_delim("results/artMS/artMS_results_OS_wide.tsv", delim = "\t") %>% 
        select(geneSys, gene, 
               log2FC.i_SSY5, adj.p.i_SSY5,
               log2FC.i_POP1, adj.p.i_POP1,
               log2FC.i_SIT4, adj.p.i_SIT4) %>%
        filter_at(vars(starts_with("adj.p.i_")), any_vars(!is.na(.))) %>%
        print()




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


# Define background (all proteins quantified across all samples)
bg <- res.m$STRING_id
string_db$set_background(background_vector = bg)
length(bg) # 2819




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
                mutate_at(vars(p_value, fdr), list(as.double)) %>%
                mutate(term = str_replace(term, "\\.", ":")) %>%
                mutate(nall = length(bg)) %>%
                mutate(nset = length(hits[[i]])) %>%
                rename(GOID = term,
                       GOterm = description,
                       nhits  = number_of_genes,
                       nrefs  = number_of_genes_in_background,
                       hits   = inputGenes,
                       hits2  = preferredNames,
                       p      = p_value,
                       fdr    = fdr) %>%
                select(category, GOID, GOterm, nall, nset, nhits, nrefs, p, fdr, hits, hits2)
}
str(enr)


# Save results
saveRDS(enr, paste0(resdir, "enrichments.RDS"))
#enr <- readRDS(paste0(resdir, "enrichments.RDS")); str(enr)


# Write individual files with enrichments
names(enr) %>% map(~ write_tsv(enr[[.]], file = paste0(resdir, ., "_enr.tsv")))


# Combine all results into a big table
enr.tb <- enr %>% 
        map(~ select(., category, GOID, GOterm, nall, nrefs, nset, nhits, fdr, hits, hits2)) %>%
        purrr::reduce(full_join, by = c("category", "GOID", "GOterm", "nall", "nrefs")) %>%
        print()

names(enr.tb)[6:ncol(enr.tb)] <- paste0(rep(names(enr), each = 5), "_", 
                                        names(enr.tb)[6:ncol(enr.tb)]) %>%
        str_replace_all("\\..*", "")
enr.tb


# Filter categories and write to file
enr.tb.GO <- enr.tb %>% filter(category %in% c("Process", "Component"))
write_tsv(enr.tb.GO, file = paste0(resdir, "combined_enr.tsv"))




# Manual selection of GO terms for plotting ------------------------------------

selectedGOs <- rbind(enr[["SSY5_up"]] %>% filter(GOID == "GO:0006520" & GOterm == "cellular amino acid metabolic process"),
                     enr[["SSY5_dn"]] %>% filter(GOID == "GO:0098798" & GOterm == "mitochondrial protein complex"),
                     enr[["POP1_up"]] %>% filter(GOID == "GO:0008652" & GOterm == "cellular amino acid biosynthetic process"),
                     enr[["POP1_dn"]] %>% filter(GOID == "GO:0005730" & GOterm == "nucleolus"),
                     enr[["SIT4_up"]] %>% filter(GOID == "GO:0045333" & GOterm == "cellular respiration"),
                     enr[["SIT4_up"]] %>% filter(GOID == "GO:0005975" & GOterm == "carbohydrate metabolic process"),
                     enr[["SIT4_up"]] %>% filter(GOID == "GO:0006979" & GOterm == "response to oxidative stress")) %>%
        mutate(group = c("SSY5_up", "SSY5_dn", "POP1_up", "POP1_dn", "SIT4_up", "SIT4_up", "SIT4_up")) %>%
        mutate(dir = str_remove(group, ".*_"),
               mutant = str_remove(group, "_..")) %>%
        select(group, mutant, dir, category, GOID, GOterm, nall, nrefs, nset, nhits, fdr, hits) %>%
        print()
write_delim(selectedGOs, paste0(resdir, "selectedGOs.tsv"), delim = "\t")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-02_SessionInfo.txt")

