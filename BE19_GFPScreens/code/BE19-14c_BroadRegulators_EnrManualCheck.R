# Libraries --------------------------------------------------------------------

library(tidyverse)
library(org.Sc.sgd.db)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()




# Get genes for each GO IDs of interest ----------------------------------------

# GO:0042254 (Ribosome biogenesis)
ribo <- AnnotationDbi::select(org.Sc.sgd.db, columns = "ENSEMBL",
                              keytype = "GOALL", keys = "GO:0042254") %>%
        filter(!is.na(ENSEMBL)) %>%
        pull(ENSEMBL) %>% unique() %>% print()

# GO:0006399 (tRNA metabolic process)
tRNA <- AnnotationDbi::select(org.Sc.sgd.db, columns = "ENSEMBL",
                              keytype = "GOALL", keys = "GO:0006399") %>%
        filter(!is.na(ENSEMBL)) %>%
        pull(ENSEMBL) %>% unique() %>% print()

# GO:0006351 (DNA-templated transcription)
trxn <- AnnotationDbi::select(org.Sc.sgd.db, columns = "ENSEMBL",
                              keytype = "GOALL", keys = "GO:0006351") %>%
        filter(!is.na(ENSEMBL)) %>%
        pull(ENSEMBL) %>% unique() %>% print()




# Generate input lists ---------------------------------------------------------

backg <- combdf.gn %>%
        pull(geneSys) %>%
        print()

broad <- combdf.gn %>%
        filter(FDR0.05_count > 7) %>%
        pull(geneSys) %>%
        print()

speci <- combdf.gn %>%
        filter(FDR0.05_count %in% c(1,2)) %>%
        pull(geneSys) %>%
        print()

backg_ribo <- intersect(backg, ribo) %>% print()
backg_tRNA <- intersect(backg, tRNA) %>% print()
backg_trxn <- intersect(backg, trxn) %>% print()

broad_ribo <- intersect(broad, ribo) %>% print()
broad_tRNA <- intersect(broad, tRNA) %>% print()
broad_trxn <- intersect(broad, trxn) %>% print()

speci_ribo <- intersect(speci, ribo) %>% print()
speci_tRNA <- intersect(speci, tRNA) %>% print()
speci_trxn <- intersect(speci, trxn) %>% print()




# Test for over- and underrepresentation ---------------------------------------

# Note that Fisher's exact test and hypergeometric test are identical

# Hypergeometric test for over/underrepresentation
# phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE)
# phyper(x, m, n, k, lower.tail = T, log.p = FALSE)

# Fisher's exact test for over/underrepresentation
# fisher.test(matrix(c(x, k-x, m-x, n-(k-x)), nrow = 2), alternative = "greater")$p.value
# fisher.test(matrix(c(x, k-x, m-x, n-(k-x)), nrow = 2), alternative = "less")$p.value



# Broad regulators and ribosome biogenesis
m <- length(backg_ribo)                    # Genes in GO term
n <- length(backg) - length(backg_ribo)    # Genes not in GO term
k <- length(broad)                         # Differentially expressed genes
x <- length(broad_ribo)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  7.098989e-08
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 1


# Broad regulators and tRNA metabolic process
m <- length(backg_tRNA)                    # Genes in GO term
n <- length(backg) - length(backg_tRNA)    # Genes not in GO term
k <- length(broad)                         # Differentially expressed genes
x <- length(broad_tRNA)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  4.606359e-05
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 0.9999957


# Broad regulators and DNA-templated transcription
m <- length(backg_trxn)                    # Genes in GO term
n <- length(backg) - length(backg_trxn)    # Genes not in GO term
k <- length(broad)                         # Differentially expressed genes
x <- length(broad_trxn)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  0.9182247
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 0.2290474


# Specific regulators and ribosome biogenesis
m <- length(backg_ribo)                    # Genes in GO term
n <- length(backg) - length(backg_ribo)    # Genes not in GO term
k <- length(speci)                         # Differentially expressed genes
x <- length(speci_ribo)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  4.333058e-05
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 0.9999791


# Specific regulators and tRNA metabolic process
m <- length(backg_tRNA)                    # Genes in GO term
n <- length(backg) - length(backg_tRNA)    # Genes not in GO term
k <- length(speci)                         # Differentially expressed genes
x <- length(speci_tRNA)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  0.01273385
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 0.9931845


# Specific regulators and DNA-templated transcription
m <- length(backg_trxn)                    # Genes in GO term
n <- length(backg) - length(backg_trxn)    # Genes not in GO term
k <- length(speci)                         # Differentially expressed genes
x <- length(speci_trxn)                    # Differentially expressed genes in GO term
phyper(x-1, m, n, k, lower.tail = F, log.p = FALSE) # Overrepresentation:  5.556805e-12
phyper(x, m, n, k, lower.tail = T, log.p = FALSE)   # Underrepresentation: 1




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-14c_SessionInfo.txt")

