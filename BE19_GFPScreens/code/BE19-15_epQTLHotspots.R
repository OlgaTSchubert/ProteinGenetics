# Libraries --------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)  # Beware of function names overlapping with tidyverse!




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/epQTLHotspots/"
dir.create(resdir)


# Import data
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()

cols_log2fc <- str_subset(names(combdf.gn), "_log2fc")
cols_qval   <- str_subset(names(combdf.gn), "_q")
experiments <- str_sub(cols_qval, start = 1, end = -3) %>% print()


# Get for each target gene the list of proteins it affects significantly
prots <- combdf.gn %>%
        select(-geneDescr, -essential) %>%
        pivot_longer(-c(geneSys, gene, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        mutate(sign = ifelse(log2fc > 0, "pos", "neg")) %>%
        filter(q < 0.05) %>%
        group_by(geneSys, gene) %>%
        summarize(FDR0.05_prots = paste(toupper(protein), collapse = ";")) %>% print()


# Add this info to the big results table
combdf.gn <- combdf.gn %>% left_join(prots, by = c("geneSys", "gene")) %>% print()


# Chromosome lengths
sgdchr <- importExtData(dataset = "SGD_chromlengths", localfile = T) %>% 
        select(-NCBI_RefSeq_accession_number) %>%
        mutate(chrlencum = lag(cumsum(length))) %>% 
        mutate(chrlencum = if_else(is.na(chrlencum), 0, chrlencum)) %>% print()


# Gene annotations
sgdpos <- importExtData(dataset = "SGD_features", localfile = T) %>%
        dplyr::rename(geneSys = Feature_name,
                      gene    = Standard_gene_name,
                      chr     = Chromosome,
                      start   = Start_coordinate, 
                      end     = Stop_coordinate) %>%
        select(geneSys, gene, chr, start, end) %>%
        left_join(sgdchr, by = c("chr" = "chromosome")) %>%
        mutate(gstart = ifelse(start < end, start + chrlencum, end + chrlencum),
               gene = ifelse(is.na(gene), geneSys, gene)) %>% 
        print()

sgd <- sgdpos %>% select(geneSys, gene) %>% print()


# Gene annotations as txdb GRanges object
txdb <- importExtData(dataset = "txdb", localfile = T)
txdb <- renameSeqlevels(txdb, c(as.character(1:16), "M"))
txdb <- txdb[names(txdb) %in% sgd$geneSys] 
txdb


# BYxRM variants
RMvar <- importExtData(dataset = "BYxRM_variants", localfile = T) %>% 
        filter(IMPACT %in% c("MODERATE", "HIGH")) %>% 
        mutate(Name = ifelse(SYMBOL == "-", Gene, SYMBOL)) %>%
        unique() %>% arrange(Name) %>% print()




# Prepare BE hotspot table for comparisons -------------------------------------

combdf.gnx <- combdf.gn %>%
        left_join(sgdpos, by = c("gene", "geneSys")) %>%
        rename(gene = "Name") %>%
        select(Name, geneSys, chr, gstart, FDR0.05_count, FDR0.05_prots) %>% 
        filter(FDR0.05_count > 0) %>%
        arrange(gstart) %>% print()




# Find overlap with eLife 2018 eQTL hotspots -----------------------------------

# Import eQTL hotspots: Albert & Bloom et al., eLife 2018
elife2018        <- importExtData(dataset = "Albert2018_hotspots", localfile = T) %>% print()
elife2018$HSID   <- seq_along(elife2018$hotspotMarker)
elife2018$HSName <- NA


# Known and newly fine-mapped hotspots in this paper
elife2018$HSName[21] <- "GIS"
elife2018$HSName[45] <- "GPA1"
elife2018$HSName[46] <- "ERC1"
elife2018$HSName[47] <- "STB5"
elife2018$HSName[69] <- "HAP1"
elife2018$HSName[77] <- "MOT3"
elife2018$HSName[83] <- "KRE33"
elife2018$HSName[84] <- "MKT1"
elife2018$HSName[88] <- "IRA2"

# Additional fine-mapped hotspots from Lutz et al., PLoS Genetics 2019 (Albert lab)
elife2018$HSName[1]  <- "OAF1"
elife2018$HSName[18] <- "RGT2"
elife2018$HSName[40] <- "OLE1"


# Print this eQTL hotspot table as reference
write_csv(elife2018, paste0(resdir, "Albert2018_hotspots.csv"))


# Generate tibble with all genes overlapping any hotspot (long format)
eHSGenes <- as_tibble(str_split(elife2018$allGenesInInterval, pattern = ";", simplify = T), .name_repair = "universal") %>%
        mutate(eHS.ID      = elife2018$HSID,
               eHS.name    = elife2018$HSName,
               eHS.genes   = elife2018$allGenesInInterval,
               eHS.targetn = elife2018$numberEQTLInHotspot,
               eHS.strVar  = elife2018$genesWithHighImpactVariants,
               eHS.medVar  = elife2018$genesWithModerateImpactVariants) %>%
        dplyr::select(eHS.ID, eHS.name, everything()) %>%
        pivot_longer(cols = -c(eHS.ID, eHS.name, eHS.genes, eHS.targetn, eHS.strVar, eHS.medVar), 
                     names_to = "eHS.pos", values_to = "eHS.gene") %>%
        mutate(eHS.pos = str_replace(eHS.pos, "...", "P")) %>%
        dplyr::filter(eHS.gene != "") %>% print()


# Combine our data with Albert et al., 2018
comb <- combdf.gnx %>%
        full_join(eHSGenes, by = c("Name" = "eHS.gene")) %>%
        dplyr::filter(!is.na(eHS.ID)) %>%
        group_by(eHS.ID) %>%
        mutate(eHS.BE.sig.prots = paste(FDR0.05_prots[!is.na(FDR0.05_prots)], collapse = ";")) %>%
        ungroup() %>% print()


# Combine significant proteins per hotspot (sort and remove duplicates)
comb$eHS.BE.sig.prots <- apply(as.data.frame(unlist(str_split(comb$eHS.BE.sig.prots, pattern = ";", simplify = T))), 
                               1, function(x) paste(sort(unique(x)), collapse = ";"))


# Add back into comb table
comb <- comb %>%
        mutate(eHS.BE.sig.n = str_count(eHS.BE.sig.prots, ";"),
               eHS.BE.sig.prots = ifelse(str_detect(eHS.BE.sig.prots, "^;.+"), 
                                         str_sub(eHS.BE.sig.prots, start = 2, end = -1),
                                         eHS.BE.sig.prots),
               eHS.BE.sig.prots = ifelse(eHS.BE.sig.prots == "", NA_character_, eHS.BE.sig.prots))




# Use BE data to fine-map eQTL hotspots ----------------------------------------

# Annotate candidates to fine-map eQTL hotspots
cand <- comb %>%
        arrange(eHS.ID) %>%
        group_by(eHS.ID) %>%
        mutate(Max_FDR0.05_count = max(FDR0.05_count, na.rm = T),
               Max_FDR0.05_count = ifelse(is.infinite(Max_FDR0.05_count), NA_integer_, Max_FDR0.05_count),
               Candidates        = as.character(ifelse(FDR0.05_count > 2, Name, NA)), 
               ImpactVarMed      = map2_lgl(eHS.medVar, Candidates, str_detect),
               ImpactVarStr      = map2_lgl(eHS.strVar, Candidates, str_detect),
               ImpactVarInCand   = ifelse(ImpactVarMed == T | ImpactVarStr == T, T, F),
               ImpactVarInCand   = replace_na(ImpactVarInCand, FALSE),
               ImpactVarInCand   = sum(ImpactVarInCand),
               Candidates        = paste(Candidates[!is.na(Candidates)], collapse = ";"),
               Candidates        = ifelse(Candidates == "", NA_character_, Candidates)) %>%
        dplyr::select(-c(Name, geneSys, chr, gstart, FDR0.05_count, FDR0.05_prots, 
                         eHS.pos, ImpactVarMed, ImpactVarStr)) %>% 
        ungroup() %>%
        unique() %>%
        mutate(Cat = case_when(!is.na(eHS.name) & str_detect(Candidates, eHS.name) ~ "Both",
                               !is.na(eHS.name) & !str_detect(Candidates, eHS.name) ~ "Diff",
                               !is.na(eHS.name) ~ "Known",
                               !is.na(Candidates) ~ "New",
                               TRUE ~ NA_character_)) %>% 
        mutate(Cat = factor(Cat, levels = c("Known", "New", "Both", "Diff"))) %>%
        mutate(Label = case_when(!is.na(Candidates) ~ Candidates,
                                 !is.na(eHS.name) ~ eHS.name,
                                 TRUE ~ NA_character_)) %>%
        mutate(PlotLabel = str_replace_all(Label, ";", "\n")) %>% 
        arrange(eHS.ID) %>% print()
write_csv(cand, paste0(resdir, "Albert2018_finemapping.csv"))




# Check if eQTL hotspots affect my proteins ------------------------------------

# Import table containing all eQTLs
eQTL <- importExtData(dataset = "Albert2018_eQTLs", localfile = T) %>% 
        dplyr::rename(geneSys = gene) %>%
        left_join(sgd, by = "geneSys") %>%
        mutate(eQ.ID = seq_along(row_number())) %>%
        select(geneSys, gene, eQ.ID, chr, CI.l, CI.r, r, LOD) %>%
        mutate(CI.l = as.numeric(str_extract(CI.l, "(?<=:).+(?=_)")),
               CI.r = as.numeric(str_extract(CI.r, "(?<=:).+(?=_)"))) %>%
        print()


# Convert eQTL table into GRanges object
eQTL.gr <- GRanges(seqnames = eQTL$chr, 
                   ranges = IRanges(eQTL$CI.l, 
                                    eQTL$CI.r, 
                                    names = eQTL$eQ.ID),
                   strand = "*",
                   eQ.ID  = eQTL$eQ.ID); eQTL.gr


# Convert hotspot table into GRanges object
eHSInt.gr <- GRanges(seqnames = elife2018$chromosome, 
                     ranges = IRanges(elife2018$bootstrapIntervalLeft, 
                                      elife2018$bootstrapIntervalRight, 
                                      names = elife2018$HSID),
                     strand = "*",
                     eHS.ID   = elife2018$HSID); eHSInt.gr


# Get eHS.ID for each eQTL
eHS.assoc <- as_tibble(as.data.frame(findOverlaps(eQTL.gr, eHSInt.gr, type = "any"))) %>%
        dplyr::rename(eQ.ID = queryHits,
                      eHS.ID  = subjectHits) %>%  print()


# Get list of proteins that are affected by each eQTL hotspot
eHS.prots <- eQTL %>%
        left_join(eHS.assoc, by = "eQ.ID") %>%
        dplyr::filter(gene %in% toupper(experiments)) %>%
        arrange(gene) %>% 
        group_by(eHS.ID) %>%
        summarise(eHS.sig.prots = paste(gene, collapse = ";")) %>%
        print()


# Combine with other information about hotspots (cand)
cand.prots <- cand %>% left_join(eHS.prots, by = "eHS.ID") %>% print()
write_csv(cand.prots, paste0(resdir, "Albert2018_finemapping.csv"))




# Find overlap with Nature 2014 pQTL hotspots ----------------------------------

# Manually import pQTL hotspots table (Extended Data Table 2)
# Convert coordinates from SacCer2 to SacCer3 using UCSC Liftover tool.
nature2014 <- tibble(pHS.ID = 1:20, 
                     pHS.chr = c(1, 2, 2, 4, 5, 5, 7, 7, 8, 8, 10, 10, 11, 12, 12, 12, 13, 14, 14, 15), 
                     pHS.cpos_SacCer2 = c(39010, 132948, 397978, 223943, 192064, 371845, 137332, 505871, 103041, 
                                          419747, 142009, 655465, 234462, 238302, 656893, 1039502, 96832, 232509, 
                                          465007, 162766),
                     pHS.cpos = c(39009, 132945, 397984, 223943, 192065, 371849, 137326, 505867, 103046, 
                                  419744, 142012, 655474, 234818, 238301, 656891, 1039504, 96832, 232508, 
                                  465005, 162767),
                     pHS.name = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                  NA, NA, "HAP4", NA, "HAP1", NA, NA, NA, "MKT1", "IRA2"),
                     pHS.targetf = c(0.31, 0.31, 0.09, 0.12, 0.16, 0.16, 0.15, 0.16, 0.19, 0.08, 
                                     0.18, 0.11, 0.16, 0.16, 0.41, 0.12, 0.31, 0.13, 0.58, 0.58)) %>% 
        left_join(sgdchr, by = c("pHS.chr" = "chromosome")) %>%
        mutate(pHS.gpos = pHS.cpos + chrlencum) %>%
        dplyr::select(pHS.ID, pHS.chr, pHS.cpos, pHS.gpos, pHS.name, pHS.targetf) %>% print()


# Import pQTL table and keep coordinates as chromosome coordinates
pQTL <- importExtData(dataset = "Albert2014_pQTLs_SacCer3", localfile = T) %>%
        dplyr::rename(geneSys = gene,
                      pQ.chrR  = chromosome,
                      pQ.peak  = peakPosition,
                      pQ.LOD   = LOD,
                      pQ.left  = X2LODIntervalLeft,
                      pQ.right = X2LODIntervalRight,
                      pQ.diff  = alleleFrequencyDifference) %>% 
        mutate(pQ.chr   = as.numeric(as.roman(str_extract(pQ.chrR, "(?<=chr).+"))),
               pQ.peak  = as.numeric(pQ.peak),
               pQ.left  = as.numeric(pQ.left),
               pQ.right = as.numeric(pQ.right)) %>%
        left_join(sgd, by = "geneSys") %>%
        select(gene, geneSys, pQ.chr, pQ.peak, pQ.left, pQ.right, pQ.LOD, pQ.diff) %>%
        arrange(pQ.chr, pQ.peak) %>%
        mutate(pQ.ID = seq_along(row_number())) %>% 
        relocate(gene, geneSys, pQ.ID, everything()) %>% print()


# Build a pQTL GenomicRanges object with chromosome coordinates
pQTL.gr <- GRanges(seqnames = pQTL$pQ.chr, 
                   ranges = IRanges(pQTL$pQ.left, 
                                    pQTL$pQ.right, 
                                    names = pQTL$pQ.ID),
                   strand = "*",
                   pQ.ID  = pQTL$pQ.ID); pQTL.gr


# Find hotspot regions where ≥12 pQTLs overlap (definition by Albert et al., 2014)
pHSInt.irl <- slice(coverage(pQTL.gr), lower = 12, rangesOnly = T); pHSInt.irl


# Merge two ranges on chromosome 8, which are very close
# (Otherwise the hotspot numbering won't match the paper.)
pHSInt.irl[[8]]
end(pHSInt.irl[[8]][1]) <- end(pHSInt.irl[[8]][2])
pHSInt.irl[[8]][2] <- NULL


# Convert into GRanges object
pHSInt.gr <- as(pHSInt.irl, "GRanges"); pHSInt.gr
pHSInt.gr$pHS.ID <- c(1:20)


# Convert into tibble
pHSInt.df <- as_tibble(as.data.frame(pHSInt.gr)) %>% 
        mutate(pHS.chr = as.integer(seqnames)) %>%
        dplyr::rename(pHS.start = start,
                      pHS.end = end) %>%
        select(pHS.ID, pHS.chr, pHS.start, pHS.end) %>% print()


# Get all genes that overlap the hotspot intervals
pHSGenes.df        <- as.data.frame(subsetByOverlaps(txdb, pHSInt.gr, type = "any"))
pHSGenes.df$pHS.ID <- as.data.frame(findOverlaps(txdb, pHSInt.gr, type = "any"))$subjectHits
head(pHSGenes.df)


# Convert to tibble, change names and add RM variant info (long format)
pHSGenes.tb <- as_tibble(pHSGenes.df) %>%
        left_join(sgd, by = c("gene_id" = "geneSys")) %>% 
        dplyr::rename(pHS.gene = gene,
                      pHS.geneSys = gene_id) %>%
        select(pHS.ID, pHS.gene, pHS.geneSys) %>%
        mutate(pHS.strVar.ind = ifelse(pHS.gene %in% RMvar[RMvar$IMPACT == "HIGH", ]$Name, pHS.gene, NA_character_),
               pHS.medVar.ind = ifelse(pHS.gene %in% RMvar[RMvar$IMPACT == "MODERATE", ]$Name, pHS.gene, NA_character_)) %>% print()


# Condense into hotspot table
pHS <- pHSGenes.tb %>%
        group_by(pHS.ID) %>%
        summarise(pHS.genes  = paste(pHS.gene, collapse = ";"),
                  pHS.strVar = paste(pHS.strVar.ind[!is.na(pHS.strVar.ind)], collapse = ";"),
                  pHS.medVar = paste(pHS.medVar.ind[!is.na(pHS.medVar.ind)], collapse = ";")) %>% 
        left_join(nature2014, by = "pHS.ID") %>% 
        left_join(pHSInt.df, by = c("pHS.ID", "pHS.chr")) %>%
        print()
write_csv(pHS, paste0(resdir, "Albert2014_hotspots.csv"))


# Combine the two tables above (long format)
pHSGenes <- pHSGenes.tb %>% 
        left_join(pHS, by = "pHS.ID") %>% 
        arrange(pHS.gpos) %>%
        select(-c(pHS.chr, pHS.cpos, pHS.gpos)) %>% 
        group_by(pHS.ID) %>%
        mutate(pHS.pos = row_number()) %>% 
        select(pHS.ID, pHS.name, pHS.genes, pHS.targetf, 
               pHS.strVar, pHS.medVar, pHS.pos, pHS.gene) %>% print()


# Combine the pQTL and BE hotspot data
combP <- combdf.gnx %>%
        full_join(pHSGenes, by = c("Name" = "pHS.gene")) %>%
        dplyr::filter(!is.na(pHS.ID)) %>%
        group_by(pHS.ID) %>%
        mutate(pHS.BE.sig.prots = paste(FDR0.05_prots[!is.na(FDR0.05_prots)], collapse = ";")) %>%
        ungroup() %>% print()


# Combine sig proteins per hotspot (sort and remove duplicates)
combP$pHS.BE.sig.prots <- apply(as.data.frame(unlist(str_split(combP$pHS.BE.sig.prots, pattern = ";", simplify = T))), 
                                1, function(x) paste(sort(unique(x)), collapse = ";"))

# Add back into comb table
combP <- combP %>%
        mutate(pHS.BE.sig.n = str_count(pHS.BE.sig.prots, ";"),
               pHS.BE.sig.prots = ifelse(str_detect(pHS.BE.sig.prots, "^;.+"), 
                                         str_sub(pHS.BE.sig.prots, start = 2, end = -1),
                                         pHS.BE.sig.prots),
               pHS.BE.sig.prots = ifelse(pHS.BE.sig.prots == "", NA_character_, pHS.BE.sig.prots))


# Create reduced table with one row per pQTL hotspot
combP.red <- combP %>%
        dplyr::select(-c(Name, geneSys, chr, gstart, FDR0.05_count, FDR0.05_prots, pHS.pos)) %>%
        arrange(pHS.ID) %>%
        unique() %>% print()




# Use BE data to fine-map pQTL hotspots ----------------------------------------

# Annotate candidates to fine-map pQTL hotspots
candP <- combP %>%
        arrange(pHS.ID) %>%
        group_by(pHS.ID) %>%
        mutate(Max_FDR0.05_count = max(FDR0.05_count, na.rm = TRUE),
               Max_FDR0.05_count = ifelse(is.infinite(Max_FDR0.05_count), NA_integer_, Max_FDR0.05_count),
               Candidates        = if_else(FDR0.05_count > 2, Name, NA_character_), 
               ImpactVarMed      = map2_lgl(pHS.medVar, Candidates, str_detect),
               ImpactVarStr      = map2_lgl(pHS.strVar, Candidates, str_detect),
               ImpactVarInCand   = ifelse(ImpactVarMed == T | ImpactVarStr == T, T, F),
               ImpactVarInCand   = replace_na(ImpactVarInCand, FALSE),
               ImpactVarInCand   = sum(ImpactVarInCand),
               Candidates        = paste(Candidates[!is.na(Candidates)], collapse = ";"),
               Candidates        = ifelse(Candidates == "", NA_character_, Candidates)) %>%
        dplyr::select(-c(Name, geneSys, chr, gstart, FDR0.05_count, FDR0.05_prots,
                         pHS.pos, ImpactVarMed, ImpactVarStr)) %>%
        ungroup() %>%
        unique() %>%
        mutate(Cat = case_when(!is.na(pHS.name) & str_detect(Candidates, pHS.name) ~ "Both",
                               !is.na(pHS.name) & !str_detect(Candidates, pHS.name) ~ "Diff",
                               !is.na(pHS.name) ~ "Known",
                               !is.na(Candidates) ~ "New",
                               TRUE ~ NA_character_)) %>% 
        mutate(Cat = factor(Cat, levels = c("Known", "New", "Both", "Diff"))) %>%
        mutate(Label = case_when(!is.na(Candidates) ~ Candidates,
                                 !is.na(pHS.name) ~ pHS.name,
                                 TRUE ~ NA_character_)) %>%
        mutate(PlotLabel = str_replace_all(Label, ";", "\n")) %>% 
        arrange(pHS.ID) %>% print()




# Check if pQTL hotspots affect my proteins ------------------------------------

# Get eHS.ID for each eQTL
pHS.assoc <- as_tibble(as.data.frame(findOverlaps(pQTL.gr, pHSInt.gr, type = "any"))) %>%
        dplyr::rename(pQ.ID = queryHits,
                      pHS.ID  = subjectHits) %>%  print()


# Get list of proteins that are affected by each eQTL hotspot
pHS.prots <- pQTL %>%
        left_join(pHS.assoc, by = "pQ.ID") %>%
        dplyr::filter(gene %in% toupper(experiments)) %>%
        arrange(gene) %>% 
        group_by(pHS.ID) %>%
        summarise(pHS.sig.prots = paste(gene, collapse = ";")) %>%
        print()


# Combine with other information about hotspots (cand)
candP.prots <- candP %>% left_join(pHS.prots, by = "pHS.ID") %>% print()
write_csv(candP.prots, paste0(resdir, "Albert2014_finemapping.csv"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-15_SessionInfo.txt")

