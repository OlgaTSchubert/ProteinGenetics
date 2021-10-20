# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/positionEffects/"
dir.create(resdir)


# Import results
combdf.gd   <- readRDS("results/processing/combdf_gd.RDS") %>% print()

cols_log2fc <- str_subset(names(combdf.gd), "_log2fc")
cols_qval   <- str_subset(names(combdf.gd), "_q")
experiments <- str_sub(cols_qval, start = 1, end = -3) %>% print()


# Import table containing protein lengths
sgdlen <- importExtData(dataset = "SGD_proteins", localfile = T) %>% 
        select(ORF, `Protein Length`) %>% 
        rename(geneSys = ORF, protlen = `Protein Length`) %>% print()


# Generate table with protein length etc. for each mutation
mutations <- left_join(combdf.gd, sgdlen, by = "geneSys") %>%
        mutate(mutLoc = as.integer(str_extract(.$mut1, "\\d+"))) %>%
        mutate(mutLocFromEnd = protlen - mutLoc) %>%
        mutate(mutLocRel = mutLoc/protlen) %>% 
        select(guide, set, geneSys, gene, protlen, mutLoc, mutLocFromEnd,
               mutLocRel, FDR0.05_count, cols_log2fc, cols_qval) %>% print()
        



# Plot dependence of number gene perturbations on position ---------------------

mutations %>% 
        mutate(sig = ifelse(FDR0.05_count > 0, T, F)) %>%
        filter(mutLocFromEnd < 220) %>% 
        mutate(bin = cut_width(mutLocFromEnd, width = 20, boundary = 0)) %>%
        group_by(bin, sig) %>%
        summarize(n = n()) %>%
        mutate(freq = n/sum(n)) %>%
        filter(sig == T) %>% 
        ungroup() %>% 
        arrange(-freq) %>%
        ggplot(aes(x = bin, y = freq)) +
        geom_bar(stat = "identity", position = position_nudge(-0.5)) +
        scale_x_discrete(limits = rev, labels = as.character(seq(200, 0, by = -20))) +
        expand_limits(x = c(-0.5, 10)) +
        xlab("Distance from 3' end of gene (codons)") +
        ylab("Ratio of significant vs all gRNAs") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "NsigRatio_BinsFromEnd.pdf"), width = 4, height = 4)




# Plot dependence of effect size on position -----------------------------------

# Distance from 3' end of gene (codons)
mutations %>%
        pivot_longer(cols = -c(guide, set, geneSys, gene, protlen, 
                               mutLoc, mutLocFromEnd, mutLocRel, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(q < 0.05) %>%
        filter(mutLocFromEnd < 200) %>% 
        mutate(Stop = ifelse(set == "neStops", "neStops", "neProvs / eProvs")) %>% #print()
        ggplot(aes(x = mutLocFromEnd, y = abs(log2fc))) +
        geom_point() +
        geom_smooth(method = "loess") +
        scale_x_reverse() +
        ylim(0, 2) +
        xlab("Distance from 3' end of gene (codons)") +
        ylab("Abs(log2FC)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none")
ggsave(paste0(resdir, "effectSize_scatter_loess_fromEnd.pdf"), width = 4, height = 4)


# Relative gene position
mutations %>%
        pivot_longer(cols = -c(guide, set, geneSys, gene, protlen, mutLoc, mutLocFromEnd, mutLocRel, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(q < 0.05) %>%
        mutate(Stop = ifelse(set == "neStops", "neStops", "neProvs / eProvs")) %>% #print()
        ggplot(aes(x = mutLocRel, y = abs(log2fc))) +
        geom_point() +
        geom_smooth(method = "loess") +
        ylim(0, 2) +
        xlab("Relative position in gene") +
        ylab("Abs(log2FC)") +
        #facet_wrap(~ Stop) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none")
ggsave(paste0(resdir, "effectSize_scatter_loess.pdf"), width = 4, height = 4)
ggsave(paste0(resdir, "effectSize_scatter_loess_facet.pdf"), width = 5, height = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-10_SessionInfo.txt")

