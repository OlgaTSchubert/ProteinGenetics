# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/abundanceEffects/"
dir.create(resdir)


# Import data
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()


# Import protein abundances 

# Yeast GFP Collection (Newman et al., Nature 2006)
newman <- importExtData(dataset = "Newman2006", localfile = T) %>% select(c(1:3))
colnames(newman) <- c("geneSys", "gene", "abundance"); newman

# Yeast proteomics (Ho et al., Cell Systems 2018)
ho <- importExtData(dataset = "Ho2018", localfile = T) %>% select(c(1, 2, 4))
colnames(ho) <- c("geneSys", "gene", "abundance"); ho




# Plot number of significant genes vs GFP-protein abundance --------------------

# Number of significant genes per protein
gn <- combdf.gn %>%
        pivot_longer(cols = -c(geneSys, gene, geneDescr, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>% 
        filter(q < 0.05) %>% 
        group_by(protein) %>% 
        summarize(n = n()) %>% 
        print()


# Newman et al.
newman %>%
        arrange(gene) %>%
        mutate(gene = str_to_title(gene),
               gene = ifelse(gene == "Rpl9a", "Rpl9A", gene)) %>%
        filter(gene %in% gn$protein) %>%
        select(gene, abundance) %>%
        mutate(abundance = log10(abundance)) %>%
        left_join(gn, by = c("gene" = "protein")) %>%
        ggplot(aes(x = abundance, y = n)) +
        geom_smooth(method = "lm") +
        geom_point() +
        geom_text(aes(label = gene), hjust = 0, nudge_x = 0.05) +
        ggpubr::stat_cor(label.x.npc = "center", label.y.npc = "top") +
        scale_x_continuous(limits = c(3.2, 5.2)) +
        xlab("Log10(abundance) by Newman et al., 2006") +
        ylab("Number of genes") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "NsigVsAbu_Newman2006.pdf"), width = 5, height = 4)


# Ho et al.
ho %>%
        arrange(gene) %>%
        mutate(gene = str_to_title(gene),
               gene = ifelse(gene == "Rpl9a", "Rpl9A", gene)) %>%
        filter(gene %in% gn$protein) %>%
        select(gene, abundance) %>%
        mutate(abundance = log10(abundance)) %>%
        left_join(gn, by = c("gene" = "protein")) %>%
        ggplot(aes(x = abundance, y = n)) +
        geom_smooth(method = "lm") +
        geom_point() +
        geom_text(aes(label = gene), hjust = 0, nudge_x = 0.05) +
        ggpubr::stat_cor(label.x.npc = "center", label.y.npc = "top") +
        scale_x_continuous(limits = c(4.6, 6.4)) +
        xlab("Log10(abundance) by Ho et al., 2018") +
        ylab("Number of genes") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "NsigVsAbu_Ho2018.pdf"), width = 5, height = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-09_SessionInfo.txt")

