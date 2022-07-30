# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/geneticIntComparison/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()

schubert <- combdf.gn %>%
        #filter(!is.na(gene)) %>%   # remove genes YGR109W-B and YIL082W-A
        #select(-geneSys, -geneDescr, -essential, -FDR0.05_count) %>%
        pivot_longer(-c(gene, geneSys, geneDescr, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        mutate(protein = toupper(protein)) %>%
        print()


# Get list of GFP proteins
prots    <- toupper(str_extract(str_subset(names(combdf.gn), "_log2fc"), ".*(?=_log2fc)")) %>% print()
protsSys <- importExtData(dataset = "SGD_features", localfile = T) %>%
        filter(Standard_gene_name %in% prots) %>%
        select("Feature_name") %>%
        pull() %>% print()


# Import genetic interaction tables from The Cell Map (https://thecellmap.org/)
gi.o <- importExtData(dataset = "Costanzo2016", localfile = T) %>% print()




# Reformat tables, get counts and combine -------------------------------------- 

gi <- gi.o %>%
        select(-direction) %>%
        # Better delete FAS2 entries because it was tested for many different alleles
        mutate(protein = ifelse(str_detect(protein, "fas2"), "FAS2", protein)) %>%
        filter(protein != "FAS2") %>%
        filter(p < 0.05) %>% 
        group_by(protein, ORF, allele) %>%
        slice(ifelse(score > 0, which.max(score), which.min(score))) %>%
        print()

gisperprot <- gi %>%
        filter(abs(score) > 0.08) %>%  # "intermediate filter"
        #filter(score > 0.16 | score < -0.12) %>%  # "stringent filter"
        group_by(protein) %>%
        summarize(n = n()) %>%
        print()

comb <- schubert %>%
        filter(q < 0.05) %>%
        group_by(protein) %>%
        summarize(n = n()) %>%
        left_join(gisperprot, suffix = c("schubert", "costanzo"), by = "protein") %>%
        print()




# Scatter plot ----------------------------------------------------------------- 

comb %>%
        ggplot(aes(x = nschubert, y = ncostanzo)) +
        geom_smooth(method = "lm") +
        geom_point() +
        ggpubr::stat_cor(color = "grey25",
                         label.x.npc = "left", label.y.npc = "top") +
        ggrepel::geom_text_repel(label = comb$protein,
                size = 2, box.padding = 0.2, point.padding = 0.05,
                segment.color = "darkgrey", segment.size = 0.5) +
        labs(x = "This study", y = "Costanzo et al., 2016") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "geneticIntPerProt.pdf"), width = 4, height = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-17_SessionInfo.txt")

