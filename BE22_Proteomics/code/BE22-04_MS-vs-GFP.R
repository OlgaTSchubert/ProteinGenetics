# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Functions
source("../BE00_gRNALibraryDesign/code/BE00-00_ImportExtData.R")


# Results directory
resdir <- "results/MS-vs-GFP/"
dir.create(resdir)


# Import annotations
muts  <- read_delim("annotations/MutationLookup.tsv", delim = "\t") %>% as_tibble %>% print()
prots <- read_delim("annotations/ProteinLookup.tsv", delim = "\t") %>% pull(gene) %>% toupper() %>% print()


# Import results from base editor screens (guide-level, not gene-level)
gfp <- readRDS("../BE19_GFPScreens/results/processing/combdf_gd.RDS") %>% 
        select(guide, gene, Eno2_log2fc:Yhb1_q) %>%
        filter(guide %in% muts$guide) %>%
        pivot_longer(-c(guide, gene), 
                     names_to = c("protein", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        mutate(protein = toupper(protein)) %>%
        rename(GFP_log2fc = log2fc,
               GFP_q = q) %>% print()


# Import and reformat artMS results
ms <- read_delim("results/artMS/artMS_results_OS.tsv", delim = "\t") %>%
        mutate(adj.p.i = ifelse(is.na(adj.p), 1E-10, adj.p.i)) %>%
        filter(gene %in% prots) %>%
        rename(protein    = gene,
               MS_log2fc  = log2FC,
               MS_q       = adj.p,
               MS_log2fci = log2FC.i,
               MS_qi      = adj.p.i) %>% 
        select(-geneSys) %>% print()


# Import gene annotations
sgd <- importExtData(dataset = "SGD_features", localfile = T) %>%
        dplyr::rename(geneSys = Feature_name, 
                      gene = Standard_gene_name) %>% 
        mutate(gene = ifelse(is.na(gene), geneSys, gene)) %>% 
        select(geneSys, gene) %>% print()




# Combine data -----------------------------------------------------------------
# There is no difference between data with and without imputation

comb <- muts %>%
        left_join(gfp, by = c("guide", "gene")) %>%
        left_join(ms, by = c("alt" = "mutant", "protein")) %>%
        select(-MS_log2fci, -MS_qi) %>%   # no difference, see below
        mutate(protein = str_to_title(protein)) %>%
        print()




# Scatter plot -----------------------------------------------------------------

# All proteins in POP1
comb %>%
        filter(gene %in% c("POP1")) %>%
        mutate(sig = ifelse(GFP_q < 0.1 & MS_q < 0.1, T, F)) %>%
        ggplot(aes(x = GFP_log2fc, y = MS_log2fc, color = sig)) +
        geom_hline(aes(yintercept = 0), color = "grey") +
        geom_vline(aes(xintercept = 0), color = "grey") +
        geom_point() +
        ggrepel::geom_text_repel(aes(label = protein, color = sig),
                                 size = 4, box.padding = 0.3, point.padding = 0.1) +
        scale_color_manual(values = c("darkgrey", "tomato")) +
        scale_x_continuous(limits = c(-2, 1), name = "GFP Log2FC") +
        scale_y_continuous(limits = c(-2, 2), name = "MS Log2FC") +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "GFP-vs-MS_scatter_POP1.pdf"), width = 3, height = 3)




# Heatmaps ---------------------------------------------------------------------

# For Tdh123 in IRA2* and RAS2*
comb %>%
        pivot_longer(-c(guide, geneSys, gene, mutation, alt, protein),
                     names_to = c("type", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        mutate(alt = factor(alt, levels = c("IRA2", "RAS22"))) %>%
        filter(!is.na(type)) %>%
        filter(!is.na(alt)) %>%
        filter(alt %in% c("IRA2", "RAS22")) %>%
        filter(protein %in% c("Tdh1", "Tdh2")) %>%
        mutate(alt = recode(alt,
                            "IRA2"  = "IRA2 W342*",
                            "RAS22" = "RAS2 Q272*")) %>%
        mutate(log2fc = case_when(log2fc >  1 ~  1,
                                  log2fc < -1 ~ -1,
                                  TRUE ~ log2fc)) %>%
        ggplot() +
        geom_tile(aes(x = protein, y = alt, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1, 1)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_grid(~ type)
ggsave(paste0(resdir, "GFP-vs-MS_heatmap_Tdh12.pdf"), width = 3, height = 1.55)


# For Tdh123 in IRA2*, RAS2-S46L and RAS2* (GFP only)
comb %>%
        pivot_longer(-c(guide, geneSys, gene, mutation, alt, protein),
                     names_to = c("type", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        mutate(alt = factor(alt, levels = c("IRA2", "RAS21", "RAS22"))) %>%
        filter(type == "GFP") %>% 
        filter(!is.na(alt)) %>%
        filter(alt %in% c("IRA2", "RAS21", "RAS22")) %>%
        filter(protein %in% c("Tdh1", "Tdh2", "Tdh3")) %>% 
        mutate(alt = recode(alt, 
                            "IRA2"  = "IRA2 W342*",
                            "RAS22" = "RAS2 Q272*",
                            "RAS21" = "RAS2 S46L")) %>%
        mutate(log2fc = case_when(log2fc >  1 ~  1,
                                  log2fc < -1 ~ -1,
                                  TRUE ~ log2fc)) %>% #print()
        ggplot() +
        geom_tile(aes(x = protein, y = alt, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1, 1)) +
        scale_y_discrete(limits = rev) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(resdir, "GFP_heatmap_discordant.pdf"), width = 3, height = 1.8)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-04_SessionInfo.txt")

