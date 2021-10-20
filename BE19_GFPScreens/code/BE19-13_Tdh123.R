# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/Tdh123/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()


# Convert results to long format
gn <- combdf.gn %>%
        select(-geneDescr) %>%
        pivot_longer(cols = -c(geneSys, gene, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(!is.na(log2fc)) %>%
        print()


# Import annotations
CCM <- importExtData(dataset = "SGD_pathways", localfile = T) %>%
        filter(!is.na(gene_name)) %>%
        filter(biochemical_pathway_common_name %in% c("glycolysis",
                                                      "gluconeogenesis",
                                                      "non-oxidative branch of the pentose phosphate pathway",
                                                      "oxidative branch of the pentose phosphate pathway",
                                                      "TCA cycle, aerobic respiration",
                                                      "glyoxylate cycle",
                                                      "superpathway of glucose fermentation")) %>% print()

Pka      <- read_delim("annotations/PKA.tsv", delim = "\t") %>% print()

Mediator <- read_delim("annotations/Mediator.tsv", delim = "\t") %>% print()
MedCdk8  <- Mediator %>% filter(set == "Cdk8") %>% print()

Diff <- c("ERG20", "TRR1", "ILV3")
Same <- c("ARO1", "ARO2", "ADK1", "QNS1", "REB1")




# Heatmap ----------------------------------------------------------------------

gn %>%
        inner_join(Pka, by = c("gene", "geneSys")) %>% 
        filter(protein %in% c("Tdh1", "Tdh2", "Tdh3")) %>%
        group_by(gene) %>%
        mutate(minq = min(q)) %>%
        ungroup() %>%
        filter(minq < 0.1) %>%  # Keep gene only if sign. for at least 1 protein
        mutate(gene    = factor(gene, levels = rev(Pka$gene)),
               set     = factor(set,  levels = c("Activators", "Inhibitors")),
               protein = factor(protein, levels = c("Tdh1", "Tdh2", "Tdh3"))) %>%
        mutate(log2fc = case_when(log2fc >  1 ~  1,
                                  log2fc < -1 ~ -1,
                                  TRUE ~ log2fc)) %>%
        ggplot() +
        geom_tile(aes(x = protein, y = gene, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1, 1)) +
        facet_grid(vars(set), space = "free", scales = "free_y") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(resdir, "PKA_heatmap.pdf"), width = 2.5, height = 2.5)




# Scatterplots -----------------------------------------------------------------

# Add annotations to data table
tdh123 <- combdf.gn %>%
        select(geneSys, gene, geneDescr, starts_with("Tdh")) %>%
        mutate(type12 = ifelse(Tdh1_log2fc * Tdh2_log2fc < 0, "diff", "same"),
               type13 = ifelse(Tdh1_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               type23 = ifelse(Tdh2_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               label = case_when(gene %in% CCM$gene_name ~ "CCM",
                                 gene %in% Same ~ "Same",
                                 gene %in% Diff ~ "Diff",
                                 gene %in% Pka$gene ~ "RasPKA",
                                 gene %in% MedCdk8$gene ~ "MedCdk8",
                                 TRUE ~ NA_character_),
               label = factor(label, levels = c("Same", "CCM", "Diff", "RasPKA", "MedCdk8"))) %>% print()


colo <- c("royalblue4", "royalblue1", "tomato4", "tomato", "goldenrod2")

tdh123 %>%
        filter(Tdh1_q < 0.1 & Tdh2_q < 0.1) %>%
        ggplot(aes(x = Tdh1_log2fc, y = Tdh2_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "gray") +
        geom_hline(yintercept = 0, color = "gray") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label))) +
        scale_color_manual(values = colo) +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 2.5,
                                 box.padding = 0.3, point.padding = 0.1,
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh1 log2FC") +
        ylab("Tdh2 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh2_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh2_log2fc")], na.rm = T))
ggsave(paste0(resdir, "Tdh1-Tdh2.pdf"), width = 3, height = 3)


tdh123 %>%
        filter(Tdh1_q < 0.1 & Tdh3_q < 0.1) %>%
        ggplot(aes(x = Tdh1_log2fc, y = Tdh3_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "grey") +
        geom_hline(yintercept = 0, color = "grey") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label))) +
        scale_color_manual(values = colo) +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 2.5,
                                 box.padding = 0.3, point.padding = 0.1,
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh1 log2FC") +
        ylab("Tdh3 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh3_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh3_log2fc")], na.rm = T))
ggsave(paste0(resdir, "Tdh1-Tdh3.pdf"), width = 3, height = 3)


tdh123 %>%
        filter(Tdh2_q < 0.1 & Tdh3_q < 0.1) %>%
        ggplot(aes(x = Tdh2_log2fc, y = Tdh3_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "grey") +
        geom_hline(yintercept = 0, color = "grey") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label))) +
        scale_color_manual(values = colo) +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 2.5, max.overlaps = 15, 
                                 box.padding = 0.3, point.padding = 0.1, 
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh2 log2FC") +
        ylab("Tdh3 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh2_log2fc", "Tdh3_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh2_log2fc", "Tdh3_log2fc")], na.rm = T))
ggsave(paste0(resdir, "Tdh2-Tdh3.pdf"), width = 3, height = 3)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-13_SessionInfo.txt")




# Scatter plots for slides (fewer annotations) ---------------------------------

# Slides directory
slidesdir <- "results/Tdh123/slides"
dir.create(slidesdir)


# Prepare tables with few annotations
tdh123_CCM <- combdf.gn %>%
        select(geneSys, gene, geneDescr, starts_with("Tdh")) %>%
        mutate(type12 = ifelse(Tdh1_log2fc * Tdh2_log2fc < 0, "diff", "same"),
               type13 = ifelse(Tdh1_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               type23 = ifelse(Tdh2_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               label = case_when(gene %in% CCM$gene_name ~ "CCM",
                                 TRUE ~ NA_character_),
               label = factor(label, levels = c("CCM"))) %>% print()

tdh123_Med <- combdf.gn %>%
        select(geneSys, gene, geneDescr, starts_with("Tdh")) %>%
        mutate(type12 = ifelse(Tdh1_log2fc * Tdh2_log2fc < 0, "diff", "same"),
               type13 = ifelse(Tdh1_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               type23 = ifelse(Tdh2_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               label = case_when(gene %in% MedCdk8$gene ~ "MedCdk8",
                                 TRUE ~ NA_character_),
               label = factor(label, levels = c("MedCdk8"))) %>% print()

tdh123_PKA <- combdf.gn %>%
        select(geneSys, gene, geneDescr, starts_with("Tdh")) %>%
        mutate(type12 = ifelse(Tdh1_log2fc * Tdh2_log2fc < 0, "diff", "same"),
               type13 = ifelse(Tdh1_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               type23 = ifelse(Tdh2_log2fc * Tdh3_log2fc < 0, "diff", "same"),
               label = case_when(gene %in% Pka$gene ~ "RasPKA",
                                 TRUE ~ NA_character_),
               label = factor(label, levels = c("RasPKA"))) %>% print()


tdh123_CCM %>%
#tdh123_Med %>%
#tdh123_PKA %>%
        filter(Tdh1_q < 0.1 & Tdh2_q < 0.1) %>%
        ggplot(aes(x = Tdh1_log2fc, y = Tdh2_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "gray") +
        geom_hline(yintercept = 0, color = "gray") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label)), color = "tomato") +
        #geom_point(color = "black") +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 4,
                                 box.padding = 0.3, point.padding = 0.1,
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh1 log2FC") +
        ylab("Tdh2 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh2_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh2_log2fc")], na.rm = T))
ggsave(paste0(slidesdir, "Tdh1-Tdh2_black.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh2_CCM.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh2_Med.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh2_PKA.pdf"), width = 3, height = 3)


#tdh123_CCM %>%
tdh123_Med %>%
#tdh123_PKA %>%
        filter(Tdh1_q < 0.1 & Tdh3_q < 0.1) %>%
        ggplot(aes(x = Tdh1_log2fc, y = Tdh3_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "grey") +
        geom_hline(yintercept = 0, color = "grey") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label)), color = "tomato") +
        #geom_point(color = "black") +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 4,
                                 box.padding = 0.3, point.padding = 0.1,
                                 ylim = c(0, 2),
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh1 log2FC") +
        ylab("Tdh3 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh3_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh1_log2fc", "Tdh3_log2fc")], na.rm = T))
ggsave(paste0(slidesdir, "Tdh1-Tdh3_black.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh3_CCM.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh3_Med.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh1-Tdh3_PKA.pdf"), width = 3, height = 3)


tdh123_CCM %>%
#tdh123_Med %>%
#tdh123_PKA %>%
        filter(Tdh2_q < 0.1 & Tdh3_q < 0.1) %>%
        ggplot(aes(x = Tdh2_log2fc, y = Tdh3_log2fc, color = label)) +
        geom_vline(xintercept = 0, color = "grey") +
        geom_hline(yintercept = 0, color = "grey") +
        geom_point(data = . %>% filter(is.na(label)), color = "gray80") +
        geom_point(data = . %>% filter(!is.na(label)), color = "tomato") +
        #geom_point(color = "black") +
        ggrepel::geom_text_repel(data = ~ subset(., !is.na(label)),
                                 aes(label = gene), size = 4, max.overlaps = 15,
                                 box.padding = 0.3, point.padding = 0.1,
                                 xlim = c(0, 2), ylim = c(0, 2),
                                 segment.color = "darkgray", segment.size = 0.3) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none") +
        xlab("Tdh2 log2FC") +
        ylab("Tdh3 log2FC") +
        coord_cartesian(xlim = range(combdf.gn[, c("Tdh2_log2fc", "Tdh3_log2fc")], na.rm = T), 
                        ylim = range(combdf.gn[, c("Tdh2_log2fc", "Tdh3_log2fc")], na.rm = T))
ggsave(paste0(slidesdir, "Tdh2-Tdh3_black.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh2-Tdh3_CCM.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh2-Tdh3_Med.pdf"), width = 3, height = 3)
ggsave(paste0(slidesdir, "Tdh2-Tdh3_PKA.pdf"), width = 3, height = 3)


