# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/specificRegulators/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()


# Labels for plots
volclabs <- read_delim("annotations/VolcanoLabels.tsv", delim = "\t") %>% print()
mediator <- read_delim("annotations/Mediator.tsv", delim = "\t") %>% print(n=30)
sps      <- read_delim("annotations/SPS.tsv", delim = "\t") %>% print()


# Convert results to long format
gn <- combdf.gn %>%
        select(-geneDescr) %>%
        pivot_longer(cols = -c(geneSys, gene, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(!is.na(log2fc)) %>%
        print()




# Volcano plots ----------------------------------------------------------------

gn %>%
        left_join(volclabs, by = c("gene", "protein")) %>%
        mutate(negLogQ = -log10(q),
               negLogQ = ifelse(negLogQ > 30, 30, negLogQ)) %>%
        arrange(-negLogQ) %>%
        ggplot(aes(x = log2fc, y = negLogQ, color = comment)) +
        geom_point(data = . %>% filter(is.na(set)), color = "gray80") +
        geom_point(data = . %>% filter(q < 0.05 & FDR0.05_count > 7), color = "gray40") +
        geom_point(data = . %>% filter(!is.na(set))) +
        ylab("-Log10(q)") +
        xlab("Log2FC") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              legend.position = "none") +
        ggrepel::geom_text_repel(
                data = . %>% filter(!is.na(set)) %>% distinct(comment, protein, .keep_all = TRUE),
                mapping = aes(x = log2fc, y = negLogQ, label = comment, color = comment),
                size = 3,
                box.padding = 0.3, point.padding = 0.1, force = 50,
                segment.color = "darkgrey", segment.size = 0.3) +
        facet_wrap(~ protein, nrow = 2)
ggsave(paste0(resdir, "VolcanoPlots.pdf"), width = 10, height = 5)




# Heatmap for Mediator ---------------------------------------------------------

gn %>%
        inner_join(mediator, by = c("gene", "geneSys")) %>%
        mutate(gene    = factor(gene2, levels = rev(pull(mediator, gene2))),
               set     = factor(set, levels = c("Cdk8", "Middle", "Head", "Tail")),
               protein = factor(protein, levels = c("Htb2", "Tdh1", "Tdh2", "Tdh3", 
                                                    "Eno2", "Fas1", "Fas2", "Rnr2", 
                                                    "Rpl9A", "Ssa1", "Yhb1"))) %>%
        mutate(log2fc = case_when(log2fc >  1.5 ~  1.5,
                                  log2fc < -1.5 ~ -1.5,
                                  TRUE ~ log2fc)) %>%
        ggplot() +
        geom_tile(aes(x = protein, y = gene, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1.5, 1.5)) +
        facet_grid(vars(set), space = "free", scales = "free_y") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(resdir, "Mediator_heatmap.pdf"), width = 4, height = 3.8)




# Heatmap for SPS sensor -------------------------------------------------------

gn %>%
        inner_join(sps, by = c("gene", "geneSys")) %>% #print(n=50)
        mutate(gene    = factor(gene, levels = rev(pull(sps, gene))),
               protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2", 
                                                    "Htb2", "Rnr2", "Rpl9A", "Ssa1", 
                                                    "Tdh1", "Tdh2", "Tdh3", "Yhb1"))) %>%
        mutate(log2fc = case_when(log2fc >  1.5 ~  1.5,
                                  log2fc < -1.5 ~ -1.5,
                                  TRUE ~ log2fc)) %>%
        ggplot() +
        geom_tile(aes(x = protein, y = gene, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1.5, 1.5)) +
        coord_fixed() +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(resdir, "SPS_heatmap.pdf"), width = 3.3, height = 1.85)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-12_SessionInfo.txt")

