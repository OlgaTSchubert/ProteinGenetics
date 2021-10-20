# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/broadRegulators/"
dir.create(resdir)


# Import results
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()


# Convert results to long format
gn <- combdf.gn %>%
        select(-geneDescr) %>%
        pivot_longer(cols = -c(geneSys, gene, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(!is.na(log2fc)) %>%
        print()




# Pie chart for fraction of genes by specificity -------------------------------

pie <- combdf.gn %>%
        filter(FDR0.05_count > 0) %>%
        mutate(FDR0.05_count = case_when(FDR0.05_count < 8 ~ as.character(FDR0.05_count),
                                         FDR0.05_count > 7 ~ "8+",
                                         TRUE ~ NA_character_)) %>%
        group_by(FDR0.05_count) %>%
        summarize(n = n()) %>%
        arrange(desc(FDR0.05_count)) %>%
        mutate(prop = n / sum(n) * 100) %>%
        mutate(ypos = cumsum(prop) - 0.5*prop ) %>% print()
write_csv(pie, paste0(resdir, "Piechart.csv"))

pie %>% ggplot(aes(x = "", y = prop, fill = FDR0.05_count)) +
        geom_bar(stat = "identity", color = "black", size = 0.3) +
        coord_polar("y", start = 0) +
        scale_fill_grey(start = 0.99, end = 0.5) +
        theme_void() +
        theme(legend.position = "none") +
        geom_text(aes(y = ypos, x = 1.3, label = FDR0.05_count), size = 4, color = "black")
ggsave(paste0(resdir, "Piechart.pdf"), width = 2, height = 2)




# Heatmap ----------------------------------------------------------------------

gn %>%
        filter(FDR0.05_count > 7) %>%
        #filter(FDR0.05_count > 7 | gene == "SAP155") %>%
        #filter(gene %in% c("POP1", "SIT4", "SAP155")) %>%
        #mutate(gene = factor(gene, levels = c("POP1", "SIT4", "SAP155"))) %>%
        mutate(log2fc = case_when(log2fc >  1.5 ~  1.5,
                                  log2fc < -1.5 ~ -1.5,
                                  TRUE ~ log2fc)) %>%
        ggplot() +
        geom_tile(aes(x = protein, y = gene, fill = log2fc)) +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC", limit = c(-1.5, 1.5)) +
        scale_y_discrete(limits = rev) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "Heatmap_29.pdf"), width = 3.3, height = 5)
ggsave(paste0(resdir, "Heatmap_29_SAP155.pdf"), width = 3.3, height = 5)
ggsave(paste0(resdir, "Heatmap_POP1_SIT4_SAP155.pdf"), width = 3.3, height = 2)




# Fraction of ess/set as function of specificity -------------------------------

# By set
combdf.gd %>%
        mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
        mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
        mutate(set = forcats::fct_rev(as.factor(set))) %>%
        ggplot(aes(x = FDR0.05_count, fill = set)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("lightsalmon2", "lightsalmon3", "lightsalmon4")) +
        labs(x = "Number of proteins affected", y = "Fraction of gRNAs") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "Fraction_set_spec_gd.pdf"), height = 3, width = 4)


# By essentiality (guide-level)
combdf.gd %>%
        #filter(set %in% c("neProvs", "eProvs")) %>%
        mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
        mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
        mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
               essential = forcats::fct_rev(as.factor(essential))) %>%
        ggplot(aes(x = FDR0.05_count, fill = essential)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
        labs(x = "Number of proteins affected", y = "Fraction of gRNAs") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "Fraction_ess_spec_gd.pdf"), height = 3, width = 4)
ggsave(paste0(resdir, "Fraction_ess_spec_gd_ProvsOnly.pdf"), height = 3, width = 4)


# By essentiality (gene-level)
combdf.gn %>%
        mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
        mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
        mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
               essential = forcats::fct_rev(as.factor(essential))) %>%
        ggplot(aes(x = FDR0.05_count, fill = essential)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
        labs(x = "Number of proteins affected", y = "Fraction of genes") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "Fraction_ess_spec_gn.pdf"), height = 3, width = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-14_SessionInfo.txt")


