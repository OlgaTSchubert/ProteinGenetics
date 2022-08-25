# Libraries --------------------------------------------------------------------

library(tidyverse)
library(ggpubr)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE21_FitnessScreen")
getwd()


# Results directory
resdir <- "results/fitness-vs-protein/"
dir.create(resdir)


# Import fitness data
fitn.gd <- readRDS("results/fitness_all4/all4_Ratios_gd_long.RDS") %>% 
        select(-set, -geneSys, -gene, -dropout) %>%
        pivot_wider(id_cols = guide, names_from = timepoint, values_from = log2fc) %>%
        select(guide, T00h, T08h, T24h) %>% 
        rename(T00h_log2fc = T00h,
               T08h_log2fc = T08h,
               T24h_log2fc = T24h) %>% print()


# Import protein data
prot.gd <- readRDS("../BE19_GFPScreens/results/processing/combdf_gd.RDS") %>%
        select(guide, set, essential, gene, geneSys, Eno2_log2fc:FDR0.05_count) %>% print()


# Combine protein and fitness data
comb.gd <- prot.gd %>%
        left_join(fitn.gd, by = "guide") %>% print()
write_csv(comb.gd, paste0(resdir, "combdf_gd_fit-prot.csv"))




# Plots ------------------------------------------------------------------------

(p1 <- comb.gd %>%
         mutate(FDR0.05_count_bin = str_c("n", FDR0.05_count)) %>%
         mutate(FDR0.05_count_bin = ifelse(FDR0.05_count_bin %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count_bin)) %>%
         mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
                essential = forcats::fct_rev(as.factor(essential))) %>%
         ggplot() +
         geom_point(aes(x = FDR0.05_count_bin, y = T24h_log2fc, color = essential),
                    position = "jitter", alpha = 0.7, show.legend = T) +
         geom_boxplot(aes(x = FDR0.05_count_bin, y = T24h_log2fc),
                      alpha = 0.8, outlier.shape = NA, show.legend = F) +
         ggpubr::stat_cor(aes(x = FDR0.05_count, y = T24h_log2fc),
                          label.x.npc = 0.38, label.y.npc = 0.82) +
         scale_color_manual(values = c("darkseagreen3", "darkseagreen4")) +
         coord_cartesian(ylim = c(-2.5, 1)) +
         labs(x = "Number of proteins affected by gRNA",
              y = "Fitness effect of gRNA (Log2FC)",
              color = "gRNA target gene") +
         theme_bw() +
         theme(#legend.position = "none",
               #legend.title = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()))
ggsave(paste0(resdir, "Fitness-vs-protSpecificity_gd_24h.pdf"), width = 5, height = 3)


(p2 <- comb.gd %>%
                mutate(FDR0.05_count_bin = str_c("n", FDR0.05_count)) %>%
                mutate(FDR0.05_count_bin = ifelse(FDR0.05_count_bin %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
                mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
                       essential = forcats::fct_rev(as.factor(essential))) %>%
                ggplot(aes(x = FDR0.05_count_bin, fill = essential)) +
                geom_bar(position = "fill") +
                scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
                scale_y_continuous(breaks = c(0, 0.5, 1)) +
                labs(x = "Number of proteins affected by gRNA", 
                     y = "Fraction of gRNAs",
                     fill = "gRNA target gene") +
                theme_bw() +
                theme(#legend.title = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()))
ggsave(paste0(resdir, "Fraction-vs-protSpecificity_gd_24h.pdf"), width = 5, height = 2)


# Combined plot
ggarrange(p1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), NULL, p2, 
          ncol = 1, heights = c(2, -0.3, 1.5), widths = c(1, 1, 1), align = "hv",
          common.legend = T, legend = "right", legend.grob = get_legend(p2))
ggsave(paste0(resdir, "Fitness-vs-protSpecificity_gd_24h_WithFraction.pdf"), width = 5, height = 4.5)




# Counts -----------------------------------------------------------------------

counts <- comb.gd %>%
        mutate(bin = case_when(FDR0.05_count <3 ~ "spec",
                               FDR0.05_count >2 ~ "broad",
                               TRUE ~ NA_character_)) %>%
        mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
               essential = forcats::fct_rev(as.factor(essential))) %>%
        filter(!is.na(bin)) %>%
        group_by(bin, essential) %>%
        summarize(n = n()) %>% 
        print()

# bin   essential        n
# broad Nonessential    67
# broad Essential      194
# spec  Nonessential 10840
# spec  Essential     5249

# Fraction of essential genes among genes that upon perturbation affect 3 or more proteins:
194/(194+67) # 74%

# Fraction of essential genes among genes that upon perturbation affect less than 3 proteins:
5249/(5249+10840) # 33%

fisher.test(matrix(counts$n, ncol = 2))
# p-value < 2.2e-16
# odds ratio = 0.1672554 

chisq.test(matrix(counts$n, ncol = 2))
# X-squared = 199.27, df = 1, p-value < 2.2e-16




# Reads per gRNA at "pre" time point -------------------------------------------
# Addressing concern by reviewer that the initial frequencies of the n6/n7/n8+ 
# bins are already skewed/lower and therefore drop out more due to stochastic
# effects.

# Import fitness data with more info
fitn.gd2 <- readRDS("results/fitness_all4/all4_Ratios_gd_allInfo.RDS") %>% print()


# Combine protein and fitness data
comb.gd2 <- prot.gd %>%
        left_join(fitn.gd2, by = "guide") %>% 
        print()

comb.gd2 %>%
                mutate(FDR0.05_count_bin = str_c("n", FDR0.05_count)) %>%
                mutate(FDR0.05_count_bin = ifelse(FDR0.05_count_bin %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count_bin)) %>%
                mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
                       essential = forcats::fct_rev(as.factor(essential))) %>%
                ggplot() +
                geom_point(aes(x = FDR0.05_count_bin, y = Tpre_A, color = essential),
                           position = "jitter", alpha = 0.7) +
                geom_point(aes(x = FDR0.05_count_bin, y = Tpre_B, color = essential),
                           position = "jitter", alpha = 0.7) +
                geom_point(aes(x = FDR0.05_count_bin, y = Tpre_C, color = essential),
                           position = "jitter", alpha = 0.7) +
                geom_point(aes(x = FDR0.05_count_bin, y = Tpre_D, color = essential),
                           position = "jitter", alpha = 0.7) +
                geom_boxplot(aes(x = FDR0.05_count_bin, y = Tpre_A),
                             alpha = 0.8, outlier.shape = NA, show.legend = F) +
                coord_cartesian(ylim = c(0, 2000)) +
                labs(x = "Number of proteins affected by gRNA",
                     y = "Number of reads per gRNA",
                     color = "gRNA target gene") +
                theme_bw() +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
ggsave(paste0(resdir, "Reads-vs-protSpecificity_gd_preGalInd_boxplotA.pdf"), width = 5, height = 3)




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-05_SessionInfo.txt")

