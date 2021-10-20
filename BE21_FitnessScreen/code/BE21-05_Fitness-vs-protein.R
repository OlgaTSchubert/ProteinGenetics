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




# Plotting fitness as a function of specificity --------------------------------

fitnessVsSpecificityPlot.gd <- function(input, log2fc, ess = F) {
        
        input <- input %>%
                mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
                mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count))

        if(ess == F) {
                input %>%
                        ggplot(aes(x = FDR0.05_count, y = {{log2fc}})) +
                        geom_point(position = "jitter", alpha = 0.3) +
                        geom_boxplot(alpha = 0.8, outlier.shape = NA) +
                        coord_cartesian(ylim = c(-2.5, 1)) +
                        labs(x = "Number of proteins affected",
                             y = "Log2FC") +
                        theme_bw() +
                        theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank())
                
        } else {
                input %>%
                        mutate(essential = ifelse(essential == TRUE, "essential", "nonessential"),
                               essential = forcats::fct_rev(as.factor(essential))) %>%
                        ggplot(aes(x = FDR0.05_count, y = {{log2fc}}, color = essential)) +
                        geom_point(position = position_jitterdodge()) +
                        geom_boxplot(alpha = 0.8, outlier.shape = NA, position = "dodge") +
                        scale_color_manual(values = c("darkseagreen3", "darkseagreen4")) +
                        coord_cartesian(ylim = c(-2.5, 1)) +
                        labs(x = "Number of proteins affected",
                             y = "Log2FC") +
                        theme_bw() +
                        theme(legend.title = element_blank(),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank())
        }
}

(p1 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T00h_log2fc, ess = F))
(p2 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T08h_log2fc, ess = F))
(p3 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T24h_log2fc, ess = F))

(p123 <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1))
ggsave(plot = p123, paste0(resdir, "Fitness-vs-protSpecificity_gd.pdf"), width = 9, height = 3)
ggsave(plot = p3, paste0(resdir, "Fitness-vs-protSpecificity_gd_24h.pdf"), width = 3, height = 3)


(p4 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T00h_log2fc, ess = T))
(p5 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T08h_log2fc, ess = T))
(p6 <- fitnessVsSpecificityPlot.gd(input = comb.gd, log2fc = T24h_log2fc, ess = T))

(p456 <- ggarrange(p4, p5, p6, ncol = 3, nrow = 1))
ggsave(plot = p456, paste0(resdir, "Fitness-vs-protSpecificity_gd_ess.pdf"), width = 14, height = 3)
ggsave(plot = p6, paste0(resdir, "Fitness-vs-protSpecificity_gd_ess_24h.pdf"), width = 4.5, height = 3)






# Combined plot ----------------------------------------------------------------

(p6a <- comb.gd %>%
        mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
        mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
        mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
               essential = forcats::fct_rev(as.factor(essential))) %>%
        ggplot(aes(x = FDR0.05_count, y = T24h_log2fc, color = essential)) +
        geom_point(position = position_jitterdodge()) +
        geom_boxplot(alpha = 0.8, outlier.shape = NA, 
                     position = "dodge", show.legend = F) +
        scale_color_manual(values = c("darkseagreen3", "darkseagreen4")) +
        coord_cartesian(ylim = c(-2.5, 2)) +
        labs(x = "Number of proteins affected",
             y = "Log2FC") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()))

(p7 <- comb.gd %>%
                mutate(FDR0.05_count = str_c("n", FDR0.05_count)) %>%
                mutate(FDR0.05_count = ifelse(FDR0.05_count %in% c("n8", "n9", "n10"), "n8+", FDR0.05_count)) %>%
                mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
                       essential = forcats::fct_rev(as.factor(essential))) %>%
                ggplot(aes(x = FDR0.05_count, fill = essential)) +
                geom_bar(position = "fill") +
                scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
                scale_y_continuous(breaks = c(0, 0.5, 1)) +
                labs(x = "Number of proteins affected", y = "Fraction") +
                theme_bw() +
                theme(legend.title = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()))

ggarrange(p6 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
          NULL, 
          p7, 
          ncol = 1, heights = c(2, -0.3, 1), widths = c(1, 1, 1), align = "hv",
          common.legend = T, legend = "right", legend.grob = get_legend(p7))
ggsave(paste0(resdir, "Fitness-vs-protSpecificity_gd_ess_24h_WithFraction.pdf"), width = 4.5, height = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-05_SessionInfo.txt")

