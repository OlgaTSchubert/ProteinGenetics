# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE21_FitnessScreen")
getwd()


# General results directory
resdir0 <- "results/"
dir.create(resdir0)


# Results directory
resdir <- "results/fitness_wES/"
dir.create(resdir)


# Read count table
res <- readRDS("data/readCounts_wES.RDS") %>% 
        select(-starts_with("T48h")) %>% print()




# Transformation ---------------------------------------------------------------

# Simplify names
res2 <- res %>%
        dplyr::rename(Tpre_A = Tpre_wES_A,
                      Tpre_B = Tpre_wES_B,
                      T00h_A = T00h_wES_A,
                      T00h_B = T00h_wES_B,
                      T08h_A = T08h_wES_A,
                      T08h_B = T08h_wES_B,
                      T24h_A = T24h_wES_A,
                      T24h_B = T24h_wES_B) %>% print()


# To normalize across samples get read counts for controls only
ctrsums <- res2 %>%
        filter(set %in% c("Ctr1", "Ctr2", "Ctr3")) %>%
        summarise_at(vars(T00h_A:Tpre_B), sum, na.rm = TRUE) %>% print()


# Normalize, log2fc, mean
res3 <- res2 %>%
        mutate(Tpre_A_norm = Tpre_A/ctrsums$Tpre_A*1E6,
               Tpre_B_norm = Tpre_B/ctrsums$Tpre_B*1E6,
               T00h_A_norm = T00h_A/ctrsums$T00h_A*1E6,
               T00h_B_norm = T00h_B/ctrsums$T00h_B*1E6,
               T08h_A_norm = T08h_A/ctrsums$T08h_A*1E6,
               T08h_B_norm = T08h_B/ctrsums$T08h_B*1E6,
               T24h_A_norm = T24h_A/ctrsums$T24h_A*1E6,
               T24h_B_norm = T24h_B/ctrsums$T24h_B*1E6) %>%
        mutate(Tpre_A_log2fc = log2(Tpre_A_norm/Tpre_A_norm),
               Tpre_B_log2fc = log2(Tpre_B_norm/Tpre_B_norm),
               T00h_A_log2fc = log2(T00h_A_norm/Tpre_A_norm),
               T00h_B_log2fc = log2(T00h_B_norm/Tpre_B_norm),
               T08h_A_log2fc = log2(T08h_A_norm/Tpre_A_norm),
               T08h_B_log2fc = log2(T08h_B_norm/Tpre_B_norm),
               T24h_A_log2fc = log2(T24h_A_norm/Tpre_A_norm),
               T24h_B_log2fc = log2(T24h_B_norm/Tpre_B_norm)) %>%
        mutate(T00h_log2fc = rowMeans(cbind(T00h_A_log2fc, T00h_B_log2fc), na.rm = T),
               T08h_log2fc = rowMeans(cbind(T08h_A_log2fc, T08h_B_log2fc), na.rm = T),
               T24h_log2fc = rowMeans(cbind(T24h_A_log2fc, T24h_B_log2fc), na.rm = T)) %>%
        select(c(guide, set, geneSys, gene, T00h_log2fc:T24h_log2fc)) %>% print()
write_excel_csv(res3, paste0(resdir, "wES_Ratios_gd_wide.csv"))


# Convert to long format
res4 <- res3 %>%
        pivot_longer(-c(guide, set, geneSys, gene),
                     names_to = c("timepoint", ".value"),
                     names_pattern = "(.+)_(.+)") %>% 
        mutate(dropout = ifelse(log2fc < -0.5, T, F)) %>% print()
saveRDS(res4, paste0(resdir, "wES_Ratios_gd_long.RDS"))


# Summarization per gene
res4_gn <- res4 %>%
        filter(set %in% c("eStops", "eProvs", "neStops", "neProvs")) %>%
        mutate(essential = ifelse(set %in% c("eStops", "eProvs"), T, F)) %>%
        select(-c(guide, set, dropout)) %>%
        group_by(geneSys, timepoint) %>%
        filter(abs(log2fc) == max(abs(log2fc))) %>%
        distinct() %>%  # because filter on max leads to duplicates when tie
        arrange(geneSys, -log2fc) %>%
        mutate(dropout = ifelse(log2fc < -0.5, T, F)) %>% print()
saveRDS(res4_gn, paste0(resdir, "wES_Ratios_gn_long.RDS"))

res4_gn %>%
        select(-dropout) %>%
        pivot_wider(names_from = timepoint, values_from = log2fc,
                    names_glue = "{timepoint}_log2fc", names_sort = T) %>%
        write_excel_csv(paste0(resdir, "wES_Ratios_gn_wide.csv"))




# Fitness jitter plots ---------------------------------------------------------

# Change names for sets and timepoints
res5 <- res4 %>% 
        mutate(set = case_when(set == "Ctr1" ~ "Non-yeast",
                               set == "Ctr2" ~ "No C",
                               set == "Ctr3" ~ "Synonymous",
                               set == "eStops"  ~ "Ess. stop",
                               set == "eProvs"  ~ "Ess. provean",
                               set == "neStops" ~ "NonEss. stop",
                               set == "neProvs" ~ "Noness. provean",
                               TRUE ~ NA_character_),
               set = factor(set, levels = c("Non-yeast", "No C", "Synonymous", 
                                            "Ess. stop", "Ess. provean", 
                                            "NonEss. stop", "Noness. provean")),
               timepoint = case_when(timepoint == "T00h" ~ "24 hours",
                                     timepoint == "T08h" ~ "32 hours",
                                     timepoint == "T24h" ~ "48 hours",
                                     TRUE ~ NA_character_)) %>% print()

res5 %>%
        filter(timepoint != "32 hours") %>%
        filter(set %in% c("Non-yeast", "No C", "Synonymous", "Ess. stop")) %>%
        ggplot(aes(x = set, y = log2fc)) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        geom_jitter(alpha = 0.1) +
        geom_violin(trim = T, alpha = 0.5) +
        labs(x = "", y = "Log2FC") +
        facet_wrap(~ timepoint) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(resdir, "wES_FitnessJitter_ESonly.pdf"), height = 4, width = 4.5)
ggsave(paste0(resdir, "wES_FitnessJitter.pdf"), height = 4, width = 6)




# Dropout fraction plot --------------------------------------------------------
# Dropout if log2fc < -0.5

res5 %>%
        filter(set %in% c("Non-yeast", "No C", "Synonymous", "Ess. stop")) %>%
        group_by(timepoint, set, dropout) %>%
        summarize(n = n()) %>%
        ggplot(aes(x = set, y = n, fill = dropout)) +
        geom_col(position = "fill") +
        labs(x = "", y = "Fraction of gRNAs") +
        facet_wrap(~ timepoint) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(resdir, "wES_DropoutFraction_ESonly.pdf"), height = 3, width = 4)
ggsave(paste0(resdir, "wES_DropoutFraction.pdf"), height = 3, width = 6)


# Numbers
res5.sum <- res5 %>%
        filter(set %in% c("Non-yeast", "No C", "Synonymous", "Ess. stop")) %>% 
        mutate(set2 = ifelse(set == "Ess. stop", "Ess.stop", "Control")) %>%
        group_by(timepoint, set2, dropout) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        filter(timepoint == "48 hours") %>%
        print()

ES.dropout <- res5.sum %>% filter(set2 == "Ess.stop", dropout == T) %>% pull(n)
ES.notdo <- res5.sum %>% filter(set2 == "Ess.stop", dropout == F) %>% pull(n)
ES.dropout/(ES.dropout+ES.notdo) # 0.59

res5.sum.red <- res5.sum %>%
        select(-"timepoint") %>%
        pivot_wider(names_from = set2, values_from = n) %>%
        print()
        
chisq.test(res5.sum.red[, 2:3])
# X-squared = 1061, df = 1, p-value < 2.2e-16


        

# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-03a_SessionInfo.txt")

