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
resdir <- "results/fitness_replComparison/"
dir.create(resdir)


# Read count table
res <- readRDS("data/readCounts.RDS") %>% 
        select(-starts_with("T48h")) %>% print()




# Transformation ---------------------------------------------------------------

# Rename the 2 replicates each of wES and woES to 4 replicates
res2 <- res %>%
        filter(set %in% c("eProvs", "neStops", "neProvs")) %>%
        dplyr::rename(Tpre_A = Tpre_wES_A,
                      Tpre_B = Tpre_wES_B,
                      Tpre_C = Tpre_woES_A,
                      Tpre_D = Tpre_woES_B,
                      T00h_A = T00h_wES_A,
                      T00h_B = T00h_wES_B,
                      T00h_C = T00h_woES_A,
                      T00h_D = T00h_woES_B,
                      T08h_A = T08h_wES_A,
                      T08h_B = T08h_wES_B,
                      T08h_C = T08h_woES_A,
                      T08h_D = T08h_woES_B,
                      T24h_A = T24h_wES_A,
                      T24h_B = T24h_wES_B,
                      T24h_C = T24h_woES_A,
                      T24h_D = T24h_woES_B) %>%
        print()


# To normalize across samples get total read counts
totsums <- res2 %>%
        summarise_at(vars(T00h_A:Tpre_D), sum, na.rm = TRUE) %>% print()


# Normalize, log2fc, mean
res3 <- res2 %>%
        mutate(Tpre_A_norm = Tpre_A/totsums$Tpre_A*1E6,
               Tpre_B_norm = Tpre_B/totsums$Tpre_B*1E6,
               Tpre_C_norm = Tpre_C/totsums$Tpre_C*1E6,
               Tpre_D_norm = Tpre_D/totsums$Tpre_D*1E6,
               T00h_A_norm = T00h_A/totsums$T00h_A*1E6,
               T00h_B_norm = T00h_B/totsums$T00h_B*1E6,
               T00h_C_norm = T00h_C/totsums$T00h_C*1E6,
               T00h_D_norm = T00h_D/totsums$T00h_D*1E6,
               T08h_A_norm = T08h_A/totsums$T08h_A*1E6,
               T08h_B_norm = T08h_B/totsums$T08h_B*1E6,
               T08h_C_norm = T08h_C/totsums$T08h_C*1E6,
               T08h_D_norm = T08h_D/totsums$T08h_D*1E6,
               T24h_A_norm = T24h_A/totsums$T24h_A*1E6,
               T24h_B_norm = T24h_B/totsums$T24h_B*1E6,
               T24h_C_norm = T24h_C/totsums$T24h_C*1E6,
               T24h_D_norm = T24h_D/totsums$T24h_D*1E6) %>%
        mutate(Tpre_A_log2fc = log2(Tpre_A_norm/Tpre_A_norm),
               Tpre_B_log2fc = log2(Tpre_B_norm/Tpre_B_norm),
               Tpre_C_log2fc = log2(Tpre_C_norm/Tpre_C_norm),
               Tpre_D_log2fc = log2(Tpre_D_norm/Tpre_D_norm),
               T00h_A_log2fc = log2(T00h_A_norm/Tpre_A_norm),
               T00h_B_log2fc = log2(T00h_B_norm/Tpre_B_norm),
               T00h_C_log2fc = log2(T00h_C_norm/Tpre_C_norm),
               T00h_D_log2fc = log2(T00h_D_norm/Tpre_D_norm),
               T08h_A_log2fc = log2(T08h_A_norm/Tpre_A_norm),
               T08h_B_log2fc = log2(T08h_B_norm/Tpre_B_norm),
               T08h_C_log2fc = log2(T08h_C_norm/Tpre_C_norm),
               T08h_D_log2fc = log2(T08h_D_norm/Tpre_D_norm),
               T24h_A_log2fc = log2(T24h_A_norm/Tpre_A_norm),
               T24h_B_log2fc = log2(T24h_B_norm/Tpre_B_norm),
               T24h_C_log2fc = log2(T24h_C_norm/Tpre_C_norm),
               T24h_D_log2fc = log2(T24h_D_norm/Tpre_D_norm)) %>%
        print()

res4_ABCD <- res3 %>%
        mutate(T00h_log2fc = rowMeans(cbind(T00h_A_log2fc, T00h_B_log2fc, T00h_C_log2fc, T00h_D_log2fc), na.rm = T),
               T08h_log2fc = rowMeans(cbind(T08h_A_log2fc, T08h_B_log2fc, T08h_C_log2fc, T08h_D_log2fc), na.rm = T),
               T24h_log2fc = rowMeans(cbind(T24h_A_log2fc, T24h_B_log2fc, T24h_C_log2fc, T24h_D_log2fc), na.rm = T)) %>%
        select(c(guide, set, geneSys, gene, T00h_log2fc:T24h_log2fc)) %>%
        print()

res4_AB_CD <- res3 %>%
        mutate(T00h_log2fc_AB = rowMeans(cbind(T00h_A_log2fc, T00h_B_log2fc), na.rm = T),
               T08h_log2fc_AB = rowMeans(cbind(T08h_A_log2fc, T08h_B_log2fc), na.rm = T),
               T24h_log2fc_AB = rowMeans(cbind(T24h_A_log2fc, T24h_B_log2fc), na.rm = T),
               T00h_log2fc_CD = rowMeans(cbind(T00h_C_log2fc, T00h_D_log2fc), na.rm = T),
               T08h_log2fc_CD = rowMeans(cbind(T08h_C_log2fc, T08h_D_log2fc), na.rm = T),
               T24h_log2fc_CD = rowMeans(cbind(T24h_C_log2fc, T24h_D_log2fc), na.rm = T)) %>%
        select(c(guide, set, geneSys, gene, T00h_log2fc_AB:T24h_log2fc_CD)) %>%
        print()


# Convert to long format and add "dropout" column
res5_ABCD <- res4_ABCD %>% 
        pivot_longer(-c(guide, set, geneSys, gene),
                     names_to = c("timepoint", ".value"),
                     names_pattern = "(.+h)_(log2fc)") %>%
        mutate(dropout = ifelse(log2fc < -0.5, T, F)) %>% 
        print()


res5_AB_CD <- res4_AB_CD %>%
        pivot_longer(-c(guide, set, geneSys, gene),
                     names_to = c("timepoint", ".value"),
                     names_pattern = "(.+h)_(log2fc_..)") %>% 
        mutate(dropout_AB = ifelse(log2fc_AB < -0.5, T, F),
               dropout_CD = ifelse(log2fc_CD < -0.5, T, F)) %>% 
        print()


# Combine tables and plot
res5_AB_CD %>%
        inner_join(res5_ABCD, by = c("guide", "set", "geneSys", "gene", "timepoint")) %>%
        #filter(dropout == T) %>%
        ggplot(aes(x = log2fc_AB, y = log2fc_CD)) +
        geom_hline(yintercept = 0, alpha = 0.3) +
        geom_vline(xintercept = 0, alpha = 0.3) +
        #geom_smooth(method = "lm", color = "black") +
        ggpubr::stat_cor(color = "grey25",
                         r.digits = 3, 
                         label.x.npc = "left", 
                         label.y.npc = "top", 
                         label.sep = "\n") +
        geom_point(alpha = 0.2) +
        labs(x = "Log2FC for replicates A, B",
             y = "Log2FC for replicates C, D") +
        lims(x = c(-5, 3), y = c(-5, 3)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "AB-vs-CD_gd.pdf"), width = 4, height = 4)
ggsave(paste0(resdir, "AB-vs-CD_gd.jpg"), width = 4, height = 4)


# Summarization per gene
res5_AB_gn <- res5_AB_CD %>%
        select(-log2fc_CD, -dropout_AB, -dropout_CD) %>%
        rename(log2fc = log2fc_AB) %>%
        mutate(essential = ifelse(set %in% c("eStops", "eProvs"), T, F)) %>%
        select(-c(guide, set)) %>%
        group_by(gene, timepoint) %>%
        filter(abs(log2fc) == max(abs(log2fc))) %>%
        distinct() %>%  # because filter on max leads to duplicates when tie
        arrange(gene, -log2fc) %>% print()

res5_CD_gn <- res5_AB_CD %>%
        select(-log2fc_AB, -dropout_AB, -dropout_CD) %>%
        rename(log2fc = log2fc_CD) %>%
        mutate(essential = ifelse(set %in% c("eStops", "eProvs"), T, F)) %>%
        select(-c(guide, set)) %>%
        group_by(gene, timepoint) %>%
        filter(abs(log2fc) == max(abs(log2fc))) %>%
        distinct() %>%  # because filter on max leads to duplicates when tie
        arrange(gene, -log2fc) %>% print()


# Plot correlation per gene
res5_AB_gn %>%
        full_join(res5_CD_gn, by = c("geneSys", "gene", "timepoint", "essential"),
                  suffix = c("_AB", "_CD")) %>%
        ggplot(aes(x = log2fc_AB, y = log2fc_CD)) +
        geom_hline(yintercept = 0, alpha = 0.3) +
        geom_vline(xintercept = 0, alpha = 0.3) +
        #geom_smooth(method = "lm", color = "black") +
        ggpubr::stat_cor(color = "grey25",
                         r.digits = 3, 
                         label.x.npc = "left", 
                         label.y.npc = "top", 
                         label.sep = "\n") +
        geom_point(alpha = 0.2) +
        labs(x = "Log2FC for replicates A, B",
             y = "Log2FC for replicates C, D") +
        lims(x = c(-5, 3), y = c(-5, 3)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "AB-vs-CD_gn.pdf"), width = 4, height = 4)
ggsave(paste0(resdir, "AB-vs-CD_gn.jpg"), width = 4, height = 4)





# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-03c_SessionInfo.txt")

