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
resdir <- "results/fitness_all4/"
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
res3a <- res2 %>%
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
        mutate(T00h_log2fc = rowMeans(cbind(T00h_A_log2fc, T00h_B_log2fc, T00h_C_log2fc, T00h_D_log2fc), na.rm = T),
               T08h_log2fc = rowMeans(cbind(T08h_A_log2fc, T08h_B_log2fc, T08h_C_log2fc, T08h_D_log2fc), na.rm = T),
               T24h_log2fc = rowMeans(cbind(T24h_A_log2fc, T24h_B_log2fc, T24h_C_log2fc, T24h_D_log2fc), na.rm = T)) %>%
        print()
saveRDS(res3a, paste0(resdir, "all4_Ratios_gd_allInfo.RDS"))


res3 <- res3a %>%
        select(c(guide, set, geneSys, gene, T00h_log2fc:T24h_log2fc)) %>%
        print()
write_excel_csv(res3, paste0(resdir, "all4_Ratios_gd_wide.csv"))


# Convert to long format
res4 <- res3 %>%
        pivot_longer(-c(guide, set, geneSys, gene),
                     names_to = c("timepoint", ".value"),
                     names_pattern = "(.+)_(.+)") %>% 
        mutate(dropout = ifelse(log2fc < -0.5, T, F)) %>% print()
saveRDS(res4, paste0(resdir, "all4_Ratios_gd_long.RDS"))


# Summarization per gene
res4_gn <- res4 %>%
        mutate(essential = ifelse(set %in% c("eStops", "eProvs"), T, F)) %>%
        select(-c(guide, set, dropout)) %>%
        group_by(gene, timepoint) %>%
        filter(abs(log2fc) == max(abs(log2fc))) %>%
        distinct() %>%  # because filter on max leads to duplicates when tie
        arrange(gene, -log2fc) %>%
        mutate(dropout = ifelse(log2fc < -0.5, T, F)) %>% print()
saveRDS(res4_gn, paste0(resdir, "all4_Ratios_gn_long.RDS"))

res4_gn %>%
        select(-dropout) %>%
        pivot_wider(names_from = timepoint, values_from = log2fc,
                    names_glue = "{timepoint}_log2fc", names_sort = T) %>%
        write_excel_csv(paste0(resdir, "all4_Ratios_gn_wide.csv"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-03b_SessionInfo.txt")

