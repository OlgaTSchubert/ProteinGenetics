# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE04_CanavanineScreen")
getwd()


# Results directory
resdir <- "results/"


# Processed reads files
samplesDF <- readRDS("results/ReadsRDS/Canavanine_Reads.RDS") %>% print()
reads     <- samplesDF$results %>% print()


# Guides table
# Add column with shorter guide sequence to account for sequencing issue
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set == "Can1Ade2") %>%
        select(guide, geneSys, gene, mut1, mut2, consequence, provean, maxAbsProvean) %>%
        mutate(guide18 = str_sub(guide, start = 1L, end = 18L)) %>% print()
        



# Read counts per guide --------------------------------------------------------

# Initiate results table (based on guides table)
res <- guides


# Add read counts to results table
for (i in seq_along(reads)) {
        
        readsX <- reads[[i]] %>%
                as_tibble() %>%
                dplyr::rename(guide18 = guide) %>%
                group_by(guide18) %>%
                summarize(n_reads = n()) 

        res <- left_join(res, readsX, by = "guide18") %>% 
                rename_at(vars(starts_with("n_reads")), ~names(reads)[i])
}


# Replace NAs with 1 in read counts columns
res <- res %>% mutate_if(is.integer, ~replace_na(., 1)) %>% print()




# Transformation ---------------------------------------------------------------

# Normalize by total read counts per sample
totsums <- summarise_at(res, vars(T00h_Ctr:T48h_Ctr), sum, na.rm = TRUE) %>% print()

res2 <- res %>%
        mutate(T00h_Ctr_norm = T00h_Ctr/totsums$T00h_Ctr*20000,
               T24h_Can_norm = T24h_Can/totsums$T24h_Can*20000,
               T24h_Ctr_norm = T24h_Ctr/totsums$T24h_Ctr*20000,
               T48h_Can_norm = T48h_Can/totsums$T48h_Can*20000,
               T48h_Ctr_norm = T48h_Ctr/totsums$T48h_Ctr*20000) %>%
        print()


# Clean up table
res3 <- res2 %>% 
        filter(geneSys == "YEL063C") %>% 
        select(-c(guide18, T00h_Ctr:T48h_Ctr)) %>%
        dplyr::rename(`Day 0_Control`    = T00h_Ctr_norm,
                      `Day 1_Control`    = T24h_Ctr_norm,
                      `Day 2_Control`    = T48h_Ctr_norm,
                      `Day 1_Canavanine` = T24h_Can_norm,
                      `Day 2_Canavanine` = T48h_Can_norm) %>%
        mutate(`Day 0_Canavanine` = `Day 0_Control`) %>%
        print()


# Define a guide as "enriched" if 2x more reads at day 2 compared to day 0
res4 <- res3 %>%  
        mutate(enriched = ifelse(`Day 2_Canavanine` / `Day 0_Canavanine` > 2, T, F)) %>%
        print()
write_excel_csv(res4, paste0(resdir, "results.csv"))




# Analyses ---------------------------------------------------------------------

# Are guides introducing stop or nonsynonymous mutations more likely to lead to 
# canavanine resistance than guides introducing synonymous mutations?

table(res4$enriched, res4$consequence)
#       nonsynonymous synonymous stop
# FALSE            54         19    2
# TRUE             11          0    4

4/(2+4)     # 0.67 ( 4 out of  6 stop guides are enriched in Can condition)
11/(54+11)  # 0.17 (11 out of 65 nonsynonymous guides are enriched in Can condition)
0/(19+0)    # 0    ( 0 out of 19 synonymous guides are enriched in Can condition)


# Are guides introducing mutations at conserved sites (Provean > 5) more likely
# to lead to canavanine resistance?
res4.ns <- res4 %>%
        filter(consequence == "nonsynonymous") %>%
        mutate(provean5 = ifelse(maxAbsProvean > 5, "hiProv", "loProv")) %>%
        print()

table(res4.ns$enriched, res4.ns$provean5)
#       hiProv loProv
# FALSE     11     43
# TRUE       7      4

7/(11+7)  # 0.39 (7 out of 18 hiProv guides are enriched in Can condition)
4/(43+4)  # 0.09 (4 out of 47 loProv guides are enriched in Can condition)


# Among the 11 nonsynonymous guides that are enriched, how many are highly 
# conserved (Provean > 5)?
sum(res4.ns$enriched == T & res4.ns$maxAbsProvean > 5)/sum(res4.ns$enriched == T)  # 0.64
sum(res4.ns$enriched == F & res4.ns$maxAbsProvean > 5)/sum(res4.ns$enriched == F)  # 0.20





# Plots ------------------------------------------------------------------------

# Stacked columns
res4 %>%
        pivot_longer(-c(guide, geneSys, gene, mut1, mut2, consequence, provean, maxAbsProvean, enriched), 
                     names_to = c("time", "condition"),
                     names_sep = "_",
                     values_to = "reads") %>%
        mutate(condition = factor(condition, levels = c("Control", "Canavanine")),
               consequence = case_when(consequence == "synonymous" ~ "Synonymous",
                                       consequence == "nonsynonymous" ~ "Missense",
                                       consequence == "stop" ~ "Stop")) %>%
        mutate(consequence = factor(consequence, levels = c("Synonymous", "Missense", "Stop"))) %>%
        mutate(time = case_when(time == "Day 0" ~ "0 hours",
                                time == "Day 1" ~ "24 hours",
                                time == "Day 2" ~ "48 hours")) %>%
        ggplot(aes(x = time, y = reads, fill = consequence)) +
        geom_col(position = "fill") +
        scale_fill_manual(values = c("#E0F3DB", "#A8DDB5", "#43A2CA")) +
        labs(x = "", y = "Fraction of reads", fill = "Mutation type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) +
        facet_wrap(~ condition, )
ggsave(paste0(resdir, "Columns3.pdf"), width = 4, height = 3)


# Stacked columns (missense guides split up by enriched and nonenriched)
res4 %>%
        pivot_longer(-c(guide, geneSys, gene, mut1, mut2, consequence, provean, maxAbsProvean, enriched), 
                     names_to = c("time", "condition"),
                     names_sep = "_",
                     values_to = "reads") %>%
        filter(!is.na(enriched)) %>%
        mutate(condition = factor(condition, levels = c("Control", "Canavanine")),
               consequence = case_when(consequence == "synonymous" ~ "Synonymous",
                                       consequence == "nonsynonymous" & enriched == T ~ "Missense*",
                                       consequence == "nonsynonymous" & enriched == F ~ "Missense",
                                       consequence == "stop" ~ "Stop"),
               consequence = factor(consequence, levels = c("Synonymous", "Missense", "Missense*", "Stop"))) %>%
        mutate(time = case_when(time == "Day 0" ~ "0 hours",
                                time == "Day 1" ~ "24 hours",
                                time == "Day 2" ~ "48 hours")) %>%
        ggplot(aes(x = time, y = reads, fill = consequence)) +
        geom_col(position = "fill") +
        scale_fill_manual(values = c("#F0F9E8", "#BAE4BC", "#7BCCC4", "2B8CBE")) +
        labs(x = "", y = "Fraction of reads", fill = "Mutation type") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) +
        facet_wrap(~ condition, )
ggsave(paste0(resdir, "Columns4.pdf"), width = 4, height = 3)




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE04-02_SessionInfo.txt")

