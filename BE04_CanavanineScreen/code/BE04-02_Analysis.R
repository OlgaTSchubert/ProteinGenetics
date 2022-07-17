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

table(res4$consequence, res4$enriched)
#               FALSE TRUE
# nonsynonymous    54   11
# synonymous       19    0
# stop              2    4


# Stats by guide count
res4.add <- res4 %>%
        mutate(consequence_stop = ifelse(consequence == "stop", T, F),
               consequence_syn = ifelse(consequence == "synonymous", T, F)) %>%
        print()

table(res4.add$consequence_stop, res4.add$enriched)
#       FALSE TRUE
# FALSE    73   11
# TRUE      2    4

fisher.test(res4.add$consequence_stop, res4.add$enriched)
# p-value = 0.006454
# odds ratio = 12.6

table(res4.add$consequence_syn, res4.add$enriched)
#       FALSE TRUE
# FALSE    56   15
# TRUE     19    0

fisher.test(res4.add$consequence_syn, res4.add$enriched)
# p-value = 0.03419
# odds ratio = 0


# Stats by read count
res4.sum <- res4 %>%
        pivot_longer(-c(guide, geneSys, gene, mut1, mut2, consequence, provean, maxAbsProvean, enriched), 
                     names_to = c("time", "condition"),
                     names_sep = "_",
                     values_to = "reads") %>%
        group_by(consequence, condition) %>%
        summarize(readsum = sum(reads)) %>%
        pivot_wider(names_from = consequence, values_from = readsum) %>%
        mutate(notstop = nonsynonymous+synonymous,
               notsyn  = nonsynonymous+stop) %>%
        print()
# condition  nonsynonymous synonymous  stop notstop notsyn
# Canavanine         32397.      4587. 7007.  36984. 39404.
# Control            23652.      6803. 1979.  30455. 25631.

chisq.test(res4.sum[, c("stop", "notstop")])
# X-squared = 1736.8, df = 1, p-value < 2.2e-16

chisq.test(res4.sum[, c("synonymous", "notsyn")])
# X-squared = 1636.8, df = 1, p-value < 2.2e-16


# Are guides introducing mutations at conserved sites (maxAbsProvean > 5) more likely
# to lead to canavanine resistance?
res4.ns <- res4 %>%
        filter(consequence == "nonsynonymous") %>%
        mutate(provean5 = ifelse(maxAbsProvean > 5, "hiProv", "loProv")) %>%
        print()

table(res4.ns$provean5, res4.ns$enriched)
#        FALSE TRUE
# hiProv    11    7
# loProv    43    4

fisher.test(res4.ns$provean5, res4.ns$enriched)
# p-value = 0.007196
# odds ratio = 0.1518771 




# Plots ------------------------------------------------------------------------

# Jitter plot
res4 %>%
        mutate(`Day 0_Control_log2fc`    = log2(`Day 0_Control`   /`Day 0_Control`),
               `Day 1_Canavanine_log2fc` = log2(`Day 1_Canavanine`/`Day 0_Control`),
               `Day 1_Control_log2fc`    = log2(`Day 1_Control`   /`Day 0_Control`),
               `Day 2_Canavanine_log2fc` = log2(`Day 2_Canavanine`/`Day 0_Control`),
               `Day 2_Control_log2fc`    = log2(`Day 2_Control`   /`Day 0_Control`)) %>%
        select(-(`Day 0_Control`:`Day 0_Canavanine`)) %>%
        pivot_longer(-c(guide, geneSys, gene, mut1, mut2, consequence, provean, 
                        maxAbsProvean, enriched), 
                     names_to = c("time", "condition", "remove"),
                     names_sep = "_",
                     values_to = "Log2FC") %>%
        select(-remove) %>%
        filter(!is.na(enriched)) %>% 
        filter(condition != "Control") %>%
        mutate(consequence = case_when(consequence == "synonymous" ~ "Synonymous",
                                       consequence == "nonsynonymous" ~ "Missense",
                                       consequence == "stop" ~ "Stop"),
               consequence = factor(consequence, levels = c("Synonymous", "Missense", "Stop"))) %>%
        mutate(time = case_when(time == "Day 0" ~ "0 hours",
                                time == "Day 1" ~ "24 hours",
                                time == "Day 2" ~ "48 hours")) %>%
        ggplot(aes(x = consequence, y = Log2FC)) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        geom_jitter(alpha = 0.3) +
        geom_violin(trim = F, alpha = 0.2) +
        labs(x = "", y = "Log2FC") +
        facet_wrap(~ time) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(resdir, "Jitter.pdf"), height = 3, width = 3)


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

