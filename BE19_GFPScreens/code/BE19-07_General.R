# Libraries --------------------------------------------------------------------

library(tidyverse)
library(egg)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/general/"
dir.create(resdir)


# Import results
comblist  <- readRDS("results/processing/comblist.RDS") %>% print()
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()

cols_log2fc <- str_subset(names(combdf.gd), "_log2fc")
cols_qval   <- str_subset(names(combdf.gd), "_q")


# Convert to long format
gd <- combdf.gd %>%
        select(guide, geneSys, gene, set, essential, all_of(cols_log2fc), 
               all_of(cols_qval), FDR0.05_count) %>%
        pivot_longer(cols = -c(guide, set, essential, geneSys, gene, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>% 
        mutate(sign = ifelse(log2fc > 0, "positive", "negative")) %>%
        mutate(specificity = case_when(FDR0.05_count < 3 ~ "specific",
                                       FDR0.05_count > 2 ~ "nonspecific",
                                       TRUE ~ NA_character_)) %>%
        filter(q < 0.05) %>% print()

gn <- combdf.gn %>%
        select(-geneDescr) %>%
        pivot_longer(cols = -c(geneSys, gene, essential, FDR0.05_count),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>% 
        mutate(sign = ifelse(log2fc > 0, "positive", "negative")) %>%
        mutate(specificity = case_when(FDR0.05_count < 3 ~ "specific",
                                       FDR0.05_count > 2 ~ "nonspecific",
                                       TRUE ~ NA_character_)) %>%
        filter(q < 0.05) %>% print()




# Stats on guide-level ---------------------------------------------------------

# Total number of guides with effect
gd %>%
        filter(FDR0.05_count > 0) %>%
        select(guide) %>%
        unique() %>%
        nrow()        # 1020


# Mean number of significant guides per protein
gd %>%        
        group_by(protein) %>% 
        summarize(n = n()) %>% #print()
        pull(n) %>% 
        mean()     # 191


# Mean number of negative/positive guides per protein
gd %>%
        group_by(protein, sign) %>%
        summarize(n = n()) %>% #print()
        group_by(sign) %>%        
        summarize(mean = mean(n))
# negative    134.
# positive     57 


# Mean number of specific/nonspecific guides per protein
gd %>%
        group_by(protein, specificity) %>%
        summarize(n = n()) %>% #print()
        group_by(specificity) %>%        
        summarize(mean = mean(n))
# nonspecific  106.
# specific      84.8


# Mean number of guides per protein by set
gd %>%
        group_by(protein, set) %>%
        summarize(n = n()) %>% #print(n=50)
        group_by(set) %>%        
        summarize(mean = mean(n))
# eProvs  126. 
# neProvs  24.2
# neStops  40.5


# Mean number of guides per protein by set and sign
gd %>%
        group_by(protein, sign, set) %>%
        summarize(n = n()) %>% #print()
        group_by(sign, set) %>%        
        summarize(mean = mean(n))
# negative eProvs  103. 
# negative neProvs  10.1
# negative neStops  20.4
# positive eProvs   22.7
# positive neProvs  14.1
# positive neStops  20.2


# Mean number of guides by set
gd.sum <- combdf.gd %>%
        mutate(sig = ifelse(FDR0.05_count > 0, T, F)) %>%
        group_by(set, sig) %>%
        summarize(n = n()) %>%
        pivot_wider(names_from = set, values_from = n) %>%
        print()
# sig   eProvs neProvs neStops
# FALSE   4861    5317    5152
# TRUE     582     157     281

chisq.test(gd.sum[, c("neProvs", "neStops")])
# X-squared = 36.957, df = 1, p-value = 1.208e-09


# Mean log2fc by set and sign
gd %>%
        group_by(set, sign) %>%
        summarize(mean = mean(log2fc),
                  median = median(log2fc))
# set     sign      mean median
# eProvs  negative -1.18 -1.13 
# eProvs  positive  1.15  0.969
# neProvs negative -1.06 -0.927
# neProvs positive  1.09  0.970
# neStops negative -1.16 -1.12 
# neStops positive  1.19  1.12 


# Mean abs(log2fc) by set
gd %>%
        group_by(set) %>%
        summarize(mean = mean(abs(log2fc)),
                  median = median(abs(log2fc)))
# set      mean median
# eProvs   1.17  1.11 
# neProvs  1.08  0.960
# neStops  1.18  1.12 


# t-test for abs(log2fc) by set
neProvs.abs <- gd %>%
        filter(set == "neProvs") %>%
        pull(log2fc) %>% abs() %>% print()

neStops.abs <- gd %>% 
        filter(set == "neStops") %>% 
        pull(log2fc) %>% abs() %>% print()

t.test(neProvs.abs, neStops.abs, alternative = "two.sided")
# p-value = 0.002807




# Stats on gene-level ----------------------------------------------------------

# Total number of genes with effect
gn %>%
        filter(FDR0.05_count > 0) %>%
        select(geneSys) %>%
        unique() %>%
        nrow()      # 710


# Mean number of significant genes per protein
gn %>%        
        group_by(protein) %>%
        summarize(n = n()) %>% #print()
        pull(n) %>% 
        mean()      # 156


# Mean number of positive/negative genes per protein
gn %>%
        group_by(protein, sign) %>%
        summarize(n = n()) %>% #print()
        group_by(sign) %>%        
        summarize(mean = mean(n))
# negative    106.
# positive     50 

binom.test(x = 106, n = (106 + 50), p = 0.5) 
# number of successes = 106, number of trials = 156, p-value = 8.633e-06
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval: 0.6001447 0.7518587
# sample estimates: probability of success: # 0.6794872 


# Mean number of positive/negative genes per protein, by protein fct
gn %>%
        mutate(proteinfct = ifelse(protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1"), "stress", "core")) %>%
        group_by(proteinfct, protein, sign) %>%
        summarize(n = n()) %>% #print()
        group_by(proteinfct, sign) %>%       
        summarize(mean = mean(n))

# proteinfct sign      mean
# core       negative 119. 
# core       positive  33.7
# stress     negative  83.5
# stress     positive  78.5

binom.test(x = 119, n = (119 + 34), p = 0.5)
# number of successes = 119, number of trials = 153, p-value = 2.963e-12
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval: 0.7035548 0.8409278
# sample estimates: probability of success: 0.7777778 

binom.test(x = 84, n = (84 + 79), p = 0.5)
# number of successes = 84, number of trials = 163, p-value = 0.7542
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval: 0.4358847 0.5942224
# sample estimates: probability of success: 0.5153374 


# Mean number of specific/nonspecific genes per protein
gn %>%
        group_by(protein, specificity) %>%
        summarize(n = n()) %>% #print()
        group_by(specificity) %>%        
        summarize(mean = mean(n))
# nonspecific  103.
# specific      53.5


# Mean number of essential/non-essential genes per protein
gn %>%
        group_by(protein, essential) %>%
        summarize(n = n()) %>% #print(n=50)
        group_by(essential) %>%        
        summarize(mean = mean(n))
# essential  mean
# FALSE      59.6
# TRUE       96.5

binom.test(x = 97, n = (60 + 97), p = sum(combdf.gn$essential)/nrow(combdf.gn))
# p-value < 2.2e-16 

table(combdf.gn$essential, ifelse(combdf.gn$FDR0.05_count > 0, "sig", "nonsig"))
#       nonsig  sig
# FALSE   3202  366
# TRUE     600  344

chisq.test(combdf.gn$essential, ifelse(combdf.gn$FDR0.05_count > 0, "sig", "nonsig"))
# X-squared = 383.98, df = 1, p-value < 2.2e-16


# Mean number of essential/non-essential genes per protein, by sign
gn %>%
        group_by(protein, sign, essential) %>%
        summarize(n = n()) %>% #print()
        group_by(sign, essential) %>%        
        summarize(mean = mean(n))
# sign     essential  mean
# negative FALSE      29.4
# negative TRUE       76.8
# positive FALSE      30.3
# positive TRUE       19.7

76.8/(76.8+29.4) # 72% of negative regulators are essential
19.7/(19.7+30.3) # 39% of positive regulators are essential
76.8/(76.8+19.7) # 80% of essential regulators have negative effect
29.4/(29.4+30.3) # 49% of nonessential regulators have negative effect




# Plot counts of significant guides/genes --------------------------------------

(pgd1 <- gd %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2",
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, sign) %>%
        summarize(n = n()) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n),
               sign = case_when(sign == "negative" ~ "Negative",
                                sign == "positive" ~ "Positive",
                                TRUE ~ NA_character_)) %>%
        ggplot(aes(x = protein, y = n, fill = sign)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of guides") +
        scale_y_continuous(limits = c(-300, 150)) +
        scale_fill_manual(values = c("#B2ABD2", "#FDB863")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
ggsave(paste0(resdir, "Nsig_gd_fill.pdf"), width = 5, height = 4)


(pgn1 <- gn %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, sign) %>%
        summarize(n = n()) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n),
               sign = case_when(sign == "negative" ~ "Negative",
                                sign == "positive" ~ "Positive",
                                TRUE ~ NA_character_)) %>%
        ggplot(aes(x = protein, y = n, fill = sign)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of genes") +
        scale_y_continuous(limits = c(-225, 125)) +
        scale_fill_manual(values = c("#B2ABD2", "#FDB863")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
ggsave(paste0(resdir, "Nsig_gn_fill.pdf"), width = 5, height = 4)




# Plot counts of significant guides/genes by specificity -----------------------

(pgd2 <- gd %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2",
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, sign, specificity) %>%
        summarize(n = n()) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n)) %>%
        mutate(specificity = case_when(specificity == "nonspecific" ~ "Nonspecific",
                                       specificity == "specific" ~ "Specific",
                                       TRUE ~ NA_character_)) %>%
        mutate(specificity = forcats::fct_rev(as.factor(specificity))) %>%
        ggplot(aes(x = protein, y = n, fill = specificity)) +
        #geom_bar(stat = "identity", position = "dodge") +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of guides") +
        scale_y_continuous(limits = c(-300, 150)) +
        scale_fill_manual(values = c("lightskyblue3", "lightskyblue4")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
#ggsave(paste0(resdir, "Nsig_spec_gd_dodge.pdf"), width = 6, height = 4)
ggsave(paste0(resdir, "Nsig_spec_gd_fill.pdf"), width = 5, height = 4)


(pgn2 <- gn %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, sign, specificity) %>%
        summarize(n = n()) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n)) %>%
        mutate(specificity = case_when(specificity == "nonspecific" ~ "Nonspecific",
                                       specificity == "specific" ~ "Specific",
                                       TRUE ~ NA_character_)) %>%
        mutate(specificity = forcats::fct_rev(as.factor(specificity))) %>%
        ggplot(aes(x = protein, y = n, fill = specificity)) +
        #geom_bar(stat = "identity", position = "dodge") +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of genes") +
        scale_y_continuous(limits = c(-225, 125)) +
        scale_fill_manual(values = c("lightskyblue3", "lightskyblue4")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
#ggsave(paste0(resdir, "Nsig_spec_gn_dodge.pdf"), width = 6, height = 4)
ggsave(paste0(resdir, "Nsig_spec_gn_fill.pdf"), width = 5, height = 4)




# Plot counts of significant guides/genes by essentiality -------------------

(pgd3 <- gd %>%
         filter(set %in% c("neProvs", "eProvs")) %>%
         mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
                essential = forcats::fct_rev(as.factor(essential))) %>%
         mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                                "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                 protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                 TRUE ~ NA_character_)) %>%
         mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                     "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                     "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
         group_by(type, protein, essential, sign) %>%
         summarize(n = n()) %>%
         mutate(n = ifelse(sign == "negative", -1*n, n)) %>%
         ggplot(aes(x = protein, y = n, fill = essential)) +
         #geom_bar(stat = "identity", position = "dodge") +
         geom_bar(stat = "identity") +
         geom_hline(yintercept = 0, color = "white") +
         ylab("Number of guides") +
         scale_y_continuous(limits = c(-250, 100)) +
         scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               legend.title = element_blank(),
               axis.text.x = element_text(angle = 30, vjust = 0.6), 
               axis.title.x = element_blank()) +
         facet_grid(~ type, space = "free", scales = "free_x"))
#ggsave(paste0(resdir, "Nsig_ess_gd-provsOnly_dodge.pdf"), width = 6, height = 4)
ggsave(paste0(resdir, "Nsig_ess_gd-provsOnly_fill.pdf"), width = 5, height = 4)


(pgn3 <- gn %>%
        mutate(essential = ifelse(essential == TRUE, "Essential", "Nonessential"),
               essential = forcats::fct_rev(as.factor(essential))) %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, essential, sign) %>%
        summarize(n = n()) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n)) %>% 
        ggplot(aes(x = protein, y = n, fill = essential)) +
        #geom_bar(stat = "identity", position = "dodge") +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of genes") +
        scale_y_continuous(limits = c(-250, 125)) +
        scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
#ggsave(paste0(resdir, "Nsig_ess_gn_dodge.pdf"), width = 6, height = 4)
ggsave(paste0(resdir, "Nsig_ess_gn_fill.pdf"), width = 5, height = 4)




# Plot counts of significant guides by set -------------------------------------

(pgd4 <- gd %>%
        mutate(set = case_when(set == "eProvs"  ~ "Conserved residue,\nessential gene",
                               set == "neProvs" ~ "Conserved residue,\nnonessential gene",
                               set == "neStops" ~ "Premature stop,\nnonessential gene",
                               TRUE ~ NA_character_)) %>%
        mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                               "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                TRUE ~ NA_character_)) %>%
        mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                    "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                    "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
        group_by(type, protein, set, sign) %>%
        summarize(n = n()) %>%
        mutate(set = as.factor(set),
              set = forcats::fct_rev(set)) %>%
        mutate(n = ifelse(sign == "negative", -1*n, n)) %>% print() %>%
        ggplot(aes(x = protein, y = n, fill = set)) +
        #geom_bar(stat = "identity", position = "dodge") +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, color = "white") +
        ylab("Number of gRNAs") +
        scale_y_continuous(limits = c(-300, 150)) +
        scale_fill_manual(values = c("lightsalmon2", "lightsalmon3", "lightsalmon4")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              #legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(margin = margin(t = 3, b = 3, unit = "pt")),
              axis.text.x = element_text(angle = 30, vjust = 0.6), 
              axis.title.x = element_blank()) +
        facet_grid(~ type, space = "free", scales = "free_x"))
#ggsave(paste0(resdir, "Nsig_set_gd_dodge.pdf"), width = 6.5, height = 4)
ggsave(paste0(resdir, "Nsig_set_gd_fill.pdf"), width = 5.5, height = 4)




# All counts plots combined ----------------------------------------------------

all <- egg::ggarrange(pgn1, pgn2, pgn3, pgd4, ncol = 1)
ggsave(plot = all, paste0(resdir, "Nsig_all.pdf"), width = 5.5, height = 12)




# Plot log2fc of significant guides by set -------------------------------------

(pgd5 <- gd %>%
         mutate(set = case_when(set == "eProvs"  ~ "Conserved residue,\nessential gene",
                                set == "neProvs" ~ "Conserved residue,\nnonessential gene",
                                set == "neStops" ~ "Premature stop,\nnonessential gene",
                                TRUE ~ NA_character_)) %>%
         mutate(set = forcats::fct_rev(as.factor(set))) %>%
         mutate(type = case_when(protein %in% c("Eno2", "Fas1", "Fas2", 
                                                "Htb2", "Rnr2", "Rpl9A", "Tdh3") ~ "Core functions",
                                 protein %in% c("Tdh1", "Tdh2", "Ssa1", "Yhb1") ~ "Stress regulated",
                                 TRUE ~ NA_character_)) %>%
         mutate(protein = factor(protein, levels = c("Eno2", "Fas1", "Fas2",
                                                     "Htb2", "Rnr2", "Rpl9A", "Tdh3",
                                                     "Tdh1", "Tdh2", "Ssa1", "Yhb1"))) %>%
         ggplot(aes(x = set, y = log2fc, fill = set)) +
         geom_violin() +
         ylab("Log2FC") +
         scale_fill_manual(values = c("lightsalmon2", "lightsalmon3", "lightsalmon4")) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               legend.position = "none",
               legend.title = element_blank(),
               axis.text.x = element_text(angle = 30, vjust = 0.6), 
               axis.title.x = element_blank()))
ggsave(paste0(resdir, "Eff_set_gd.pdf"), width = 3, height = 3.5)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-07_SessionInfo.txt")

