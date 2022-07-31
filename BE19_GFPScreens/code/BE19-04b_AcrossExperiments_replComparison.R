# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/processing_replComparison/"


# Processing directory (4-replicate analysis)
procdir <- "results/processing_replComparison/"


# List of all results tables (BCstats.csv)
BCstats.files <- list.files(path = procdir, pattern = "BCstats.csv", 
                            recursive = T, full.names = T); BCstats.files


# List of experiment names (proteins)
experiments <- BCstats.files %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_processing_..../BCstats.csv)") %>%
        unique() %>% print()


# Import all results tables and store in a list
comblist.ABCD <- BCstats.files %>%
        .[matches(".*ABCD.*", vars = .)] %>%
        map(~ read_csv(.))
names(comblist.ABCD) <- experiments

comblist.EFGH <- BCstats.files %>%
        .[matches(".*EFGH.*", vars = .)] %>%
        map(~ read_csv(.))
names(comblist.EFGH) <- experiments


# Import results (8-replicate analysis)
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()




# Combine guide results --------------------------------------------------------

# Combine results into a table
combdf.gd.ABCD <- map(comblist.ABCD, ~ select(., guide, log2fc, q)) %>%
        purrr::reduce(full_join, by = "guide") %>% print()

combdf.gd.EFGH <- map(comblist.EFGH, ~ select(., guide, log2fc, q)) %>%
        purrr::reduce(full_join, by = "guide") %>% print()


# Add protein names to column headers
nam1.ABCD  <- rep(experiments, each = 2)
nam2.ABCD  <- rep(c("log2fc", "q"), times = length(experiments))
names(combdf.gd.ABCD)[2:(1+length(experiments)*2)] <- paste0(nam1.ABCD, "_", nam2.ABCD)

nam1.EFGH  <- rep(experiments, each = 2)
nam2.EFGH  <- rep(c("log2fc", "q"), times = length(experiments))
names(combdf.gd.EFGH)[2:(1+length(experiments)*2)] <- paste0(nam1.EFGH, "_", nam2.EFGH)


# Combine all guide results into  one table
comb <- combdf.gd %>%
        select(guide, Eno2_log2fc:Yhb1_q) %>% 
        left_join(combdf.gd.ABCD, by = c("guide"), suffix = c("_all", "_ABCD")) %>% 
        left_join(combdf.gd.EFGH, by = c("guide")) %>% 
        rename_with(~ str_replace(., "_log2fc$", "_log2fc_EFGH")) %>%
        rename_with(~ str_replace(., "_q$", "_q_EFGH")) %>% 
        pivot_longer(cols = -guide, 
                     names_to = c("protein", ".value", "replicates"),
                     names_pattern = "(.*)_(.*)_(.*)") %>% 
        pivot_wider(id_cols = c(guide, protein),
                    names_from = replicates,
                    values_from = c(log2fc, q)) %>%
        print()




# Number of sign flips ---------------------------------------------------------

comb %>%
        mutate(signflip = ifelse(log2fc_ABCD*log2fc_EFGH < 0, T, F)) %>%
        filter(q_all < 0.05 & !is.na(signflip)) %>%
        summarize(fraction = sum(signflip)/n())

# fraction = 0.236 (23.6%)


comb %>%
        mutate(signflip = ifelse(log2fc_ABCD*log2fc_EFGH < 0, T, F)) %>%
        filter(q_all < 0.05 & !is.na(signflip)) %>%
        group_by(protein) %>%
        summarize(fraction = sum(signflip)/n())

# protein  n
# Eno2     0.00444
# Fas1     0      
# Fas2     0.182  
# Htb2     0.199  
# Rnr2     0.175  
# Rpl9A    0.275  
# Ssa1     0.292  
# Tdh1     0.485  
# Tdh2     0.417  
# Tdh3     0.163  
# Yhb1     0.393  




# Correlation plots ------------------------------------------------------------

comb %>%
        filter(q_all < 0.05) %>% 
        #filter(q_ABCD < 0.05 & q_EFGH < 0.05) %>%
        #filter(q_ABCD < 0.05 | q_EFGH < 0.05) %>%
        ggplot(aes(x = log2fc_ABCD, y = log2fc_EFGH)) +
        geom_hline(yintercept = 0, alpha = 0.3) +
        geom_vline(xintercept = 0, alpha = 0.3) +
        geom_smooth(method = "lm", color = "black",
                    formula = y ~ x-1) +  # force through origin
        ggpubr::stat_cor(color = "grey25",               # numbers not for line forced through origin
                         r.digits = 3, 
                         label.x.npc = "left", 
                         label.y.npc = "top", 
                         label.sep = "\n") +
        geom_point(alpha = 0.2) +
        labs(x = "Log2FC for replicates A, B, C, D",
             y = "Log2FC for replicates E, F, G, H") +
        lims(x = c(-4, 4), y = c(-4, 4)) +
        #facet_wrap(~ protein) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())

ggsave(paste0(resdir, "ABCD-vs-EFGH_rep8-q0.05.jpg"), width = 4, height = 4)
ggsave(paste0(resdir, "ABCD-vs-EFGH_rep8-q0.05.pdf"), width = 4, height = 4)
ggsave(paste0(resdir, "ABCD-vs-EFGH_rep8-q0.05_facet.jpg"), width = 7, height = 6)
ggsave(paste0(resdir, "ABCD-vs-EFGH_rep8-q0.05_facet.pdf"), width = 7, height = 6)
#ggsave(paste0(resdir, "ABCD-vs-EFGH_either-q0.05.pdf"), width = 4, height = 4)
#ggsave(paste0(resdir, "ABCD-vs-EFGH_both-q0.05.pdf"), width = 4, height = 4)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-04b_SessionInfo.txt")

