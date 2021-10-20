# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/overlappingGuides/"
dir.create(resdir)


# Import results
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()


# Guide table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set %in% c("eProvs", "neProvs", "neStops")) %>% print()




# Identify gRNAs introducing the same mutations --------------------------------

mutations2 <- guides %>%
        select(guide, geneSys, gene, prot.start, mut1, mut2, mut3) %>%
        arrange(geneSys, prot.start) %>%
        group_by(geneSys) %>%
        mutate(locDiff = prot.start - lag(prot.start)) %>%
        mutate(sel1 = ifelse(mut1 == lag(mut1) |
                             mut1 == lag(mut2) |
                             mut1 == lag(mut3) |
                             mut2 == lag(mut1) |
                             mut2 == lag(mut2) |
                             mut2 == lag(mut3) |
                             mut3 == lag(mut1) |
                             mut3 == lag(mut2) |
                             mut3 == lag(mut3), T, F)) %>%
        mutate(sel2 = ifelse(lead(sel1) == T, T, F)) %>%
        filter(sel1 == T | sel2 == T) %>% 
        ungroup() %>% print()

sel1 <- mutations2 %>%
        filter(sel2 == T) %>%
        mutate(idx1 = row_number()) %>%
        select(guide, idx1) %>% print()

sel2 <- mutations2 %>%
        filter(sel1 == T) %>%
        mutate(idx2 = row_number()) %>%
        select(guide, idx2) %>% print()



# Plot log2fc from the two guides targeting the same mutation
combdf.gd %>%
        select(-c(PAMstrand:essential, geneDescr, FDR0.05_count)) %>%
        pivot_longer(cols = -c(guide, set, geneSys, gene),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        full_join(sel1, by = "guide") %>%
        full_join(sel2, by = "guide") %>%
        mutate(xy = ifelse(!is.na(idx1), "x", ifelse(!is.na(idx2), "y", NA)),
               id = ifelse(!is.na(idx1), idx1, idx2)) %>%
        filter(!is.na(xy)) %>% 
        select(-idx1, -idx2) %>%
        pivot_wider(id_cols = -c(geneSys, gene),
                    names_from = xy,
                    values_from = c(guide, set, log2fc, q)) %>%
        filter(!is.na(q_x) & !is.na(q_y)) %>%
        mutate(qmin = ifelse(q_x < q_y, q_x, q_y),
               qmax = ifelse(q_x > q_y, q_x, q_y)) %>%
        filter(qmin < 0.05) %>%
        #filter(qmax < 0.05) %>%
        #pull(id) %>% length()                # to get the number of data points
        #pull(id) %>% unique() %>% length()   # to get the number of guide pairs
        ggplot(aes(x = log2fc_x, y = log2fc_y)) +
        geom_smooth(method = "lm") +
        ggpubr::stat_cor(color = "grey25", 
                         r.digits = 3, 
                         label.x.npc = "left", 
                         label.y.npc = "top", 
                         label.sep = "\n") +
        geom_point() + 
        labs(x = "gRNA 1  Log2FC", y = "gRNA 2  Log2FC") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "OverlappingGuides_scatter_qmin0.05.pdf"), width = 3, height = 3)       
ggsave(paste0(resdir, "OverlappingGuides_scatter_qmax0.05.pdf"), width = 3, height = 3)       




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-06_SessionInfo.txt")

