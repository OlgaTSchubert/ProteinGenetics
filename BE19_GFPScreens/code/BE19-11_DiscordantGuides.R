# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/discordantGuides/"
dir.create(resdir)


# Import results
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()




# Identify discordant guides ---------------------------------------------------
# Check all guides with q < 0.1

discr <- combdf.gd %>%
        select(-FDR0.05_count) %>%
        pivot_longer(cols = ends_with("_log2fc") | ends_with("_q"), 
                     names_to = c("protein", ".value"), 
                     names_pattern = ("(.*)_(.*)")) %>%
        filter(!is.na(log2fc)) %>%
        filter(q < 0.1) %>%
        group_by(geneSys, protein) %>%
        mutate(nneg = sum(log2fc < 0),
               npos = sum(log2fc > 0)) %>%
        ungroup() %>%
        filter((npos != 0) & (nneg != 0)) %>%
        arrange(gene, protein) %>%
        print(n=50)
write_csv(discr, paste0(resdir, "discordant_q0.1.csv"))


combdf.gd %>%
        filter(geneSys %in% discr$geneSys) %>% print() %>%
        write_csv(paste0(resdir, "discordant_q0.1_big.csv"))




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-11_SessionInfo.txt")

