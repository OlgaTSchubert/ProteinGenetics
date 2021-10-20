# Libraries ---------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/general_Cytoscape/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS") %>% print()




# Reformat results for Cytoscape -----------------------------------------------

combdf.gn.cy <- combdf.gn %>%
        select(-geneDescr, -essential, -starts_with("FDR")) %>%
        pivot_longer(cols = -c("geneSys", "gene"),
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        filter(q < 0.05) %>%
        relocate(protein, starts_with("gene")) %>%
        mutate(interaction = "cr") %>%
        group_by(geneSys) %>%
        mutate(intcountnode = n(),
               intcountedge = n(),
               intcountrev = 11-n()) %>% print()
write_delim(combdf.gn.cy, paste0(resdir, "combdf_gn_cytoscape.tsv"), delim = "\t")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-08_SessionInfo.txt")

