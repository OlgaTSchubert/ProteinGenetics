# Libraries --------------------------------------------------------------------

library(tidyverse)
library(heatmaply)
library(GGally)
library(corrplot)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Results directory
resdir <- "results/correlations/"
dir.create(resdir)


# Import results
comblist  <- readRDS("results/processing/comblist.RDS") %>% print()
combdf.gd <- readRDS("results/processing/combdf_gd.RDS") %>% print()




# Replicate correlations (heatmaply with dendrogram) ---------------------------

# For pdf file saving to work, needed orca (installed through conda)

for (i in 1:length(comblist)) {
        print(i)
        readcols <- str_which(names(comblist[[i]]), "^..\\..$")
        cr <- cor(comblist[[i]][, readcols], method = "pearson", use = "pairwise.complete.obs")
        heatmaply_cor(cr, 
                      limits = c(0, 1), 
                      main = names(comblist)[i],
                      hclustfun = hclust, 
                      hclust_method = "ward.D",
                      plot_method = "ggplot", 
                      dendrogram = "both", 
                      show_dendrogram = c(T, F),
                      file = paste0(resdir, names(comblist)[i],"_replicates.pdf"))
}
dev.off()




# Experiment correlations (ggpairs) --------------------------------------------

# Table with p-values
gd.p <- select(combdf.gd, str_which(names(combdf.gd), "_q")) %>% print()


# Table with ratios
gd.r <- select(combdf.gd, str_which(names(combdf.gd), "_log2fc"))
gd.r <- replace(gd.r, abs(gd.r) > 5, NA_real_) %>%
        rename_all(list(~ str_extract(., ".*(?=_log2fc)"))) %>% print()


# Table of significant ratios only, rest filled with NA
gd.rsig <- gd.r %>% replace(gd.p > 0.05, NA_real_) %>% print()


# Pairwise correlations
pdf(paste0(resdir, "all_proteins_sigGuides_ggpairs.pdf"))
ggpairs(gd.rsig, 
        axisLabels = "none", 
        switch = "y",
        upper = list(continuous = wrap("smooth", alpha = 0.3, size = 0.2)),
        lower = list(continuous = wrap("points", alpha = 0.3, size = 0.2))) + 
        theme(panel.grid.major = element_blank(), 
              panel.background = element_rect(fill = "grey95"))
dev.off()




# Corrplots & corr-p-value histograms ------------------------------------------

correlationPlot <- function(input, filename) {
        
        gd.corrs <- cor(input, use = "pairwise.complete.obs", method = "pearson")
        gd.corrs %>% as_tibble(rownames = "Protein")
        write_csv(as.data.frame(gd.corrs), paste0(resdir, filename, ".csv"))
        
        gd.corrp <- Hmisc::rcorr(as.matrix(input), type = "pearson")$P
        
        pdf(paste0(resdir, filename, ".pdf"))
        corrplot(gd.corrs, p.mat = gd.corrp, type = "lower", diag = F, 
                 insig = "p-value", sig.level = "-1")
        hist(gd.corrs, xlim = c(-1, 1), breaks = 20)
        hist(gd.corrp, xlim = c(0, 1), breaks = 20)
        dev.off()
        
}

correlationPlot(gd.rsig, "all_proteins_sigGuides_corrplot")




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-05_SessionInfo.txt")

