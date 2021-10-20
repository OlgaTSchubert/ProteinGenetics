# Libraries --------------------------------------------------------------------

library(tidyverse)
library(qvalue)
library(poolr)
library(GGally)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")
source("code/BE19-00_Functions.R")


# Reads directory
readsdir <- "~/BigDataFiles/2021_ProteinGenetics/BE19/ReadsRDS/"


# Data directory
datadir <- "data/"
dir.create(datadir)


# Results directory
resdir <- "results/"
dir.create(resdir)


# Processing directory
procdir <- "results/processing/"
dir.create(procdir)


# Copy read counts files into local data directory (datadir)
countsfiles0 <- list.files(readsdir, pattern = ".*_ReadCounts.RDS", full.names = T); countsfiles0
for(f in seq_along(countsfiles0)) {
        file.copy(from = countsfiles0[f], to = datadir, overwrite = F)
}
rm(countsfiles0)


# List of input files
countsfiles <- list.files(datadir, full.names = T); countsfiles


# List of experiment names (proteins)
experiments <- countsfiles %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_ReadCounts.RDS)") %>%
        str_sort() %>%
        print()


# Guide table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set %in% c("eProvs", "neProvs", "neStops")) %>%
        select(c("guide", "geneSys", "gene", "set")) %>%
        print()


# Gene annotations from SGD
sgd <- importExtData(dataset = "SGD_features", localfile = T) %>% 
        select(c("Feature_name", "Standard_gene_name", "Description")) %>% 
        rename(geneSys   = Feature_name,
               gene      = Standard_gene_name,
               geneDescr = Description) %>%
        print()




# Replicates to be excluded from analysis --------------------------------------

badreps <- tibble(experiment = c("Eno2", "Fas1", "Fas2", "Htb2", "Rnr2", 
                                 "Rpl9A", "Ssa1", "Tdh1", "Tdh2", "Tdh3", "Yhb1"),
                  replicate = list("A", c("A", "F"), c("A", "D"), "A", NA, 
                                   NA, NA, NA, NA, NA, NA)) #%>% print()




# Differential analysis of BC counts -------------------------------------------

for(exp in seq_along(experiments)) {

        print(experiments[exp])
        
        
        # Create results directory for experiment
        diffdir <- paste0(procdir, experiments[exp], "_processing/")
        dir.create(diffdir)
        
        
        # Convert read count table to long format
        rcounts <- readCountsLong(counts.file = countsfiles[exp], guide.info = guides)


        # Get unique BC counts per guide and sample
        BCcountsX <- rcounts %>%
                filter(n_reads > 1) %>%  # Only keep BCs with >1 read
                group_by(sample, guide) %>%
                summarize(nBCs = n()) %>%
                ungroup() #%>% print()

        guidesToKeep <- BCcountsX %>%
                filter(str_detect(sample, "un")) %>%
                group_by(guide) %>%
                summarize(moreBCs = sum(nBCs > 1)) %>%  # Keep only guides with >1 BC
                filter(moreBCs > 4) #%>% print()        # in more than 4 "un" samples

        BCcounts <- BCcountsX %>%
                filter(guide %in% guidesToKeep$guide) %>%
                pivot_wider(names_from = sample, values_from = nBCs ) %>%
                mutate_all(list(~ replace_na(., 0))) %>%
                pivot_longer(cols = -guide, names_to = "sample", values_to = "nBCs") %>%
                separate(sample, into = c("tail", "repl"), remove = F) #%>% print()
        
        
        # Remove replicates specified above to be excluded
        excl <- badreps %>% filter(experiment == experiments[exp]) %>% 
                pull("replicate") %>% .[[1]]
        if(!NA %in% excl) { BCcounts <- BCcounts[!BCcounts$repl %in% excl, ] }

        
        # Save BCcounts        
        saveRDS(BCcounts, paste0(diffdir, "BCcounts.RDS"))
        
        
        
        
        # Differential analysis ------------------------------------------------
        
        # Check for overdispersion of counts (mean vs variance scatter plot)
        BCcounts %>%
                group_by(guide, tail) %>%
                summarize(ave = mean(nBCs), 
                          vari = (sd(nBCs))^2) %>%
                ggplot(aes(x = ave, y = vari)) +
                geom_point() +
                #geom_hex() +
                scale_x_continuous(trans = "log10") +
                scale_y_continuous(trans = "log10") +
                labs(x = "Average BC count", y = "Variance",
                     title = "Variance by average BC counts (across replicates)",
                     subtitle = "Each point corresponds to a guide in a group (hi, lo, un).")
        ggsave(paste0(diffdir, "BCcounts_mean-var.pdf"))
        
        
        # Run Poisson generalized linear model (GLM) with log link function
        BCstats <- glmStats(BCcounts)


        # Convert to tibble
        BCstats <- BCstats[[1]] %>% as_tibble()
        
        
        # Convert the fold-change estimate from ln to log2
        BCstats <- BCstats %>% mutate(log2fc = log2(exp(1)^estimate))
        
        
        # For guides where no barcode was detected in any sample of one tail
        # the estimate is ~20/-20. Replace these values with 3/-3.
        BCstats <- BCstats %>% mutate(log2fc = case_when(log2fc < -5 ~ -3,
                                                         log2fc >  5 ~  3,
                                                         TRUE ~ log2fc))

        # Get q-values
        BCstats$q.value <- qvalue(BCstats$p.value)$qvalue
        
        
        saveRDS(BCstats, paste0(diffdir, "BCstats.RDS"))

}




# Further processing of BC stats -----------------------------------------------

for(exp in seq_along(experiments)) {
        
        print(experiments[exp])
        
        # Read in the BCstats.RDS file
        diffdir  <- paste0(procdir, experiments[exp], "_processing/")
        BCstats  <- readRDS(paste0(diffdir, "BCstats.RDS"))
        BCcounts <- readRDS(paste0(diffdir, "BCcounts.RDS"))
        
        # Generate wide version of BCcounts, needed at the end
        BCcounts.wide <- BCcounts %>%
                pivot_wider(id_cols = "guide", names_from = sample, values_from = nBCs) #%>% print()
        
        
        # Plot p- and q-value distributions ------------------------------------
        
        # Plot p-value distribution
        pdf(paste0(diffdir, "pvalueHist.pdf"), width = 4, height = 4)
        hist(BCstats$p.value, breaks = 100)
        dev.off()


        # Plot volcano plot
        pdf(paste0(diffdir, "volcanoPlots.pdf"), width = 4, height = 4)
        plot(BCstats$log2fc, -log10(BCstats$p.value))
        dev.off()


        # Generate diagnostic qvalue plots
        pdf(paste0(diffdir, "qvaluePlots.pdf"), width = 5, height = 5)
        plot(qvalue(BCstats$p.value))
        dev.off()


        # Plot q-value distribution - not clear what to expect here
        pdf(paste0(diffdir, "qvalueHist.pdf"), width = 4, height = 4)
        hist(BCstats$q.value, breaks = 100)
        dev.off()
        
        
        # Check if p/q-values depend on number of BCs per guide
        BCstats.nBCs <- BCcounts %>%
                group_by(guide) %>%
                summarize(gnBCs = sum(nBCs)) %>%
                left_join(BCstats, by = c("guide" = "gseq")) %>% print()
        
        BCstats.nBCs %>%
                ggplot(aes(x = gnBCs, y = p.value)) +
                geom_point() +
                scale_x_continuous(trans = "log10") +
                labs(x = "BCs per guide", y = "P-value (Hi/Lo)")
        ggsave(paste0(diffdir, "pvalue-nBC.pdf"), width = 4, height = 4)
        
        BCstats.nBCs %>%
                ggplot(aes(x = gnBCs, y = q.value)) +
                geom_point() +
                scale_x_continuous(trans = "log10") +
                labs(x = "BCs per guide", y = "Q-value (Hi/Lo)")
        ggsave(paste0(diffdir, "qvalue-nBC.pdf"), width = 4, height = 4)

        
        
        # Get gene-level stats -------------------------------------------------
        
        geneStats <- BCstats %>%
                left_join(guides, by = c("gseq" = "guide")) %>%
                group_by(geneSys) %>%
                summarize(guidesPerGene   = n(),
                          log2fc.max      = max(log2fc),
                          q.of.log2fc.max = q.value[which.max(abs(log2fc))],
                          q.min           = min(q.value),
                          log2fc.of.q.min = log2fc[which.min(q.value)],
                          p.fisher        = poolr::fisher(p.value)$p) %>%
                mutate(q.fisher = qvalue(p.fisher)$qvalue) %>%
                mutate(q.comb   = ifelse(q.min < q.fisher, 
                                         q.min, q.fisher)) %>%
                select(-p.fisher) %>%
                relocate(log2fc.of.q.min, .before = q.min) %>%
                arrange(q.comb) %>% 
                print()
        
        
        # Compare different types of p-value combinations        
        pl <- geneStats %>%
                select(starts_with("q.")) %>%
                ggpairs(aes(alpha = 0.3))
        ggsave(paste0(diffdir, "pvalue-comb_comparison.pdf"), plot = pl)

        
        
        # Compile final guide results table ------------------------------------
        
        # Cleanup BCstats
        BCstats <- BCstats %>%
                select(gseq, log2fc, q.value) %>%
                rename(guide = gseq,
                       q     = q.value) %>% print()
        
        # Combine all the results into a single table
        BCstats.wide <- guides %>%
                inner_join(BCcounts.wide, by = "guide") %>%
                left_join(BCstats, by = "guide") %>%
                left_join(geneStats, by = "geneSys") %>% 
                relocate(c(set, guidesPerGene, geneSys, gene), .after = last_col()) %>%
                left_join(sgd, by = c("geneSys", "gene")) %>% 
                print()
        saveRDS(BCstats.wide, paste0(diffdir, "BCstats.wide.RDS"))
        write_excel_csv(BCstats.wide, paste0(diffdir, "BCstats.csv"))
}




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-03_SessionInfo.txt")

