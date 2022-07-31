# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Reads directory
readsdir <- "~/BigDataFiles/2021_ProteinGenetics/BE19/ReadsRDS/"


# Data directory
datadir <- "data/"


# Results directory
resdir <- "results/"


# Processing directory
procdir <- "results/processing_normTest/"
dir.create(procdir)


# List of input files
rcountsfiles  <- list.files(datadir, full.names = T) %>% print()
bccountsfiles <- list.files(path = "results/processing", pattern = "BCcounts.RDS", 
                            recursive = T, full.names = T) %>% print()


# List of experiment names (proteins)
experiments <- rcountsfiles %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_ReadCounts.RDS)") %>%
        str_sort() %>% print()




# Replicates to be excluded from analysis --------------------------------------

badreps <- tibble(experiment = c("Eno2", "Fas1", "Fas2", "Htb2", "Rnr2", 
                                 "Rpl9A", "Ssa1", "Tdh1", "Tdh2", "Tdh3", "Yhb1"),
                  replicate = list("A", c("A", "F"), c("A", "D"), "A", NA, 
                                   NA, NA, NA, NA, NA, NA)) %>% print()




# Compare normalization by total BC counts vs read counts ----------------------

stats_perGuide  <- list()
stats_perSample <- list()

for(exp in seq_along(experiments)) {
        
        print(experiments[exp])
        #exp = 1
        
        
        # Get BC counts per guide
        BCcounts   <- readRDS(bccountsfiles[exp]) %>%
                select(-tail, -repl) #%>% print()
        
        
        # Get read counts per guide
        ReadCounts <- readRDS(rcountsfiles[exp]) %>% 
                as_tibble() %>%
                select(-barcode) %>%
                group_by(guide) %>%
                summarize(across(where(is.numeric), sum)) %>%
                pivot_longer(-guide, names_to = "sample", values_to = "nreads") %>%
                mutate(sample = str_replace(sample, "_", ".")) #%>% print()
        
        
        # Remove replicates specified above to be excluded
        excl <- badreps %>% filter(experiment == experiments[exp]) %>%
                pull("replicate") %>% .[[1]]
        if(!NA %in% excl) { ReadCounts <- ReadCounts[!(str_extract(ReadCounts$sample, "(?<=..\\.).") %in% excl), ] }
        if(!NA %in% excl) { BCcounts <- BCcounts[!(str_extract(BCcounts$sample, "(?<=..\\.).") %in% excl), ] }
        
        
        # Get BC counts vs read counts per guide
        comb <- BCcounts %>%
                left_join(ReadCounts, by = c("guide", "sample")) #%>% print()
        
        stats_perGuide[[exp]] <- comb
        
        
        # Get read counts per sample
        ReadCounts_tot <- ReadCounts %>%
                group_by(sample) %>%
                summarize(nreads = sum(nreads)) #%>% print()
        
        
        # Get BC counts per sample
        BCcounts_tot <- BCcounts %>%
                group_by(sample) %>%
                summarize(nBCs = sum(nBCs)) #%>% print()
        
        
        # Get BC counts vs read counts per sample
        comb_tot <- BCcounts_tot %>%
                left_join(ReadCounts_tot, by = "sample") #%>% print()
        
        stats_perSample[[exp]] <- comb_tot
}


names(stats_perGuide)  <- experiments
names(stats_perSample) <- experiments


# Plot read counts vs BC counts per sample
stats_perSample %>%
        bind_rows(.id = "protein") %>%
        ggplot(aes(x = nBCs, y = nreads)) +
                geom_point(alpha = 0.3) +
                ggpubr::stat_cor(color = "grey25",
                                 r.digits = 3, 
                                 label.x.npc = "left", 
                                 label.y.npc = "top", 
                                 label.sep = "\n") +
        scale_x_continuous(labels = function(x) format(x, scientific = T)) +
        facet_wrap(~ protein) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.margin = margin(5, 10, 5, 5))
ggsave(paste0(procdir, "nBCs-vs-nreads_perSample.jpg"), width = 4, height = 4)
ggsave(paste0(procdir, "nBCs-vs-nreads_perSample.pdf"), width = 4, height = 4)
ggsave(paste0(procdir, "nBCs-vs-nreads_perSample_facet.jpg"), width = 7, height = 6)
ggsave(paste0(procdir, "nBCs-vs-nreads_perSample_facet.pdf"), width = 7, height = 6)


# Plot read counts vs BC counts per guide
stats_perGuide %>%
        bind_rows(.id = "protein")  %>%
        ggplot(aes(x = nBCs, y = nreads)) +
        geom_hex(bins = 80) +  # geom_point takes very long to compute and plot
        scale_fill_gradient(low = "grey", high = "black", trans = "log10",
                            name = "Count per bin") +
        ggpubr::stat_cor(color = "grey25",
                         r.digits = 3, 
                         label.x.npc = "left", 
                         label.y.npc = "top", 
                         label.sep = "\n") +
        facet_wrap(~ protein) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
ggsave(paste0(procdir, "nBCs-vs-nreads_perGuide.jpg"), width = 4.5, height = 4)
ggsave(paste0(procdir, "nBCs-vs-nreads_perGuide.pdf"), width = 4.5, height = 4)
ggsave(paste0(procdir, "nBCs-vs-nreads_perGuide_facet.jpg"), width = 7, height = 6)
ggsave(paste0(procdir, "nBCs-vs-nreads_perGuide_facet.pdf"), width = 7, height = 6)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-03c_SessionInfo.txt")

