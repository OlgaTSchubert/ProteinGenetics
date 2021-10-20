# Libraries --------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(stringdist)
library(Biostrings)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Reads directory
#readsdir <- "data/ReadsRDS/"
readsdir <- "~/BigDataFiles/2021_ProteinGenetics/BE19/ReadsRDS/"


# Guide table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
        filter(set %in% c("eProvs", "neProvs", "neStops")) %>%
        print()


# List of input files (read counts per barcode)
readsfiles <- list.files(readsdir, pattern = ".*_Reads.RDS", full.names = T); readsfiles


# List of experiment names (proteins)
experiments <- readsfiles %>%
        str_extract(pattern = "(?<=/)[^/]*(?=_Reads.RDS)") %>%
        str_sort() %>%
        print()




# Barcode correction -----------------------------------------------------------

edit.dist <- 5  # Barcodes are collapsed if edit distance less than this

for(exp in seq_along(experiments)) {

        print(experiments[exp])
        pb <- txtProgressBar(min = 1, max = nrow(guides), style = 3)
        
        samples     <- readRDS(paste0(readsdir, experiments[exp], "_Reads.RDS"))
        sample.seqs <- samples$results
        

        # For every sample, make lists of guides with their barcodes (list of list)
        barcodes <- lapply(sample.seqs, function(x) {
                x$gmatch <- match(x$guide, guides$guide)
                x        <- x[order(x$gmatch, decreasing = F), ]
                x        <- x[!is.na(x$gmatch), ]
                return(split(as.vector(x$barcode), as.vector(x$guide)))
        })

        
        samps          <- names(sample.seqs)
        bcs.corr.list  <- list()
        
        for(gd in 1:nrow(guides)){
                
                setTxtProgressBar(pb, gd)
                
                guide <- guides$guide[gd]
                
                # For the selected guide, get all bc reads from all samples
                bcs.list <- lapply(barcodes, function(x) x[[guide]])
                bcs      <- unlist(bcs.list)
                
                # Proceed only if at least 3 reads per guide
                if(length(bcs) < 3) {next;}
                
                # Get RLE for barcodes
                bcs.order <- order(as.vector(bcs))
                bcs.rle   <- rle(bcs[bcs.order])
                
                # Proceed only if at least 2 barcodes per guide
                if(length(bcs.rle$values) == 1) {next;}
                
                # Combine into dataframe
                bcs.df <- data.frame(barcode = bcs.rle[[2]], 
                                     barcode.count = bcs.rle[[1]], 
                                     bgroup = "", stringsAsFactors = F)
                bcs.df <- bcs.df[order(bcs.df$barcode.count, decreasing = T), ]
                
                # Get edit distance and collapse barcodes if it is less than 
                # edit.dist to most abundant barcode
                edist <- stringdist(bcs.df$barcode[1], bcs.df$barcode)
                bcs.df$bgroup[edist < edit.dist] <- bcs.df$barcode[1]
                
                for(i in 2:nrow(bcs.df)) {
                        
                        nbc            <- bcs.df$bgroup
                        nbc[nbc == ""] <- bcs.df$barcode[nbc == ""]
                        
                        if(bcs.df$bgroup[i] != "") {next;}
                        
                        edist <- stringdist(bcs.df$barcode[i], nbc)
                        bcs.df$bgroup[edist < edit.dist] <- bcs.df$barcode[i]
                }
                
                # Look up barcode group for each read and collapse into table
                bcs.corr <- rbindlist(lapply(bcs.list, function(x) {
                                m <- match(x, bcs.df$barcode)
                                data.frame(barcode = bcs.df$bgroup[m]) }), 
                        idcol = "sample")
                
                # Account for cases where no barcodes at all in one or more samples
                bcs.corr$sample <- as.factor(bcs.corr$sample)
                missing.levels  <- samps[!(samps %in% bcs.corr$sample)]
                if(length(missing.levels) > 0) {
                        levels(bcs.corr$sample) <- c(levels(bcs.corr$sample), 
                                                     missing.levels)
                }
                
                # Create read counts table and add to results list
                bcs.corr.w <- bcs.corr %>% 
                        group_by(sample, barcode) %>% 
                        tally() %>% 
                        spread(., sample, n, 0, drop = F)
                bcs.corr.w <- bcs.corr.w[, match(c("barcode", samps), 
                                                 colnames(bcs.corr.w))]
                bcs.corr.list[[guide]] <- bcs.corr.w
        }
        
        close(pb)
        
        # Convert results list into table
        read.counts <- rbindlist(bcs.corr.list, idcol = "guide")
        
        # Fix sample names
        rep <- samples$Replicate[match(colnames(read.counts)[-c(1,2)], samples$SeqID)]
        grp <- samples$Group[match(colnames(read.counts)[-c(1,2)], samples$SeqID)]
        colnames(read.counts)[-c(1,2)] <- paste0(grp, "_", rep)
        
        saveRDS(read.counts, paste0(readsdir,  experiments[exp], "_ReadCounts.RDS"))
        
        # To avoid memory issues, delete some of the large files
        rm(read.counts, samples, sample.seqs)
}




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-02_SessionInfo.txt")

