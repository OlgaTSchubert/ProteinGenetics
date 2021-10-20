# Libraries --------------------------------------------------------------------

library(openxlsx)
library(ShortRead)
library(data.table)
library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE21_FitnessScreen")
getwd()


# Data directory
#datadir <- "data/"
datadir <- "~/BigDataFiles/2021_ProteinGenetics/BE21/"


# Processed reads directory
readsdir <- paste0(datadir, "ReadsRDS/")
dir.create(readsdir)


# Sample annotation
annot <- read.xlsx("BE21_Samples.xlsx") %>% print()
annot <- split(annot, annot$Experiment) %>% print()


# Preceding sequences of guides and barcodes
pregd <- DNAString("GCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATC")
prebc <- DNAString("CAACCTTACCAGAGGGCGCCCCAGCTGGCAATTCCGG")




# Extracting guide and barcode sequences ---------------------------------------

# Delete existing files in the ReadsRDS directory before running this code!

error.rate <- 0.05
nbuffer    <- 1e7

for(exp in seq_along(annot)) {
        
        print(names(annot)[exp])
        
        samples      <- annot[[exp]]        # Need this for handling of file names
        samplesDF    <- DataFrame(samples)  # Data will be collected here ($results)
        sample.seqs  <- list()
        
        
        for(sam in 1:nrow(samples)){
                
                print(samples$Sample[sam])
                
                # Note that there might be multiple R1 and R2 lanes/files per sample
                fastq_filesR1s <- samples[sam, str_which(names(samples), "FileR1")]
                fastq_filesR2s <- samples[sam, str_which(names(samples), "FileR2")]
                
                seqs <- list()
                iter <- 1
                
                for(ln in seq_along(fastq_filesR1s)) {
                        
                        in.file1  <- paste0(datadir, fastq_filesR1s[ln])
                        in.file2  <- paste0(datadir, fastq_filesR2s[ln])
                        
                        fi1 <- FastqStreamer(in.file1, nbuffer, readerBlockSize = 1e9, verbose = F)
                        fi2 <- FastqStreamer(in.file2, nbuffer, readerBlockSize = 1e9, verbose = F)
                        
                        repeat {
                                rfq1 <- yield(fi1) 
                                rfq2 <- yield(fi2) 
                                
                                if(length(rfq1) == 0 ) { break }
                                
                                
                                # Read1 (guide)
                                cread <- sread(rfq1)
                                
                                findStagger <- which.isMatchingStartingAt(pregd, cread, 
                                                                          starting.at = 1:60, 
                                                                          max.mismatch = round(nchar(pregd) * error.rate), 
                                                                          follow.index = T, 
                                                                          with.indels = F)
                                
                                r1.trimmed <- narrow(cread, start = findStagger)
                                q1.trimmed <- narrow(quality(rfq1), start = findStagger)
                                findStagger[width(r1.trimmed) < (nchar(pregd) + 20)] <- NA
                                
                                r1.trimmed <- r1.trimmed[!is.na(findStagger)]
                                q1.trimmed <- q1.trimmed[!is.na(findStagger)]
                                guide      <- narrow(r1.trimmed, start = nchar(pregd) + 1,
                                                     end = nchar(pregd) + 20)
                                
                                rfq1   <- rfq1[!is.na(findStagger)]
                                rfq2   <- rfq2[!is.na(findStagger)]
                                
                                
                                # Read2 (barcode)
                                cread2 <- sread(rfq2)
                                
                                barcodeStartIndex <- which.isMatchingEndingAt(prebc, cread2, 
                                                                              ending.at = 1:50, 
                                                                              max.mismatch = round(nchar(prebc) * error.rate), 
                                                                              follow.index = T, 
                                                                              with.indels = F)
                                
                                clipped.barcodes <- (barcodeStartIndex + 20) > width(cread2)
                                barcodeStartIndex[clipped.barcodes == T] <- NA
                                
                                barcode    <- narrow(cread2, start = barcodeStartIndex + 1, end = barcodeStartIndex + 20)
                                barcode    <- barcode[!is.na(barcodeStartIndex)]
                                
                                r2.trimmed <-  cread2[!is.na(barcodeStartIndex)]
                                q2.trimmed <-  quality(rfq2)[!is.na(barcodeStartIndex)]
                                r1.trimmed <-  r1.trimmed[!is.na(barcodeStartIndex)]
                                q1.trimmed <-  q1.trimmed[!is.na(barcodeStartIndex)]
                                guide      <-  guide[!is.na(barcodeStartIndex)]
                                rfq2       <-  rfq2[!is.na(barcodeStartIndex)]
                                rfq1       <-  rfq1[!is.na(barcodeStartIndex)]
                                
                                print(iter)
                                iter <- iter + 1
                                
                                # Combine guide and barcode sequences into a table and append to list
                                seqs[[as.character(iter)]] <- DataFrame(guide = guide, barcode = barcode)
                        }
                        close(fi1)
                        close(fi2)
                }
                sample.seqs[[samplesDF$Sample[sam]]] <- do.call("rbind", seqs)
        }
        samplesDF$results <- sample.seqs
        saveRDS(samplesDF, file = paste0(readsdir, names(annot)[exp], "_Reads.RDS"))
}




# Quick check if read numbers are as expected ----------------------------------

samplesDF
samplesDF$results
sapply(samplesDF$results, nrow)
barplot(sapply(samplesDF$results, nrow))




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE21-01_SessionInfo.txt")

