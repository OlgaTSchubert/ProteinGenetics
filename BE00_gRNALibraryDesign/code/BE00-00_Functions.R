# Build whole genome DNA string set with reverse complement
buildSacCer3_genome_dict <- function(sacCer3, chrs) {
        
        sc3.dna <- list()
        
        for(chr in  chrs ) { 
                sc3.dna[[chr]]     <- DNAString(sacCer3[[chr]])
                rc.name            <- paste0(chr, "_rc")
                sc3.dna[[rc.name]] <- reverseComplement(DNAString(sacCer3[[chr]]))
        }
        
        sc3.set <- DNAStringSet(sc3.dna)
}




# Build pdict for matching guide seed sequences to the genome
# Includes GG of PAM sequence and 1-12 nts upstream of the PAM
build_4base_pdict <- function(guide.seq) { 
        
        guide.seeds.N13.pdict <- list()
        
        for(base in c("A", "C", "G", "T") ) {
                gtable.seq      <- sapply(guide.seq, s2c)[9:23,]
                gtable.seq[13,] <- base
                all.guide.seeds <- as.vector(apply(gtable.seq,2, c2s))
                seed.set        <- DNAStringSet(unlist(all.guide.seeds))
                guide.seeds.N13.pdict[[base]] <- PDict(seed.set)
        }
        
        return(guide.seeds.N13.pdict)
}




# Locate PAM sequences and predict base editor effect at every targetable PAM site
# If returnIdentical then only output gRNAs that target a region without Cs
get_gRNAs_and_predict_edit <- function(nchr, sacCer3, txdb, sc3.set , returnIdentical = F) {
        
        PAMp_string <- DNAString("GG")
        PAMm_string <- reverseComplement(PAMp_string)
        chr.name    <- names(sacCer3)[nchr]
        
        PAMp <- matchPattern(PAMp_string, sacCer3[[nchr]], fixed = T)
        PAMm <- matchPattern(PAMm_string, sacCer3[[nchr]], fixed = T)
        
        gPAM <- GRanges(seqnames <- names(sacCer3)[nchr], 
                        ranges   <- IRanges(start = c(start(ranges(PAMp))-1, start(ranges(PAMm))), 
                                            end = c(end(ranges(PAMp)), end(ranges(PAMm))+1)),
                        strand   <- c(rep("+", length(PAMp)), rep("-", length(PAMm))))
        
        genome(gPAM) <- "sacCer3"
        
        gPAM$guide23 <- rep("", length(gPAM))
        
        too.close.to.ends <- which( ((start(gPAM)-20)<0 ) | ((end(gPAM)+20) > length(sacCer3[[nchr]])) )
        
        gPAM              <- gPAM[-too.close.to.ends]
        
        # For pos strand: Extract targeting seq and PAM
        soi <- as.vector(strand(gPAM) == "+")
        gPAM$guide23[soi] <- getSeq(sacCer3,chr.name, 
                                    start = start(gPAM)[soi]-20,
                                    end = start(gPAM)[soi]+2,
                                    as.character = T)
        
        # For pos strand: Extract extended targeting seq and PAM (required for efficiency prediction)
        gPAM$guide30[soi] <- getSeq(sacCer3,chr.name, 
                                    start = start(gPAM)[soi]-20-4,
                                    end = start(gPAM)[soi]+2+3,
                                    as.character = T)
        
        # For neg strand: Extract targeting seq and PAM
        soi <- as.vector(strand(gPAM) == "-")
        gPAM$guide23[soi] <- as.character(reverseComplement(DNAStringSet(getSeq(sacCer3,chr.name, 
                                                                                start = end(gPAM)[soi]-2,
                                                                                end = end(gPAM)[soi]+20,
                                                                                as.character = T) )))
        
        # For neg strand: Extract extended targeting seq and PAM (required for efficiency prediction)
        gPAM$guide30[soi] <- as.character(reverseComplement(DNAStringSet(getSeq(sacCer3,chr.name, 
                                                                             start = end(gPAM)[soi]-2-3,
                                                                             end = end(gPAM)[soi]+20+4,
                                                                             as.character = T) )))
        gPAM$PAMstrand <- strand(gPAM)
        
        # Modify gPAM such that ranges match where the sequences start and stop
        start(ranges(gPAM[strand(gPAM) == "+"])) <- start(ranges(gPAM[strand(gPAM) == "+"]))-20
        end(ranges(gPAM[strand(gPAM) == "-"]))   <- end(ranges(gPAM[strand(gPAM) == "-"]))+20
        
        # Edit all Cs to Ts for positions 4-8 
        gPAMseq <- t(sapply(gPAM$guide23, s2c))
        for(p in 4:8) { gPAMseq[gPAMseq[,p] == "C", p] = "T" }
        
        # Turn back into strings
        editedSequences <- apply(gPAMseq, 1, c2s)
        
        # Reverse complement 
        editedSequences.stringset     <- DNAStringSet(editedSequences)
        editedSequences.stringset.fwd <- editedSequences.stringset
        editedSequences.stringset.fwd[as.vector(strand(gPAM) == "-")] <- reverseComplement(editedSequences.stringset.fwd[as.vector(strand(gPAM) == "-")])
        
        # Change to forward strand
        strand(gPAM)[strand(gPAM) == "-"] <- "+"
        
        # Add expected edited sequence
        gPAM$editedSequenceAsIs <- editedSequences.stringset
        gPAM$editedSequenceFwd  <- editedSequences.stringset.fwd
        
        # Remove guides that do not introduce any changes
        identical <- which(gPAM$guide23 == gPAM$editedSequenceAsIs)
        if(returnIdentical) {
                gPAM <- gPAM[identical]
        } else {
                gPAM <- gPAM[-identical]
        }
        
        pdict_for_guide_seqs <- build_4base_pdict(gPAM$guide23)
        
        matches.list <- lapply(pdict_for_guide_seqs, function(seed.dict) {
                vcountPDict(seed.dict, sc3.set, max.mismatch = 0, with.indels = F)  
        })

        gPAM$matches <- rowSums(sapply(matches.list, rowSums))
        
        return(gPAM)
}




# Get list of provean scores per gRNA
get_provean_scores <- function(coding.effects, provean.scores) {
        
        provean.effects <- list()
        
        pb <- txtProgressBar(min = 1, max = nrow(coding.effects), style = 3)
        
        for(i in 1:nrow(coding.effects))  {
                
                setTxtProgressBar(pb,i)
                
                if(coding.effects$CONSEQUENCE[i] != "synonymous") {
                        
                        goi   <- coding.effects$GENEID[i]
                        refAA <- as.matrix(coding.effects$REFAA[i])
                        varAA <- as.matrix(coding.effects$VARAA[i])
                        pdiff <- which(refAA != varAA)
                        loi   <- coding.effects$PROTEINLOC[[i]][1]+pdiff-1
                        varoi <- varAA[pdiff] 
                        lmat  <- cbind(loi, match(varoi, colnames(provean.scores[[goi]])))
                        
                        prot.seq.rownames <- rownames(provean.scores[[goi]])
                        
                        if(!is.null(prot.seq.rownames)){
                                provean.effects[[as.character(i)]]  <- c(provean.scores[[goi]][lmat] )
                        } else { provean.effects[[as.character(i)]] <- numeric(0) }
                        
                } else {
                        provean.effects[[as.character(i)]] <- numeric(0)
                }
        }
        close(pb)
        
        return(provean.effects)
}




get_max_absolute_provean <- function(coding.effects) {
        
        hasProveanScores <- which(sapply(coding.effects$provean, length) > 0)
        maxAbsProvean    <- sapply(coding.effects[hasProveanScores,]$provean, function(x) max(abs(x), na.rm = T))
        maxAbsProvean[!is.finite(maxAbsProvean)] <- NA
        
        coding.effects$maxAbsProvean <- NA
        coding.effects$maxAbsProvean[hasProveanScores] <- maxAbsProvean
        
        return(coding.effects)
}




guideAAmutations <- function(guide, geneSys, prot.start, refaa, varaa) {
        
        # Generate tibble with mutation information for each guide
        muts <- tibble(guide = guide, geneSys = geneSys, prot.start = prot.start, refaa = refaa, varaa = varaa) %>%
                
                # For each guide, find positions of AAs that change and store in list column
                mutate(localAAPos = map2(strsplit(refaa, ""), strsplit(varaa, ""), ~ which(.x != .y)),
                       protAAPos  = map2(localAAPos, prot.start, ~ .x + .y[[1]])) %>%
                
                # Make all elements of the list column the same length
                mutate(localAAPos = map(.$localAAPos, `length<-`, max(lengths(.$localAAPos))),
                       protAAPos  = map(.$protAAPos, `length<-`, max(lengths(.$protAAPos)))) %>%
                
                # Extract information on the first AA change
                mutate(aa1pos = map_int(protAAPos, ~ .x[[1]]),
                       aa1in  = map2_chr(strsplit(refaa, ""), localAAPos, ~ .x[.y[[1]]]),
                       aa1out = map2_chr(strsplit(varaa, ""), localAAPos, ~ .x[.y[[1]]]),
                       mut1   = str_c(aa1in, aa1pos-1, aa1out)) %>%
                
                # Same for the second AA change, if applicable
                mutate(aa2pos = map_int(protAAPos, ~ .x[[2]]),
                       aa2in  = map2_chr(strsplit(refaa, ""), localAAPos, ~ .x[.y[[2]]]),
                       aa2out = map2_chr(strsplit(varaa, ""), localAAPos, ~ .x[.y[[2]]]),
                       mut2   = str_c(aa2in, aa2pos-1, aa2out)) %>%
                
                # Same for the third AA change, if applicable
                mutate(aa3pos = map_int(protAAPos, ~ .x[[3]]),
                       aa3in  = map2_chr(strsplit(refaa, ""), localAAPos, ~ .x[.y[[3]]]),
                       aa3out = map2_chr(strsplit(varaa, ""), localAAPos, ~ .x[.y[[3]]]),
                       mut3   = str_c(aa3in, aa3pos-1, aa3out)) %>%
                
                # Keep only summary notation for the three positions
                select(c(guide, geneSys, prot.start, mut1, mut2, mut3)) #%>% print()
        
        
        # Return only summary notation for the three positions
        return(muts)
        
}

