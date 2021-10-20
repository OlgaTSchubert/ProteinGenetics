# Libraries --------------------------------------------------------------------

library(Biostrings)
library(tidyverse)



# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE00_gRNALibraryDesign/")
getwd()


# Results directory
resdir <- "results/"


# Guides table
guidesAll   <- readRDS(paste0(resdir, "guidesAll.RDS")) %>% print()
codingAll   <- guidesAll %>% filter(type == "coding") %>% print()
noCsAll     <- guidesAll %>% filter(type == "noCs") %>% print()
nonYeastAll <- guidesAll %>% filter(type == "nonYeast") %>% print()


# Primer binding sites for oligo subsets (by Sri Kosuri)
skpF2  <- "CGCGTCGAGTAGGGT"
skpF3  <- "CGATCGCCCTTGGTG"
skpF4  <- "GGTCGAGCCGGAACT"
skpF5  <- "TCCCGGCGTTGTCCT"
skpF6  <- "CGCAGGGTCCAGAGT"

skpR2r <- "GCCGTGTGAAGCTGG"
skpR3r <- "GGTTTAGCCGGCGTG"
skpR4r <- "GGATGCGCACCCAGA"
skpR5r <- "GCTCCGTCACTGCCC"
skpR6r <- "GTTCGCGCGAAGGAA"  # Forms strong primer dimer

skpR2  <- as.character(reverseComplement(DNAString(skpR2r)))
skpR3  <- as.character(reverseComplement(DNAString(skpR3r)))
skpR4  <- as.character(reverseComplement(DNAString(skpR4r)))
skpR5  <- as.character(reverseComplement(DNAString(skpR5r)))
skpR6  <- as.character(reverseComplement(DNAString(skpR6r)))


# Restriction site to remove primer binding sites after PCR
MluI <- "ACGCGT"


# Homologous sequences for Gibson assembly (30 nt)
gibsonL <- "AAACTTCTCCGCAGTGAAAGATAAATGATC"
gibsonR <- "GTTTTAGAGCTAGAAATAGCAAGTTAAAAT"




# Defining guide subsets -------------------------------------------------------

# Prepare the following subsets of guides:
#   5442 guides (all) for stops in essential genes
#   5500 guides for stops in non-essential genes
#   5500 guides for highly conserved non-syn. mutations in essential genes
#   5500 guides for highly conserved non-syn. mutations in non-essential genes
#    198 guides (all) for ADE2 and CAN1 genes 

# Add different types of control oligos to "stps in essential genes"-subset
#    500 guides not targeting yeast genome
#    500 guides not containing a C in target window
#    500 guides introducing synonymous mutations


# Number of guides per set
n.ctr    <- 500
n.guides <- 5500




# Control guides not targeting yeast genome ------------------------------------

set.seed(1)
nonYeast <- nonYeastAll %>%
        filter(U6 == F) %>%
        sample_n(n.ctr) %>% 
        mutate(set = "Ctr1") %>% print(n = 500)




# Control guides without Cs in editing window ----------------------------------

set.seed(2)
noCs <- noCsAll %>%
        filter(U6 == F,
               matches == 1) %>% #print() # total: 288,746
        sample_n(n.ctr) %>% 
        mutate(set = "Ctr2") %>% print(n = 500)




# Control guides introducing synonymous changes --------------------------------

set.seed(3)
synon <- codingAll %>%
        filter(consequence == "synonymous",
               U6 == F,
               matches == 1,
               overlap == F) %>% #print() # total: 80,084
        sample_n(n.ctr) %>% 
        mutate(set = "Ctr3") %>% print(n = 500)



# Guides introducing essential stops -------------------------------------------

eStops <- codingAll %>%
        filter(essential == T,
               consequence == "stop",
               U6 == F,
               matches == 1,
               overlap == F) %>% 
        mutate(set = "eStops") %>% print()




# Guides introducing non-essential stops ---------------------------------------

set.seed(5)
neStops <- codingAll %>%
        filter(essential == F,
               consequence == "stop",
               U6 == F,
               matches == 1,
               overlap == F) %>%  #print() # total: 19,772
        sample_n(n.guides) %>% 
        mutate(set = "neStops") %>% print()




# Guides mutating positions with high Provean in essential genes ---------------

set.seed(6)
eProvs <- codingAll %>%
        filter(essential == T,
               consequence == "nonsynonymous",
               maxAbsProvean > 5,
               U6 == F,
               matches == 1,
               overlap == F) %>% #print() # total: 8,955
        sample_n(n.guides) %>% 
        mutate(set = "eProvs") %>% print()




# Guides mutating positions with high Provean in non-essential genes -----------

set.seed(7)
neProvs <- codingAll %>%
        filter(essential == F,
               consequence == "nonsynonymous",
               maxAbsProvean > 5,
               U6 == F,
               matches == 1,
               overlap == F) %>% #print() # total: 33,357
        sample_n(n.guides) %>% 
        mutate(set = "neProvs") %>% print()




# All guides targeting in CAN1 and ADE2 ----------------------------------------
# Note that in the original code there was an error and 28 guides for CAN1 and
# ADE2 were selected that either contain a U6 terminator and/or have multiple 
# matches across the yeast genome.

Can1Ade2 <- codingAll %>%
        filter(geneSys %in% c("YEL063C", "YOR128C"),
               U6 == F,
               matches == 1,
               overlap == F) %>%
        mutate(set = "Can1Ade2") %>% print()




# Concatenate all selected guides ----------------------------------------------

guidesSel <- bind_rows(nonYeast,
                    noCs,
                    synon,
                    eStops,
                    neStops,
                    eProvs,
                    neProvs,
                    Can1Ade2) %>% print()



# Note that because set.seed is not working consistently across R versions, the
# guide selection above does not replicate the guide selection we did previously.
# The following code regenerates the original guide library that was synthesized.

guidesOrig <- read_delim("annotations/synthesizedGuides.tsv", delim = "\t") %>%
        left_join(guidesAll, by = "guide") %>% print()




# Compile oligos to be synthesized ---------------------------------------------

#   - Concatenate primer binding sites, restriction sites, gibson homology arms
#     and guide sequence
#   - Reverse complement if more A's than T's (this is the case for most oligos)
#   - Remove oligos that contain unwanted MluI restriction sites

#guides <- guidesSel %>%
guides <- guidesOrig %>%
        mutate(oligo = case_when(set == "Ctr1"     ~ paste0(skpF2, MluI, gibsonL, guide, gibsonR, MluI, skpR2),
                                 set == "Ctr2"     ~ paste0(skpF2, MluI, gibsonL, guide, gibsonR, MluI, skpR2),
                                 set == "Ctr3"     ~ paste0(skpF2, MluI, gibsonL, guide, gibsonR, MluI, skpR2),
                                 set == "eStops"   ~ paste0(skpF2, MluI, gibsonL, guide, gibsonR, MluI, skpR2),
                                 set == "neStops"  ~ paste0(skpF3, MluI, gibsonL, guide, gibsonR, MluI, skpR3),
                                 set == "eProvs"   ~ paste0(skpF4, MluI, gibsonL, guide, gibsonR, MluI, skpR4),
                                 set == "neProvs"  ~ paste0(skpF5, MluI, gibsonL, guide, gibsonR, MluI, skpR5),
                                 set == "Can1Ade2" ~ paste0(skpF6, MluI, gibsonL, guide, gibsonR, MluI, skpR6))) %>%
        mutate(As = str_count(oligo, "A"),
               Ts = str_count(oligo, "T"),
               oligo.rc = ifelse(As > Ts, as.character(reverseComplement(DNAStringSet(oligo))), oligo)) %>%
        mutate(MluIcount = str_count(oligo.rc, MluI)) %>%
        filter(MluIcount == 2) %>% print()


# Remove 28 guides erroneously added to Can1Ade2 set in original selection
guides <- guides %>% filter(U6 == F, matches < 2) %>% print()




# Save files -------------------------------------------------------------------

guides %>%
        select(-oligo, -As, -Ts, -oligo.rc, -MluIcount) %>%
        saveRDS(paste0(resdir, "guides.RDS"))

guides %>%
        select(-oligo, -As, -Ts, -oligo.rc, -MluIcount) %>%
        write_delim(paste0(resdir, "guides.tsv"), delim = "\t")

guides %>%
        select(oligo.rc) %>%
        write_delim(paste0(resdir, "oligos.tsv"), delim = "\t", col_names = F)




# Numbers for paper ------------------------------------------------------------

guides # 23,541 guides

table(guides$set)
# Can1Ade2  169
# Ctr1      492
# Ctr2      498
# Ctr3      500
# eProvs   5484
# eStops   5430
# neProvs  5483
# neStops  5485

guides %>% filter(set %in% c("neStops", "eProvs", "neProvs"))    # 16,452 guides

guides %>% filter(set %in% c("neStops", "eProvs", "neProvs")) %>%
        pull(geneSys) %>% unique() %>% length()                  #  4,592 genes

# All genes annotated in SGD: 6600
4592/6600  # 70% of annotated yeast genes targeted by our guides




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE00-02_SessionInfo.txt")

