# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE02_GenomicAmpliconSeq/")
getwd()


# Results directory
resdir <- "results/"
dir.create(resdir)


# Tested gRNAs
gds <- read_tsv("annotations/guides.tsv") %>% print()


# CRISPResso2 results file
resfile <- "Alleles_frequency_table_around_sgRNA_.*txt"


# Files to be imported
files <- tibble(file = list.files(pattern = resfile, 
                                  recursive = T, full.names = T)) %>%
    pull(file) %>%
    print()


# gRNA plasmid names
guides <- str_extract(files, "(?<=\\.\\/).*(?=/CRISPResso_on_)") %>% print()




# Processing data --------------------------------------------------------------

# Import all results tables and store in a list
comblist <- files %>%
    map(~ read_tsv(., col_names = T))
names(comblist) <- guides
comblist


# Collapse the list into a single table
comb <- comblist %>%
    enframe(name = "guide", value = "value") %>%
    unnest(cols = c(value)) %>%
    filter(`%Reads` > 1) %>%
    print()
write_csv(comb, paste0(resdir, "Alleles_frequency_table_around_sgRNA_combined.csv"))


# Manually add Classification column and populate with either of the following:
#   Wildtype
#   Full window only
#   Full window and outside
#   Partial window only
#   Partial window and outside
#   Outside only

comb2 <- read_csv(paste0(resdir, "Alleles_frequency_table_around_sgRNA_combined_manualClassification.csv")) %>% 
    select(guide, `%Reads`, Classification) %>%
    print()


# Find % of reads missing and classify as "Other"
missing_reads <- comb2 %>%
    group_by(guide) %>%
    summarize(`%Reads` = 100-sum(`%Reads`)) %>%
    mutate(Classification = "Other") %>%
    print()


# Final table with classification and labels for plotting
comb3 <- bind_rows(comb2, missing_reads) %>% 
    left_join(gds, by = "guide") %>%
    mutate(Classification = factor(Classification, levels = c("Wildtype",
                                                              "Other",
                                                              "Outside only",
                                                              "Partial window and outside",
                                                              "Partial window only", 
                                                              "Full window and outside", 
                                                              "Full window only"))) %>%
    mutate(guide_label = paste0(sequence, " (", guide, ")"),
           guide_label = factor(guide_label, levels = rev(c("ACGTCCAAAATTGAATGACT (pOS09)", 
                                                            "GATACGTTCTCTATGGAGGA (pLK78)", 
                                                            "AGTTACCCAAAGTGTTCCTG (pOS13)", 
                                                            "AAACCAATACATGTAACCAT (pOS08)",
                                                            "CTCCAATAACGGAATCCAAC (pOS10)",
                                                            "TCCTGCCCAGGCCGCTGAGC (pOS12)", 
                                                            "GCCCATTTTTCGGCGTACAA (pOS14)", 
                                                            "GAACCAGAACTCTGACAGTT (pOS11)")))) %>% print()



# Stacked bar chart ------------------------------------------------------------

comb3 %>% 
    filter(Classification != "Other") %>%  # Note that if this filter is active, only alleles with >1% frequency are considered
    ggplot(aes(x = guide_label, y = `%Reads`, fill = Classification)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(type = "seq", palette = 7, direction = 1) +
    labs(y = "Fraction of reads") +
    coord_flip() +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE))
ggsave(paste0(resdir, "StackedBar.pdf"), width = 8, height = 3)




# Write session info -----------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE02-02b_SessionInfo.txt")

