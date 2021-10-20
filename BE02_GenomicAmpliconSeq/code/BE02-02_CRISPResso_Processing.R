# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE02_GenomicAmpliconSeq/")
getwd()


# Results directory
resdir <- "results/"
dir.create(resdir)


# CRISPResso2 results file
resfile <- "Quantification_window_nucleotide_frequency_table.txt"


# Files to be imported
files <- tibble(file = list.files(pattern = resfile, 
                                  recursive = T, full.names = T)) %>%
        pull(file) %>%
        print()


# gRNA plasmid names
guides <- str_extract(files, "(?<=\\.\\/).*(?=/CRISPResso_on_)") %>% print()


# DNA strand targeted by gRNA
strand <- tibble(guide = guides, strand = c("f", "r", "r", "r", "r", "f", "f", "r")) %>% print()

rev <- strand %>%
        filter(strand == "r") %>%
        pull(guide) %>%
        print()

fwd <- strand %>%
        filter(strand == "f") %>%
        pull(guide) %>%
        print()



# Processing data --------------------------------------------------------------

# Import all results tables and store in a list
comblist <- files %>%
        map(~ read_tsv(., col_names = F))
names(comblist) <- guides
comblist


# Get total number of reads per amplicon
totreads <- comblist %>%
        map(~ filter(., !is.na(X1)) %>%
                    pull(., X2) %>%
                    as.numeric(.) %>% 
                    sum(.)) %>% print()


# Get sequence of each amplicon in tibble
nts <- comblist %>%
        map(~ select(., -X1) %>%
                    slice(., 1) %>%
                    unlist(., use.names = FALSE)) %>% print()

seq <- nts %>%
        map(~ tibble(pos = 1:length(.), orig = .)) %>% print()


# Combine into table
comblist2 <- comblist %>%
        map2(totreads, 
             ~ filter(.x, X1 %in% c("A", "C", "G", "T")) %>%
                     pivot_longer(., -X1, names_to = "pos", values_to = "n") %>%
                     mutate(., pos = str_sub(pos, start = 2),
                               pos = as.numeric(pos)-1) %>%
                     mutate(., n   = as.numeric(n),
                               f   = n/.y*100) %>%
                     rename(., edit = X1)) %>% print()


# Add original base
comblist3 <- comblist2 %>%
        map2(seq, 
             ~ left_join(.x, .y, by = "pos") %>%
                     select(., pos, orig, edit, n, f)) %>% print()


# For guides targeting reverse strand, reverse position and rev-complement bases
comblist4 <- comblist3 %>%
        map_at(rev, ~ mutate(., pos = 41-pos) %>%
                       mutate(., orig = case_when(orig == "A" ~ "T",
                                                  orig == "C" ~ "G",
                                                  orig == "G" ~ "C",
                                                  orig == "T" ~ "A")) %>%
                       mutate(., edit = case_when(edit == "A" ~ "T",
                                                  edit == "C" ~ "G",
                                                  edit == "G" ~ "C",
                                                  edit == "T" ~ "A"))) %>% print()


# Add column for edit name
comblist5 <- comblist4 %>%
        map(~ mutate(., type = paste0(orig, "-to-", edit))) %>% print()


# Collapse the list into a single table
comb <- comblist5 %>%
        enframe(name = "guide", value = "value") %>%
        unnest(cols = c(value)) %>%
        filter(!str_detect(type, "NA-to")) %>% 
        mutate(pos = pos-30) %>%  # Change coordinates to PAM = 0
        print()




# Get editing efficiency per gRNA ----------------------------------------------

maxeff <- comb %>%
        filter(type == c("C-to-T")) %>%
        group_by(guide) %>%
        summarize(fmax = max(f)) %>% print()
write_csv(maxeff, paste0(resdir, "MaxEfficiency.csv"))

# guide fmax
# pLK78 48.4 
# pOS08 85.9 
# pOS09 67.9 
# pOS10 28.2 
# pOS11 36.6 
# pOS12 58.9 
# pOS13 83.2 
# pOS14  6.52



# Pie chart for fraction of gRNAs by editing efficiency ------------------------

pie <- comb %>%
        filter(type == c("C-to-T")) %>%
        group_by(guide) %>%
        summarize(fmax = max(f)) %>%
        mutate(eff = case_when(fmax > 75 ~ "75-100%",
                               fmax > 50 ~ "50-75%",
                               fmax > 25 ~ "25-50%",
                               TRUE ~ "0-25%")) %>%
        group_by(eff) %>%
        summarize(n = n()) %>% 
        arrange(desc(eff)) %>%
        mutate(prop = n / sum(n) * 100) %>% print()
write_csv(pie, paste0(resdir, "Piechart.csv"))

pie %>% ggplot(aes(x = "", y = prop, fill = eff)) +
        geom_bar(stat = "identity", color = "black", size = 0.3) +
        coord_polar("y", start = 0) +
        scale_fill_grey(start = 0.99, end = 0.6) +
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill = "Editing\nefficiency") +
        theme_void()
ggsave(paste0(resdir, "Piechart.pdf"), width = 2.5, height = 2)




# Plot conversion rate by position ---------------------------------------------

# Overview
comb %>%
        ggplot() +
        geom_rect(aes(xmin = -20.5, xmax = -0.5, ymin = 0, ymax = 100),
                  fill = "gray90") +
        geom_rect(aes(xmin = -12.5, xmax = -17.5, ymin = 0, ymax = 100),
                  fill = "gray85") +
        geom_point(aes(x = pos, y = f, color = guide)) + 
        xlim(-32, 12) +
        ylim(0, 100) +
        facet_wrap(~type) +
        labs(x = "Position relative to PAM", 
             y = "Conversion rate (%)",
             color = "gRNA") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "EditingWindow_all.pdf"), width = 8, height = 8)


# Selected only (color by type)
comb %>%
        filter(orig != edit) %>%
        filter(pos > -25 & pos < 5) %>%
        ggplot() +
        geom_rect(aes(xmin = -20.5, xmax = -0.5, ymin = 0, ymax = 100),
                  fill = "gray90") +
        geom_rect(aes(xmin = -12.5, xmax = -17.5, ymin = 0, ymax = 100),
                  fill = "gray85") +
        annotate(geom = "text", label = "PAM", x = 0, y = 5, hjust = 0,
                 color = "gray30", size = 3) +
        geom_point(aes(x = pos, y = f, color = type), 
                   position = position_jitter(width = 0.1)) + 
        ylim(0, 100) +
        labs(x = "Position relative to PAM", 
             y = "Conversion rate (%)") +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
ggsave(paste0(resdir, "EditingWindow.pdf"), width = 4, height = 3.1)
ggsave(paste0(resdir, "EditingWindow_wider.pdf"), width = 5, height = 3.1)




# Write session info -----------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE02-02_SessionInfo.txt")

