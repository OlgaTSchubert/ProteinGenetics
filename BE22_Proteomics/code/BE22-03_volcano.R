# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Functions
source("../BE19_GFPScreens/code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/volcano/"
dir.create(resdir)


# STRING-enrichments directory
stringdir <- "results/STRING-enrichments/"


# Import artMS results
d0  <- read.delim2("results/artMS/results_adjpvalue/results-log2fc-long.txt") %>% 
        as_tibble() %>% print()


# Import annotations
annot <- read_tsv("annotations/volcanoLabels.tsv") %>% print()


# Import gene annotations
sgd <- importExtData(dataset = "SGD_features", localfile = T) %>%
        dplyr::rename(geneSys = Feature_name, 
                      gene = Standard_gene_name) %>% 
        mutate(gene = ifelse(is.na(gene), geneSys, gene)) %>% 
        select(geneSys, gene) %>% print()




# Reformat data ----------------------------------------------------------------

# With imputed missing values
d <- d0 %>%
        mutate(mutant = str_sub(Comparison, 1, -4)) %>%
        mutate(iPvalue = ifelse(imputed == "yes", 1E-10, iPvalue)) %>%
        left_join(sgd, by = c("Protein" = "geneSys")) %>%
        select(Protein, gene, mutant, iLog2FC, iPvalue) %>%
        rename(geneSys = Protein,
               log2fc  = iLog2FC,
               q       = iPvalue) %>%
        mutate(gene    = str_to_title(gene),
               log2fc  = as.numeric(log2fc),
               q       = as.numeric(q),
               sig     = (abs(log2fc) > 1 & q < 0.05)) %>%
        print()

# Without imputed missing values - if using these, also do enrichments on non-imputed data
# d <- d0 %>%
#         filter(imputed == "no") %>%
#         mutate(mutant = str_sub(Comparison, 1, -4)) %>%
#         left_join(sgd, by = c("Protein" = "geneSys")) %>%
#         select(Protein, gene, mutant, log2FC, adj.pvalue) %>%
#         rename(geneSys = Protein,
#                log2fc  = log2FC,
#                q       = adj.pvalue) %>%
#         mutate(log2fc  = as.numeric(log2fc),
#                q       = as.numeric(q),
#                sig     = (abs(log2fc) > 1 & q < 0.05)) %>%
#         print()




# Data points to label in plots ------------------------------------------------

volcanoLabels <- function(annot, mut, input_up, input_dn) {
        
        enr_up  <- annot %>% filter(mutation == mut, direction == "up") %>% pull(term)
        enr_dn  <- annot %>% filter(mutation == mut, direction == "dn") %>% pull(term)
        
        prio_up <- annot %>% filter(mutation == mut, direction == "up") %>% select(term, priority)
        prio_dn <- annot %>% filter(mutation == mut, direction == "dn") %>% select(term, priority)
        
        labs_up <- tibble()
        labs_dn <- tibble()
        
        if(length(enr_up) != 0) {
                for (i in 1:length(enr_up)) {
                        labs_up <- rbind(labs_up, 
                                         tibble(geneSys = input_up %>% filter(term == enr_up[i]) %>% pull(inputGenes) %>% strsplit(",") %>% unlist(),
                                                term    = input_up %>% filter(term == enr_up[i]) %>% pull(term),
                                                descr   = input_up %>% filter(term == enr_up[i]) %>% pull(description),
                                                fdr     = input_up %>% filter(term == enr_up[i]) %>% pull(fdr)))
                }
                labs_up <- left_join(labs_up, prio_up, by = "term")
        } else { print("No enriched functions among upregulated genes.") }
        
        if(length(enr_dn) != 0) {
                for (i in 1:length(enr_dn)) {
                        labs_dn <- rbind(labs_dn, 
                                         tibble(geneSys = input_dn %>% filter(term == enr_dn[i]) %>% pull(inputGenes) %>% strsplit(",") %>% unlist(),
                                                term    = input_dn %>% filter(term == enr_dn[i]) %>% pull(term),
                                                descr   = input_dn %>% filter(term == enr_dn[i]) %>% pull(description),
                                                fdr     = input_dn %>% filter(term == enr_dn[i]) %>% pull(fdr)))
                }
                labs_dn <- left_join(labs_dn, prio_dn, by = "term")
        } else { print("No enriched functions among downregulated genes.") }
        
        labs <- rbind(labs_up, labs_dn)
        return(labs)
}


POP1_labs   <- volcanoLabels(annot = annot, mut = "POP1",   
                             input_up = read_tsv(paste0(stringdir, "enriched_POP1_up.tsv")),   
                             input_dn = read_tsv(paste0(stringdir, "enriched_POP1_dn.tsv")))
SAP155_labs <- volcanoLabels(annot = annot, mut = "SAP155",
                             input_up = read_tsv(paste0(stringdir, "enriched_SAP155_up.tsv")),
                             input_dn = read_tsv(paste0(stringdir, "enriched_SAP155_dn.tsv")))
SIT4_labs   <- volcanoLabels(annot = annot, mut = "SIT4",
                             input_up = read_tsv(paste0(stringdir, "enriched_SIT4_up.tsv")),
                             input_dn = read_tsv(paste0(stringdir, "enriched_SIT4_dn.tsv")))
SSY5_labs   <- volcanoLabels(annot = annot, mut = "SSY5",
                             input_up = read_tsv(paste0(stringdir, "enriched_SSY5_up.tsv")),
                             input_dn = read_tsv(paste0(stringdir, "enriched_SSY5_dn.tsv")))




# Volcano plots ----------------------------------------------------------------

volcanoPlot <- function(dat, mut, labs, hits = c(), outname) {
        
        # Remove all but one annotations per gene
        labs <- labs %>%
                arrange(priority) %>% 
                distinct(geneSys, .keep_all = T)
        
        (p <- dat %>%
                        left_join(labs, by = "geneSys") %>%
                        filter(mutant == mut) %>%
                        mutate(lab = ifelse(is.na(term), NA, paste0(descr, " (FDR = ", fdr, ")")),
                               lab = factor(lab),
                               lab = fct_reorder(lab, priority)) %>%
                        mutate(negLogQ = -log10(q),
                               negLogQ = ifelse(negLogQ > 10, 10, negLogQ)) %>%
                        arrange(-priority) %>% 
                        ggplot(aes(x = log2fc, y = negLogQ)) +
                        geom_point(data = . %>% filter(is.na(lab)), color = "gray80") +
                        geom_point(data = . %>% filter(!is.na(lab)), aes(color = lab)) +
                        ggrepel::geom_text_repel(data = . %>% filter(gene %in% hits),
                                                 aes(label = gene), color = "black",
                                                 size = 4, box.padding = 1, max.overlaps = Inf) +
                        ylab("-Log10(adj. p)") +
                        xlab("Log2FC") +
                        labs(color = "Enriched functional categories") +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank()))
        ggsave(plot = p, filename = outname, width = 7, height = 2.5)
        return(p)
}


(POP1_plot   <- volcanoPlot(dat = d, mut = "POP1",   labs = POP1_labs,
                            outname = paste0(resdir, "POP1_volcano.pdf")))
(SAP155_plot <- volcanoPlot(dat = d, mut = "SAP155", labs = SAP155_labs,
                            outname = paste0(resdir, "SAP155_volcano.pdf")))
(SIT4_plot   <- volcanoPlot(dat = d, mut = "SIT4",   labs = SIT4_labs,
                            outname = paste0(resdir, "SIT4_volcano.pdf")))
(SSY5_plot   <- volcanoPlot(dat = d, mut = "SSY5",   labs = SSY5_labs,
                            hits = c("Yhb1", "Ssy1", "Ptr3", "Stp1", "Stp2"),
                            outname = paste0(resdir, "SSY5_volcano.pdf")))


# Combined plot for POP1, SAP155, SIT4
egg::ggarrange(plots = list(POP1_plot, SAP155_plot, SIT4_plot),
          ncol = 1, heights = c(1, 1, 1),
          labels = c("POP1 H642Y", "SAP155 W702*", "SIT4 Q184*"))
ggsave(paste0(resdir, "Combined_volcano.pdf"), width = 8, height = 10)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-03_SessionInfo.txt")

