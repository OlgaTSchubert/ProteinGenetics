# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE22_Proteomics/")
getwd()


# Results directory
resdir <- "results/volcanoes/"
dir.create(resdir)


# Enrichments directory
enrdir <- "results/enrichments/"


# Import and reformat artMS results
d <- read_delim("results/artMS/artMS_results_OS.tsv", delim = "\t") %>%
        mutate(adj.p.i = ifelse(is.na(adj.p), 1E-10, adj.p.i)) %>%
        select(geneSys, gene, mutant, log2FC.i, adj.p.i) %>%
        mutate(gene = str_to_title(gene)) %>%
        print()


# Import annotations
annot <- read_tsv(paste0(enrdir, "selectedGOs.tsv")) %>% print()



# Data points to label in plots ------------------------------------------------

volcanoLabels <- function(annot, mut, input_up, input_dn) {
        
        enr_up  <- annot %>% filter(mutant == mut, dir == "up") %>% pull(GOterm)
        enr_dn  <- annot %>% filter(mutant == mut, dir == "dn") %>% pull(GOterm)
        
        labs_up <- tibble()
        labs_dn <- tibble()
        
        if(length(enr_up) != 0) {
                for (i in 1:length(enr_up)) {
                        labs_up <- rbind(labs_up, 
                                         tibble(geneSys = input_up %>% filter(GOterm == enr_up[i]) %>% pull(hits) %>% strsplit(",") %>% unlist(),
                                                GOID    = input_up %>% filter(GOterm == enr_up[i]) %>% pull(GOID),
                                                GOterm  = input_up %>% filter(GOterm == enr_up[i]) %>% pull(GOterm),
                                                nall    = input_up %>% filter(GOterm == enr_up[i]) %>% pull(nall),
                                                nrefs   = input_up %>% filter(GOterm == enr_up[i]) %>% pull(nrefs),
                                                nset    = input_up %>% filter(GOterm == enr_up[i]) %>% pull(nset),
                                                nhits   = input_up %>% filter(GOterm == enr_up[i]) %>% pull(nhits),
                                                fdr     = input_up %>% filter(GOterm == enr_up[i]) %>% pull(fdr)))
                }
        } else { print("No enriched functions among upregulated genes.") }
        
        if(length(enr_dn) != 0) {
                for (i in 1:length(enr_dn)) {
                        labs_dn <- rbind(labs_dn, 
                                         tibble(geneSys = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(hits) %>% strsplit(",") %>% unlist(),
                                                GOID    = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(GOID),
                                                GOterm  = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(GOterm),
                                                nall    = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(nall),
                                                nrefs   = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(nrefs),
                                                nset    = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(nset),
                                                nhits   = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(nhits),
                                                fdr     = input_dn %>% filter(GOterm == enr_dn[i]) %>% pull(fdr)))
                }
        } else { print("No enriched functions among downregulated genes.") }
        
        labs <- rbind(labs_up, labs_dn)
        return(labs)
}


POP1_labs <- volcanoLabels(annot = annot, mut = "POP1",
                           input_up = read_tsv(paste0(enrdir, "POP1_up_enr.tsv")),
                           input_dn = read_tsv(paste0(enrdir, "POP1_dn_enr.tsv"))) %>% print()
SIT4_labs <- volcanoLabels(annot = annot, mut = "SIT4",
                           input_up = read_tsv(paste0(enrdir, "SIT4_up_enr.tsv")),
                           input_dn = read_tsv(paste0(enrdir, "SIT4_dn_enr.tsv"))) %>% print()
SSY5_labs <- volcanoLabels(annot = annot, mut = "SSY5",
                           input_up = read_tsv(paste0(enrdir, "SSY5_up_enr.tsv")),
                           input_dn = read_tsv(paste0(enrdir, "SSY5_dn_enr.tsv"))) %>% print()




# Volcano plots ----------------------------------------------------------------

volcanoPlot <- function(dat, mut, nbg, labs, hits = c(), outname) {
        
        (p <- dat %>%
                        left_join(labs, by = "geneSys") %>%
                        filter(mutant == mut) %>%
                        mutate(lab = ifelse(is.na(GOID), NA, 
                                            paste0(GOterm, " (", GOID, ")\n", 
                                                   "(hits = ", nhits, 
                                                   ", expected = ", format(nrefs/nall*nset, digits = 1), 
                                                   ", FDR = ", format(fdr, digits = 2), ")\n")),
                               lab = factor(lab)) %>%
                        mutate(negLogQ = -log10(adj.p.i),
                               negLogQ = ifelse(negLogQ > 10, 10, negLogQ)) %>%
                        ggplot(aes(x = log2FC.i, y = negLogQ)) +
                        geom_point(data = . %>% filter(is.na(lab)), color = "gray80") +
                        geom_point(data = . %>% filter(!is.na(lab)), aes(color = lab)) +
                        ggrepel::geom_text_repel(data = . %>% filter(gene %in% hits),
                                                 aes(label = gene), color = "black",
                                                 size = 4, box.padding = 1, max.overlaps = Inf) +
                        ylab("-Log10(adj. p)") +
                        xlab("Log2FC") +
                        labs(color = "Selected enriched GO terms") +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank()))
        ggsave(plot = p, filename = outname, width = 6.5, height = 2.5)
        return(p)
}


(POP1_plot <- volcanoPlot(dat = d, mut = "POP1", nbg = 2819, labs = POP1_labs,
                          outname = paste0(resdir, "POP1_volcano.pdf")))
(SIT4_plot <- volcanoPlot(dat = d, mut = "SIT4", nbg = 2819, labs = SIT4_labs,
                          outname = paste0(resdir, "SIT4_volcano.pdf")))
(SSY5_plot <- volcanoPlot(dat = d, mut = "SSY5", nbg = 2819, labs = SSY5_labs,
                          hits = c("Yhb1", "Ssy1", "Ptr3", "Stp1", "Stp2"),
                          outname = paste0(resdir, "SSY5_volcano.pdf")))


# Combined plot for POP1, SAP155, SIT4
all <- egg::ggarrange(plots = list(POP1_plot, SIT4_plot))
ggsave(all, filename = paste0(resdir, "Combined_volcano.pdf"), width = 6.5, height = 6)




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE22-03_SessionInfo.txt")

