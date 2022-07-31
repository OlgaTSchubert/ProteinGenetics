# Libraries --------------------------------------------------------------------

library(tidyverse)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE19_GFPScreens/")
getwd()


# Functions
source("code/BE19-00_ImportExtData.R")


# Results directory
resdir <- "results/litComparisons/"
dir.create(resdir)


# Import results
combdf.gn <- readRDS("results/processing/combdf_gn.RDS")

schubert <- combdf.gn %>%
        filter(!is.na(gene)) %>%   # remove genes YGR109W-B and YIL082W-A
        select(-geneSys, -geneDescr, -essential, -FDR0.05_count) %>%
        pivot_longer(-gene,
                     names_to = c("protein", ".value"),
                     names_sep = "_") %>%
        mutate(protein = toupper(protein)) %>%
        print()


# Get list of GFP proteins
prots    <- toupper(str_extract(str_subset(names(combdf.gn), "_log2fc"), ".*(?=_log2fc)")) %>% print()
protsSys <- importExtData(dataset = "SGD_features", localfile = T) %>%
        filter(Standard_gene_name %in% prots) %>%
        select("Feature_name") %>%
        pull() %>% print()




# Import external data ---------------------------------------------------------

kemmeren.o <- importExtData(dataset = "Kemmeren2014", localfile = T) %>% print()
oduibhir.o <- importExtData(dataset = "ODuibhir2014", localfile = T) %>% print()
reimand.o  <- importExtData(dataset = "Reimand2010",  localfile = T) %>% print()
stefely.o  <- importExtData(dataset = "Stefely2016",  localfile = T) %>% print()




# Process Kemmeren data --------------------------------------------------------

kemmeren.o <- kemmeren.o %>%
        select(-c(1, 2, tail(names(.), 9))) %>% # remove control samples at the end
        select(-str_which(names(.), "_1")) %>%  # remove rows containing "A value"
        filter(geneSymbol %in% prots) %>%
        rename(protein = geneSymbol) %>%
        mutate_at(vars(-protein), as.numeric) %>%
        print()

kemmeren.r <- kemmeren.o %>%
        select(-c(str_which(names(.), "_2"))) %>%
        rename_at(vars(-1), list(~ toupper(str_extract(., ".*(?=-del)")))) %>%
        pivot_longer(cols = -protein, names_to = "gene", values_to = "log2fc") %>% 
        print()

kemmeren.p <- kemmeren.o %>%
        select(c(1, str_which(names(.), "_2"))) %>% 
        rename_at(vars(-1), list(~ toupper(str_extract(., ".*(?=-del)")))) %>%
        pivot_longer(cols = -protein, names_to = "gene", values_to = "q") %>% 
        print()

kemmeren.comb <- kemmeren.r %>%
        left_join(kemmeren.p) %>%
        print()


# Fix mediator subunit names
kemmeren.comb <- kemmeren.comb %>%
        mutate(gene = ifelse(gene == "CYCC",  "SSN8",  gene),
               gene = ifelse(gene == "MED13", "SSN2",  gene),
               gene = ifelse(gene == "CDK8",  "SSN3",  gene),
               gene = ifelse(gene == "MED12", "SRB8",  gene),
               gene = ifelse(gene == "MED1",  "MED1",  gene),
               gene = ifelse(gene == "MED9",  "CSE2",  gene),
               gene = ifelse(gene == "MED31", "SOH1",  gene),
               gene = ifelse(gene == "MED2",  "MED2",  gene),
               gene = ifelse(gene == "MED3",  "PGD1",  gene),
               gene = ifelse(gene == "MED5",  "NUT1",  gene),
               gene = ifelse(gene == "MED15", "GAL11", gene),
               gene = ifelse(gene == "MED16", "SIN4",  gene),
               gene = ifelse(gene == "MED18", "SRB5",  gene),
               gene = ifelse(gene == "MED20", "SRB2",  gene)) %>%
        print()




# Process O'Duibhir data -------------------------------------------------------

oduibhir.r <- oduibhir.o %>%
        select(-1) %>%
        filter(commonName %in% prots) %>%
        rename(protein = commonName,           # Need to correct a few of the headers
               `mf(alpha)1.del.vs..wt` = mf.alpha.1.del.vs..wt,
               `mf(alpha)2.del.vs..wt` = mf.alpha.2.del.vs..wt,
               `arg5,6.del.vs..wt`     = arg5.6.del.vs..wt,
               `ymr031c.del.vs..wt`    = ymr031c.vs..wt,
               `yil014c-a.del.vs..wt`  = yil014c.a.del.vs..wt,
               `yol086w-a.del.vs..wt`  = yol086w.a.del.vs..wt,
               `ydr034w-b.del.vs..wt`  = ydr034w.b.del.vs..wt,
               `yal044w-a.del.vs..wt`  = yal044w.a.del.vs..wt) %>%
        rename_at(vars(-1), list(~toupper(str_extract(., ".*?(?=\\.del)")))) %>%
        pivot_longer(cols = -protein, names_to = "gene", values_to = "log2fc") %>%
        print()

oduibhir.comb <- oduibhir.r %>%
        left_join(kemmeren.p) %>%
        print()




# Process Reimand data ---------------------------------------------------------

reimand.r <- reimand.o[["ratios"]] %>%
        mutate(gene = toupper(Name)) %>%
        select(gene, protsSys) %>%
        rename_all(list(~ c("gene", prots))) %>%
        pivot_longer(-gene, names_to = "protein", values_to = "log2fc") %>%
        print()

reimand.p <- reimand.o[["pvals"]] %>%
        mutate(gene = toupper(Name)) %>%
        select(gene, protsSys) %>%
        rename_all(list(~ c("gene", prots))) %>%
        pivot_longer(-gene, names_to = "protein", values_to = "q") %>%
        print()

reimand.comb <- reimand.r %>%
        left_join(reimand.p) %>%
        print()



# Process Stefely data ---------------------------------------------------------

stefely.r <- stefely.o[["ratios"]] %>% 
        filter(`Molecule Type` == "Protein") %>%
        select(-`Molecule Type`) %>%
        mutate(protein = str_extract(`Molecule Name`, ".*(?= \\()")) %>%
        select(protein, everything(), -`Molecule Name`) %>%
        filter(protein %in% prots) %>%       # HTB2 not in this dataset!
        pivot_longer(-protein, names_to = "gene", values_to = "log2fc") %>%
        mutate(gene = str_sub(gene, start = 2L, end = -1L)) %>%
        print()

stefely.p <- stefely.o[["pvals"]] %>%
        filter(`Molecule Type` == "Protein") %>%
        select(-`Molecule Type`) %>%
        mutate(protein = str_extract(`Molecule Name`, ".*(?= \\()")) %>%
        select(protein, everything(), -`Molecule Name`) %>%
        filter(protein %in% prots) %>%       # HTB2 not in this dataset!
        pivot_longer(-protein, names_to = "gene", values_to = "q") %>%
        mutate(gene = str_sub(gene, start = 2L, end = -1L)) %>%
        print()

stefely.comb <- stefely.r %>%
        left_join(stefely.p) %>%
        print()




# Combine all datasets ---------------------------------------------------------

all.comb <- schubert %>%
        left_join(kemmeren.comb, by = c("protein", "gene"), suffix = c(".schubert", ".kemmeren")) %>%
        #left_join(oduibhir.comb, by = c("protein", "gene")) %>%
        #rename(log2fc.oduibhir = log2fc, q.oduibhir = q) %>%
        #left_join(reimand.comb, by = c("protein", "gene")) %>%
        #rename(log2fc.reimand = log2fc, q.reimand = q) %>%
        #left_join(stefely.comb, by = c("protein", "gene")) %>%
        #rename(log2fc.stefely = log2fc, q.stefely = q) %>%
        print()




# Plotting ---------------------------------------------------------------------

litCompScatter <- function(data, r1, r2, p1, p2, name1, name2, 
                           p = 0.05, r = 0.5, xli = 4, yli = 4, 
                           labels = c("Consistent", "Opposite"), 
                           genestocolor) {
        
        # To not add labels, use "labels = NULL"
        
        dat <- data %>%
                select(gene, protein, r1, r2, p1, p2) %>%
                rename_at(vars(contains(".")), list(~ c("r1", "r2", "p1", "p2"))) %>% 
                filter(!(is.na(r1) | is.na(r2))) %>%
                mutate(Class.pr = as.factor(case_when((p1 < 0.05 & abs(r1) > 0.5) & !(p2 < 0.05 & abs(r2) > 0.5) ~ "Study 1 only",
                                                      !(p1 < 0.05 & abs(r1) > 0.5) & (p2 < 0.05 & abs(r2) > 0.5) ~ "Study 2 only",
                                                      ((p1 < 0.05 & r1 > 0.5) & (p2 < 0.05 & r2 > 0.5)) |
                                                              ((p1 < 0.05 & r1 < -0.5) & (p2 < 0.05 & r2 < -0.5)) ~ "Consistent",
                                                      ((p1 < 0.05 & r1 > 0.5) & (p2 < 0.05 & r2 < -0.5)) |
                                                              ((p1 < 0.05 & r1 < -0.5) & (p2 < 0.05 & r2 > 0.5)) ~ "Opposite",
                                                      TRUE ~ "No effect"))) %>%
                mutate(Class.pr = factor(Class.pr, levels = c("No effect", 
                                                              "Study 1 only", 
                                                              "Study 2 only", 
                                                              "Opposite", 
                                                              "Consistent")))
        write_csv(dat, paste0(resdir, name1, "-", name2, ".csv")) # Factor levels lost if pipe into write_csv

        
        dat %>%
                ggplot() +
                geom_vline(xintercept = 0, color = "grey") +
                geom_hline(yintercept = 0, color = "grey") +
                geom_point(data = . %>% filter(Class.pr %in% c("No effect", "Study 1 only", "Study 2 only")),
                           aes(x = r1, y = r2), color = "grey80") +
                geom_point(data = . %>% filter(Class.pr == "Consistent"),
                           aes(x = r1, y = r2), color = "royalblue3") +
                geom_point(data = . %>% filter(Class.pr == "Opposite"),
                           aes(x = r1, y = r2), color = "cornflowerblue") +
                {if (!is.null(labels)) geom_point(data = . %>% filter(gene %in% genestocolor),
                           aes(x = r1, y = r2), color = "gold")} +
                geom_smooth(data = . %>% filter(Class.pr %in% c("Consistent", "Opposite")),
                            aes(x = r1, y = r2), color = "grey25", 
                            method = "lm") +
                ggpubr::stat_cor(data = . %>% filter(Class.pr %in% c("Consistent", "Opposite")),
                                 aes(x = r1, y = r2), color = "grey25",
                                 label.x = -xli+0.5, label.y = yli-0.5) +
                #label.x.npc = "left", label.y.npc = "top") +
                {if (!is.null(labels)) ggrepel::geom_text_repel(
                        data = . %>% filter(Class.pr %in% labels),
                        mapping = aes(x = r1, y = r2, label = paste0(gene, "\n(", protein, ")")),
                        size = 2, box.padding = 0.2, point.padding = 0.05,
                        segment.color = "darkgrey", segment.size = 0.5)} +
                labs(x = paste0(name1, " et al."), y = paste0(name2, " et al.")) +
                coord_cartesian(xlim = c(-xli, xli), ylim = c(-yli, yli)) +
                theme_bw() +
                theme(legend.position = "none",
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                theme(legend.position = "none")
        ggsave(paste0(resdir, name1, "-", name2, ".pdf"), width = 4, height = 4)
        
        return(dat)
}



# chromatinOrganization = c("ADA2", "APC9", "ARP8", "ASF1", "BRE1", "BRE2", "CAC2", 
#                           "CBF1", "CDC73", "CHD1", "CPR1", "CTI6", "CTR9", "DOC1", 
#                           "DPB4", "EAF1", "ELP4", "FKH1", "GCN5", "HDA1", "HDA1", 
#                           "HDA2", "HDA3", "HFI1", "HOS2", "HOS2", "HOS3", "HPC2", 
#                           "HST3", "IES2", "IOC4", "ISW1", "ISW2", "ITC1", "LEO1", 
#                           "LGE1", "MEC3", "MGA2", "MRC1", "MSH2", "MSI1", "MSN2", 
#                           "NAP1", "NGG1", "NUP133", "NUP170", "NUP60", "RAD54", 
#                           "REG1", "RSC1", "RSC2", "RTT106", "SEM1", "SET1", "SET2", 
#                           "SET3", "SET4", "SGF11", "SGF29", "SIF2", "SIF2", "SNF2", 
#                           "SNF5", "SNT1", "SNT1", "SPT21", "SPT3", "SPT7", "SPT8", 
#                           "SWD1", "SWI3", "TUP1", "VPS71", "YAF9")

litCompScatter(data = all.comb, 
               r1 = "log2fc.schubert", r2 = "log2fc.kemmeren",
               p1 = "q.schubert", p2 = "q.kemmeren", 
               name1 = "This study", name2 = "Kemmeren", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = NULL)

litCompScatter(data = all.comb, 
               r1 = "log2fc.schubert", r2 = "log2fc.oduibhir",
               p1 = "q.schubert", p2 = "q.oduibhir", 
               name1 = "This study", name2 = "Oduibhir", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = c("Opposite"))

litCompScatter(data = all.comb, 
               r1 = "log2fc.schubert", r2 = "log2fc.reimand",
               p1 = "q.schubert", p2 = "q.reimand", 
               name1 = "This study", name2 = "Reimand", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = c("Consistent", "Opposite"))

litCompScatter(data = all.comb, 
               r1 = "log2fc.schubert", r2 = "log2fc.stefely",
               p1 = "q.schubert", p2 = "q.stefely", 
               name1 = "This study", name2 = "Stefely", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = c("Consistent", "Opposite"))

litCompScatter(data = all.comb, 
               r1 = "log2fc.kemmeren", r2 = "log2fc.oduibhir",
               p1 = "q.kemmeren", p2 = "q.oduibhir", 
               name1 = "Kemmeren", name2 = "Oduibhir", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = NULL)

litCompScatter(data = all.comb, 
               r1 = "log2fc.kemmeren", r2 = "log2fc.reimand",
               p1 = "q.kemmeren", p2 = "q.reimand", 
               name1 = "Kemmeren", name2 = "Reimand", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = c("Consistent", "Opposite"))

litCompScatter(data = all.comb, 
               r1 = "log2fc.kemmeren", r2 = "log2fc.stefely",
               p1 = "q.kemmeren", p2 = "q.stefely", 
               name1 = "Kemmeren", name2 = "Stefely", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = NULL)




# Mediator heatmap

mediator <- read_delim("annotations/Mediator.tsv", delim = "\t") %>% print()

all.comb %>% 
        pivot_longer(-c(gene, protein),
                     names_to = c(".value", "study"),
                     names_sep = "\\.") %>%
        filter(q < 0.05 & abs(log2fc) > 0.5) %>%
        filter(study %in% c("schubert", "kemmeren")) %>%
        inner_join(mediator, by = "gene") %>%
        mutate(gene    = factor(gene2, levels = rev(pull(mediator, gene2))),
               set     = factor(set, levels = c("Cdk8", "Middle", "Head", "Tail")),
               protein = factor(protein, levels = c("HTB2", "TDH1", "TDH2", "TDH3",
                                                    "ENO2", "FAS1", "FAS2", "RNR2",
                                                    "RPL9A", "SSA1", "YHB1"))) %>%
        ggplot(aes(x = protein, y = gene, fill = log2fc)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", name = "Log2FC") +
        facet_grid(set ~ study, space = "free", scales = "free_y") +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(resdir, "mediator_heatmap.pdf"), width = 4, height = 4)
        




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-16_SessionInfo.txt")

