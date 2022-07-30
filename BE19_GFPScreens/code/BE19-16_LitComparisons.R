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
        left_join(stefely.comb, by = c("protein", "gene"), suffix = c(".schubert", ".stefely")) %>%
        left_join(kemmeren.comb, by = c("protein", "gene")) %>%
        rename(log2fc.kemmeren = log2fc, q.kemmeren = q) %>%
        left_join(oduibhir.comb, by = c("protein", "gene")) %>%
        rename(log2fc.oduibhir = log2fc, q.oduibhir = q) %>%
        left_join(reimand.comb, by = c("protein", "gene")) %>%
        rename(log2fc.reimand = log2fc, q.reimand = q) %>%
        #pivot_longer(-c(gene, protein), names_to = c(".value", "study"),
        #             names_pattern = "(.*)\\.(.*)") %>%
        print()




# Plotting ---------------------------------------------------------------------

litCompScatter <- function(data, r1, r2, p1, p2, name1, name2, 
                           p = 0.05, r = 0.5, xli = 4, yli = 4, 
                           labels = c("Consistent", "Opposite")) {
        
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
                           aes(x = r1, y = r2), color = "red3") +
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
        
        ggsave(paste0(resdir, name1, "-", name2, ".pdf"), width = 4, height = 4)
        
        return(dat)
}


litCompScatter(data = all.comb, 
               r1 = "log2fc.schubert", r2 = "log2fc.kemmeren",
               p1 = "q.schubert", p2 = "q.kemmeren", 
               name1 = "This study", name2 = "Kemmeren", 
               p = 0.05, r = 0.5, xli = 4, yli = 4, 
               labels = c("Opposite"))

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




# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "code/BE19-16_SessionInfo.txt")

