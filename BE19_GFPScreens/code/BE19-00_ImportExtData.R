# Function for import of external data -----------------------------------------

# Function calls:
# source("code/BE19-00_importExtData.R")
# importExtData(dataset = "SGD_features",             localfile = T) %>% print()
# importExtData(dataset = "SGD_chromlengths",         localfile = T) %>% print()
# importExtData(dataset = "SGD_proteins",             localfile = T) %>% print()
# importExtData(dataset = "SGD_pathways",             localfile = F) %>% print()
# importExtData(dataset = "txdb",                     localfile = T) %>% print()
# importExtData(dataset = "txdb0",                    localfile = T) %>% print()
# importExtData(dataset = "BYxRM_variants",           localfile = T) %>% print()
# importExtData(dataset = "Albert2014_pQTLs_SacCer3", localfile = T) %>% print()
# importExtData(dataset = "Albert2018_eQTLs",         localfile = T) %>% print()
# importExtData(dataset = "Albert2018_hotspots",      localfile = T) %>% print()
# importExtData(dataset = "Newman2006",               localfile = T) %>% print()
# importExtData(dataset = "Ho2018",                   localfile = T) %>% print()
# importExtData(dataset = "Kemmeren2014",             localfile = T) %>% print()
# importExtData(dataset = "ODuibhir2014",             localfile = T) %>% print()
# importExtData(dataset = "Reimand2010",              localfile = T) %>% print()
# importExtData(dataset = "Stefely2016",              localfile = T) %>% print()
# importExtData(dataset = "Costanzo2016",             localfile = T) %>% print()



importExtData <- function(dataset, localfile = FALSE) {
        
        
        if(dataset == "SGD_features") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
                        url.headers  <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.README"
                        file         <- "~/Desktop/SGD_features.tab"
                        file.headers <- "~/Desktop/SGD_features.README"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        download.file(url.headers, file.headers, method = "auto", quiet = FALSE, mode = "w")
                        
                } else {
                        
                        file         <- "../ExternalData/SGD_20180426/SGD_features.tab"
                        file.headers <- "../ExternalData/SGD_20180426/SGD_features.README"
                }
                
                ds      <- read_tsv(file, col_names = F) %>%
                        filter(X2 == "ORF") %>% 
                        filter(!is.na(X9))       # remove 2-micron ORFs
                headers <- read_tsv(file.headers, col_names = F) %>% 
                        filter(str_detect(X1, "\\d\\.\\s\\s")) %>%
                        mutate(X1 = str_extract(X1, "(?<=\\d\\.\\s\\s).+(?=\\s\\()"),
                               X1 = str_replace(X1, "^\\s", ""),
                               X1 = str_replace_all(X1, "\\s", "_")) %>% 
                        pull(X1)
                names(ds) <- headers
                
                
                
                
        } else if(dataset == "SGD_chromlengths") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/chromosome_length.tab"
                        url.headers  <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/chromosome_length.README"
                        file         <- "~/Desktop/chromosome_length.tab"
                        file.headers <- "~/Desktop/chromosome_length.README"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        download.file(url.headers, file.headers, method = "auto", quiet = FALSE, mode = "w")
                        
                } else {
                        
                        file         <- "../ExternalData/SGD_20180426/chromosome_length.tab"
                        file.headers <- "../ExternalData/SGD_20180426/chromosome_length.README"
                }
                
                ds      <- read_tsv(file, col_names = F)
                headers <- read_tsv(file.headers, col_names = F) %>%
                        filter(str_detect(X1, "\\d\\)\\s")) %>%
                        mutate(X1 = str_extract(X1, "(?<=\\d\\)\\s).+"),
                               X1 = str_replace_all(X1, "\\s", "_")) %>% 
                        pull(X1)
                names(ds) <- headers
                
                
                
                
        } else if(dataset == "SGD_proteins") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/calculated_protein_info/archive/protein_properties.20150126.tab.zip"
                        #url          <- "http://sgd-archive.yeastgenome.org/curation/calculated_protein_info/protein_properties.tab"
                        #url.headers  <- "http://sgd-archive.yeastgenome.org/curation/calculated_protein_info/protein_properties.README"
                        file         <- "~/Desktop/protein_properties.tab.zip"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        unzip(file, exdir = "~/Desktop/")
                        
                } else {
                        
                        file         <- "../ExternalData/SGD_20180426/protein_properties.tab"
                        file.headers <- "../ExternalData/SGD_20180426/protein_properties.README"
                }
                
                ds      <- read_tsv(file, col_names = T)
                
                
                
                
        } else if(dataset == "SGD_pathways") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/literature/biochemical_pathways.tab"
                        url.headers  <- "http://sgd-archive.yeastgenome.org/curation/literature/biochemical_pathways.README"
                        file         <- "~/Desktop/biochemical_pathways.tab"
                        file.headers <- "~/Desktop/biochemical_pathways.README"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        download.file(url.headers, file.headers, method = "auto", quiet = FALSE, mode = "w")
                        
                } else {
                        
                        file         <- "../ExternalData/SGD_20180426/biochemical_pathways.tab"
                        file.headers <- "../ExternalData/SGD_20180426/biochemical_pathways.README"
                }
                
                ds      <- read_tsv(file, col_names = F)
                headers <- read_tsv(file.headers, col_names = F) %>%
                        filter(str_detect(X1, "\\d\\)\\s")) %>% 
                        mutate(X1 = str_extract(X1, "(?<=\\d\\)\\s).+"),
                               X1 = str_replace_all(X1, "\\s", "_"),
                               X1 = str_replace_all(X1, "_\\(optional\\)", "")) %>%
                        pull(X1)
                names(ds) <- headers
                
                
                
                
        } else if(dataset == "txdb0") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz"
                        zip          <- "~/Desktop/saccharomyces_cerevisiae.gff.gz"
                        file         <- "~/Desktop/saccharomyces_cerevisiae.gff"
                        download.file(url, zip, method = "auto", quiet = FALSE, mode = "w")
                        R.utils::gunzip(filename = zip, destname = file)                        
                        
                        
                } else { file <- "../ExternalData/SGD_20180426/saccharomyces_cerevisiae.gff" }
                
                ds   <- GenomicFeatures::makeTxDbFromGFF(file)
                
                
                
                
        } else if(dataset == "txdb") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz"
                        zip          <- "~/Desktop/saccharomyces_cerevisiae.gff.gz"
                        file         <- "~/Desktop/saccharomyces_cerevisiae.gff"
                        download.file(url, zip, method = "auto", quiet = FALSE, mode = "w")
                        R.utils::gunzip(filename = zip, destname = file)                        
                        
                        
                } else { file <- "../ExternalData/SGD_20180426/saccharomyces_cerevisiae.gff" }
                
                ds   <- GenomicFeatures::makeTxDbFromGFF(file)
                ds   <- GenomicFeatures::genes(ds)
                
                
                
                
        } else if(dataset == "BYxRM_variants") {
                
                # GitHub of eLife 2018 paper: https://github.com/elifesciences-publications/eQTL_BYxRM/blob/master/code/additional_code_from_Frank/hotspotAnalyses/hotspotAnalyses_SVD_hotspots_170501_4repo.R
                # Downloaded from: https://drive.google.com/drive/folders/0B4EjgO02Xr8yN2lpMFpVbmNVWVk
                # Note that during import, quote = "" is important because of the gene "IMP2'" (' is part of the gene name)
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else { file <- "../ExternalData/RM_SNPs/gdata_42k_VEP_wExtraMarkersAttached.txt" }
                
                ds <- read_tsv(file, quote = "") %>%
                        dplyr::rename(Variant = "#Uploaded_variation") %>%
                        mutate(SYMBOL = str_replace(SYMBOL, "1-Oct", "OCT1"),
                               SYMBOL = str_replace(SYMBOL, "\"DUR1,2\"", "DUR1,2"))
                
                
                
                
        } else if(dataset == "Albert2014_pQTLs_SacCer3") {
                
                # Albert et al., Nature 2014
                # Genome coordinates converted from SacCer2 to SacCer3
                # File provided by Frank (email 7/5/2020)
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else {
                        
                        load("../ExternalData/Albert_Nature_2014/R_XpQTLDistantSacCer3.RData")
                        ds   <- as_tibble(XpQTLDistantsacCer3)
                        rm(XpQTLDistantsacCer3)
                }
                
                
                
                
        } else if(dataset == "Albert2018_eQTLs") {
                
                # Albert & Bloom et al., eLife 2018
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else { file <- "../ExternalData/Albert_eLife_2018/elife-35471-data4-v2.xlsx" }
                
                ds   <- readxl::read_xlsx(file)
                
                
                
                
        } else if(dataset == "Albert2018_hotspots") {
                
                # Albert & Bloom et al., eLife 2018
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else { file <- "../ExternalData/Albert_eLife_2018/elife-35471-data8-v2.xlsx" }
                
                ds   <- readxl::read_xlsx(file)
                
                
                
                
        } else if(dataset == "Newman2006") {
                
                # Newman et al., Nature 2006 - Yeast GFP collection abundances
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else { file <- "../ExternalData/Newman_Nature_2006/Newman_Nature_2006_TableS1.csv" }
                
                ds   <- readr::read_csv(file, col_names = T, na = c("","NA"))
                
                
                
                
        } else if(dataset == "Ho2018") {
                
                # Ho et al., Cell Systems 2018 - Proteome absolute abundances
                
                if(localfile == FALSE) {
                        
                        writeLines("There is no download information.")
                        
                } else { file <- "../ExternalData/Ho_CellSystems_2018/1-s2.0-S240547121730546X-mmc5.xlsx" }                
                
                ds   <- readxl::read_excel(file, skip = 2, col_names = T, na = c("","NA"))
                
                
                
                
        } else if(dataset == "Kemmeren2014") {
                
                # Kemmeren et al., Cell 2014 - YKO transcriptomics (microarray)
                # 1484 YKO mutants with (putative) role in regulation of gene expression.
                # Cultures grown in SC media to exponential phase.
                
                if(localfile == FALSE) {
                        
                        url  <- "https://deleteome.holstegelab.nl/data/downloads/deleteome_all_mutants_ex_wt_var_controls.txt"
                        file <- "~/Desktop/deleteome_all_mutants_ex_wt_var_controls.txt"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        
                } else { file <- "../ExternalData/Kemmeren_Cell_2014/deleteome_all_mutants_ex_wt_var_controls.txt" }                
                
                ds   <- readr::read_delim(file, delim = "\t")
                
                
                
                
        } else if(dataset == "ODuibhir2014") {
                
                # O'Duibhir et al., MSB 2014 - YKO transcriptomics (microarray)
                # Reanalysis of Kemmeren  et al., see above, corrected for ESR-like signature.
                
                if(localfile == FALSE) {
                        
                        url  <- "https://deleteome.holstegelab.nl/data/downloads/deleteome_all_mutants_svd_transformed.txt"
                        file <- "~/Desktop/deleteome_all_mutants_svd_transformed.txt"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        
                } else { file <- "../ExternalData/ODuibhir_MSB_2014/deleteome_all_mutants_svd_transformed.txt" }
                
                ds   <- readr::read_delim(file, delim = "\t")
                
                
                
                
        } else if(dataset == "Reimand2010") {
                
                # Reimand et al., NAR 2010 - YKO transcriptomics (microarray)
                # This is a re-analysis of data by Hu et al., Nat Methods 2006.
                # Targeted 269 non-essential yeast transcription factors.
                # Cultures grown in YPD media to mid-log phase.

                if(localfile == FALSE) {
                        
                        writeLines("Automatic download doesn't consistently work because link keeps changing.\nDownload manually from https://academic.oup.com/nar/article/38/14/4768/2409503\nUnzip and add header 'Name' to first column.")
                        
                } else {
                        
                        file.r.mod <- "../ExternalData/Reimand_NAR_2010/SuppTable1C_KOTargetGenes_Matrix_LogFoldChange_mod.txt"
                        file.p.mod <- "../ExternalData/Reimand_NAR_2010/SuppTable1B_KOTargetGenes_Matrix_PValue_mod.txt"
                        ds.r <- readr::read_delim(file.r.mod, delim = " ", trim_ws = F)
                        ds.p <- readr::read_delim(file.p.mod, delim = " ", trim_ws = F)
                        ds   <- list(ratios = ds.r, pvals = ds.p)
                        
                }
                
                
                
                
        } else if(dataset == "Stefely2016") {
                
                # Stefely et al., Nature Biotech 2016 - YKO proteomics
                # Focus on genes involved in mitochondrial processes.
                # Cultures grown in YPD ("fermentation condition") to exponential phase.
                
                if(localfile == FALSE) {
                        
                        url  <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3683/MediaObjects/41587_2016_BFnbt3683_MOESM21_ESM.xlsx"
                        file <- "~/Desktop/41587_2016_BFnbt3683_MOESM21_ESM.xlsx"
                        download.file(url, file, method = "auto", quiet = FALSE, mode = "w")
                        
                } else { file <- "../ExternalData/Stefely_NatBiotech_2016/41587_2016_BFnbt3683_MOESM21_ESM.xlsx" }                
                
                ds.r <- readxl::read_xlsx(file, sheet = "KO vs WT_Ferm (∆LFQ)", trim_ws = T)
                ds.p <- readxl::read_xlsx(file, sheet = "KO vs WT_Ferm (P-Value)", trim_ws = T)
                ds   <- list(ratios = ds.r, pvals = ds.p)
                
 
                
                
        } else if(dataset == "Costanzo2016") {
                
                # Downloaded from The Cell Map (https://thecellmap.org/)
                
                if(localfile == FALSE) {
                        
                        writeLines("Get interaction tables manually from The Cell Map (https://thecellmap.org/).
In the network view, enter a protein name (e.g., 'tdh1') and click 'View as Table'. Now, in the search field, add all other proteins.
Note that there are no results for Fas1 and Rnr2, and multiple results for Fas2.
When all proteins are added click 'Download all'.
Open the file and save as xlsx.")
                        
                } else {
                        
                        file <- "../ExternalData/Costanzo_Science_2016/tcm-ENO2_RPL9A_TDH1_fas2-5011_fas2-5012_YHB1_SSA1_HTB2_fas2-5022_fas2-5003_fas2-5001_fas2-5008_TDH2_TDH3-220730.xlsx"
                        
                        ds <- file %>% 
                                readxl::excel_sheets() %>% 
                                set_names() %>%
                                map(readxl::read_excel, path = file) %>%
                                keep(str_detect(names(.), "GI scores")) %>% 
                                map(~ select(., -`...6`)) %>%
                                map(~ rename_with(., ~ c("pos_ORF", "pos_allele", "pos_score", "pos_p", "pos_annot",
                                                         "neg_ORF", "neg_allele", "neg_score", "neg_p", "neg_annot"))) %>%
                                bind_rows(.id = "protein") %>% 
                                mutate(protein = str_replace(protein, " GI scores", "")) %>%
                                pivot_longer(-protein, names_to = c("direction", ".value"),
                                             names_sep = "_")
                }
                
                
                
        } else { writeLines( "There is no dataset with this name." ) }
        
        return(ds)
}

