# Import external data

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
                
                return(ds)
        
                
                
        } else if(dataset == "txdb") {
                
                if(localfile == FALSE) {
                        
                        url          <- "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz"
                        zip          <- "~/Desktop/saccharomyces_cerevisiae.gff.gz"
                        file         <- "~/Desktop/saccharomyces_cerevisiae.gff"
                        download.file(url, zip, method = "auto", quiet = FALSE, mode = "w")
                        R.utils::gunzip(filename = zip, destname = file)                        
                        
                } else { file <- "../ExternalData/SGD_20180426/saccharomyces_cerevisiae.gff" }
                
                ds   <- GenomicFeatures::makeTxDbFromGFF(file)
        }

}

