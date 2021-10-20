# Functions --------------------------------------------------------------------

readCountsLong <- function(counts.file, guide.info) {
        
        # Function to convert read counts file into long format with some
        # additional columns from guides.RDS
        
        # Read in data (1 row per guide-barcode combination)
        dat <- readRDS(counts.file)
        
        # Add more information (join command too slow)
        reads.w <- as_tibble(cbind(guide.info[match(dat$guide, guide.info$guide),], dat[, -1])) %>%
                mutate(BCID = paste0(.$guide, '_', .$barcode))
        
        # Convert to long format
        reads <- reads.w %>%
                gather(names(dat)[-c(1,2)], key = "sample", value = "n_reads") %>%
                separate(sample, into = c("tail", "repl"), remove = T) %>%
                unite(tail, repl, col = "sample", sep = ".", remove = F) %>%
                arrange(sample, BCID)
        
        return(reads)
}



# ------------------------------------------------------------------------------

glmStats <- function(BCCs = BCcounts) {
        
        # Evaluate for each guide if it is significantly different in the high 
        # vs the low GFP sample using a generalized linear model (glm) framework.
        
        # Load libraries
        library(data.table)
        library(broom)
        
        # Determine number of replicates
        reps <- length(unique(BCCs$repl))
        
        # Convert sample annotations to factors and set "un" as control
        BCCs$repl <- as.factor(BCCs$repl)
        BCCs$tail <- as.factor(BCCs$tail)
        BCCs$tail <- relevel(BCCs$tail, "un")
        
        # Build model matrix for hi-vs-lo
        gseq <- unique(BCCs$guide)[1]
        obs  <- BCCs[BCCs$guide == gseq, ]
        mm   <- model.matrix(~ obs$repl + obs$tail)
        colnames(mm) <- str_replace_all(colnames(mm), "obs\\$", "")
        mm[mm[, (reps+2)] == 1, (reps+1)] <- -1
        mm <- mm[, -(reps+2)]
        colnames(mm)[(reps+1)] <- "hilo"
        
        # Define offset (normalization factor) for model
        total.counts <- BCCs %>%
                group_by(sample) %>%
                summarize(totcounts = sum(nBCs))
        
        # Initiate empty list to collect results for each guide
        statsHiLo   <- list()
        
        # Intitiate progress bar
        pb <- txtProgressBar(min = 1, max = length(unique(BCCs$guide)), style = 3)

        # Loop through data, one guide at a time
        for(gID in seq_along(unique(BCCs$guide))) {

                setTxtProgressBar(pb, gID)
                
                gseq <- unique(BCCs$guide)[gID]
                obs  <- BCCs[BCCs$guide == gseq, ]
                
                # Use poisson model for error distribution and log as link function
                # For normalization across samples: offset(log())
                # To remove intercept: -1 (needed if intercept already in mm)
                
                modelHiLo <- glm(nBCs ~ mm -1 + offset(log(total.counts$totcounts)), 
                                 family = poisson(link = log), data = obs)
                statsHiLo[[as.character(gseq)]] <- tidy(modelHiLo)[(-1:-reps), ]
        }
        close(pb)
        
        statsHiLo   <- rbindlist(statsHiLo, idcol = "gseq")
        statsHiLo   <- split(statsHiLo, statsHiLo$term)
        statslist   <- list(HiLo = statsHiLo$mmhilo)
        
        return(statslist)
}

