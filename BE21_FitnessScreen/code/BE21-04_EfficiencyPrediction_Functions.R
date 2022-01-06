featurization <- function(sequences, string, seq = TRUE, seqorder = 2, pos = TRUE, posorder = 2) {
        
        # The functions countpattern and findpositions are separate below

        total    <- 0
        
        features <- data.frame(1:length(sequences))
        
        colnames(features)[length(features)] <- "Serial"
        
        if (seq == TRUE) {
                for (s in 1:seqorder) {
                        permu = gtools::permutations(n = length(string), r = s, v = string, repeats.allowed = TRUE)
                        for (i in 1:length(permu[, 1])) {
                                temp = countpattern(sequence = sequences, pattern = paste(permu[i,], collapse = ""))
                                if (sum(temp) > 0) {
                                        features = data.frame(features, temp)
                                        colnames(features)[length(features)] = paste(permu[i,], collapse = "")
                                }
                                total = total + 1
                        }
                        cat(s, " order seq. features:", length(permu[, 1]), ":total features = ", total, "\n")
                }
        }
        
        if (pos == TRUE) {
                minlength = min(unlist(lapply(sequences, function(s) { nchar(toString(s)) })))
                for (p in 1:posorder) {
                        permu = gtools::permutations(n = length(string), r = p, v = string, repeats.allowed = TRUE)
                        ps = 0
                        for (i in 1:(length(permu[, 1]))) {
                                for (j in 1:(minlength - length(permu[i, ]) + 1)) {
                                        temp = findposition(sequence = sequences, pattern = paste(permu[i,], collapse = ""), j)
                                        if (sum(temp) > 0) {
                                                features = data.frame(features, temp)
                                                colnames(features)[length(features)] = paste0(paste(permu[i,], collapse = ""), "_", j)
                                        }
                                        ps = ps + 1
                                        total = total + 1
                                }
                        }
                        cat(p, "order pos. features:", ps, ":total features = ",  total, "\n")
                }
        }
        
        features$Serial <- NULL

        return(features)
}




findposition <- function(sequence, pattern, position) {
        # Function from CRISPRpred, Rahman and Rahman, PLOS ONE 2017
        unlist(lapply(sequence,function(s)
                if(position %in% gregexpr(pattern = pattern, toString(s))[[1]]){1}else{0}
        ))
}




countpattern <- function(sequence, pattern) {
        # Function from CRISPRpred, Rahman and Rahman, PLOS ONE 2017
        unlist(lapply(sequence,function(s)
                if(-1 %in% gregexpr(pattern = pattern, toString(s))[[1]]){0}else{length(
                        gregexpr(pattern = pattern, toString(s))[[1]])}
        ))
}




