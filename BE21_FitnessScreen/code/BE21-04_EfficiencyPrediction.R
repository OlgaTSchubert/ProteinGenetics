# Libraries --------------------------------------------------------------------

library(tidyverse)
library(glmnet)




# Definitions ------------------------------------------------------------------

# Working directory
setwd("./BE21_FitnessScreen")
getwd()


# Functions
source("code/BE21-00_Functions.R")


# Results directory
resdir <- "results/efficiencyPrediction/"
dir.create(resdir)


# Guides table
guides <- readRDS("../BE00_gRNALibraryDesign/results/guides.RDS") %>%
    filter(set %in% c("eStops", "Ctr1", "Ctr2", "Ctr3")) %>%
    filter(!is.na(guide30)) %>%
    select(guide, guide30, set) %>%
    print()


# Fitness screen results
res <- readRDS("results/fitness_wES/wES_Ratios_gd_long.RDS") %>% 
    filter(set %in% c("eStops", "Ctr2")) %>%
    filter(timepoint == "T24h") %>%
    filter(guide %in% guides$guide) %>%
    print()




# Featurize guides -------------------------------------------------------------

# Generate matrix with a row per guide and a column for each sequence feature (takes ~1h to run)
feats <- featurization(guides$guide30, c("A","T","C","G"), seqorder = 4, posorder = 4)
dim(feats)
feats[1:10, 1:10]

saveRDS(feats, file = paste0(resdir, "featureMatrix.RDS"))
#feats <- readRDS(paste0(resdir, "featureMatrix.RDS"))


# Convert to tibble and add guide sequences
features <- as_tibble(feats) %>%
    mutate(guide = guides$guide, .before = "A") %>%
    filter(guide %in% res$guide) %>%
    print()




# Define test and training sets ------------------------------------------------

set.seed(1)
testids     <- sort(sample(1:nrow(features), 1000))
train       <- res[!(1:nrow(res) %in% testids), ]
test        <- res[1:nrow(res) %in% testids, ]
train.feats <- features[!(1:nrow(res) %in% testids), ]
test.feats  <- features[ 1:nrow(res) %in% testids, ]

# Binarize data for logistic lasso
train <- train %>% mutate(dropout = ifelse(train$log2fc < -0.5, T, F)) %>% print()
test  <- test  %>% mutate(dropout = ifelse(test$log2fc < -0.5, T, F)) %>% print()




# Ordinary lasso, continuous outcome -------------------------------------------

# Train ordinary lasso model on log2fc
olasso <- cv.glmnet(x = data.matrix(train.feats), y = train$log2fc)
saveRDS(olasso, file = paste0(resdir, "lasso_model_ordinary.RDS"))
olasso
plot(olasso)


# Keep only features with non-zero beta
olasso_coef <- tibble(motif = names(coef(olasso)[, 1]), beta = coef(olasso)[, 1]) %>%
    filter(motif != "(Intercept)",
           motif != "guide") %>%
    arrange(-abs(beta)) %>%
    filter(beta != 0) %>% print()
write_delim(olasso_coef, paste0(resdir, "olasso_coef.tsv"), delim = "\t")


# Test ordinary lasso model on log2fc
olasso.p <- predict(olasso, newx = data.matrix(test.feats))
cor.test(olasso.p, test$log2fc)$estimate^2  # 0.37




# Logistic lasso, binary outcome (classifier) ----------------------------------

# Train logistic lasso model
llasso <- cv.glmnet(x = data.matrix(train.feats), y = train$dropout, 
                    family = "binomial", type.measure = "auc")
saveRDS(llasso, file = paste0(resdir, "lasso_model_logistic.RDS"))
llasso
plot(llasso)


# Keep only features with non-zero beta
llasso_coef <- tibble(motif = names(coef(llasso)[, 1]), beta = coef(llasso)[, 1]) %>% 
    filter(motif != "(Intercept)",
           motif != "guide") %>%
    arrange(-abs(beta)) %>%
    filter(beta != 0) %>% print()
write_delim(llasso_coef, paste0(resdir, "llasso_coef.tsv"), delim = "\t")


# Test logistic lasso model
llasso.p <- predict(llasso, newx = data.matrix(test.feats))
cor.test(llasso.p, test$log2fc)$estimate^2  # 0.35


# Contingency table
test2 <- test %>%
    mutate(llassoPred = as.vector(llasso.p)) %>%
    mutate(llassoBinary = ifelse(llassoPred > 0, T, F)) %>% print()
    
ctab <- table(test2$dropout, test2$llassoBinary) %>% print()
#       FALSE TRUE
# FALSE   320  132
# TRUE    127  421


# Chi-square test
chisq.test(ctab)$p.value  # 6e-51


# Odds ratio
(ctab[1,1] / ctab[1,2]) / (ctab[2,1] / ctab[2,2])  # 8.0




# Plotting ---------------------------------------------------------------------

test %>%
    ggplot(aes(x = log2fc, y = olasso.p)) +
    geom_point() +
    labs(x = "Log2FC", y = "Lasso model") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0(resdir, "olasso_prediction.pdf"), height = 3, width = 3)


test %>%
    ggplot(aes(x = dropout, y = llasso.p)) +
    geom_jitter() +
    geom_boxplot(alpha = 0.8) +
    labs(x = "Dropout", y = "Lasso model") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0(resdir, "llasso_prediction.pdf"), height = 3, width = 3)



olasso_coef
llasso_coef




# Session info -----------------------------------------------------------------

writeLines(capture.output(devtools::session_info()), "code/BE21-04_SessionInfo.txt")

