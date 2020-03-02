library(ranger)
library(dplyr)
library(readr)
library(feather)
library(Pomona)

set.seed(1)

# perform variable selection within random forests on fc minhashes (k 31,
# scaled 2000, with hashes that were present only once across all samples
# removed). Uses the Pomona package vita implementation, which wraps the ranger
# package. Saves output as RDS for faster loading of output objects into
# subsequent R sessions.


## read in data ------------------------------------------------------------

## format hash table (samples x features)

fc <- read_feather(snakemake@input[['feather']]) # read in hash abund table
#fc <- read_feather("outputs/hash_tables/normalized_abund_hashes_wide.feather")
fc <- as.data.frame(fc)                         # transform to dataframe
fc$sample <- gsub("_filt_named\\.sig", "", fc$sample)
rownames(fc) <- fc$sample                       # set sample as rownames
fc <- fc[ , -ncol(fc)]                         # remove the samples column

## read in study metadata
info <- read_tsv(snakemake@input[['info']])
#info <- read_tsv("inputs/working_metadata_fecal_calprotectin.tsv")
## set validation cohorts (e.g. rm time series from training and testing)
info_validation <- info %>%
  filter(set == "valid") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
fc_validation <- fc[rownames(fc) %in% info_validation$library_name, ]

info_train <- info %>%
  filter(set == "train") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
fc_train <- fc[rownames(fc) %in% info_train$library_name, ]

info_test <- info %>%
  filter(set == "test") %>%
  mutate(library_name = gsub("-", "\\.", library_name))
fc_test <- fc[rownames(fc) %in% info_test$library_name, ]


## make classification vector
## match order of to fc
info_train <- info_train[match(rownames(fc_train), info_train$library_name), ]
## make diagnosis var
calprotectin_train <- info_train$fecal_calprotectin_log

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
fc_vita <- var.sel.vita(x = fc_train, y = calprotectin_train, p.t = 0.05,
                         ntree = 5000, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = 10, method = "ranger",
                         type = "regression")
saveRDS(fc_vita, snakemake@output[["vita_rf"]])

# write files -------------------------------------------------------------

## write predictive hashes
var <- fc_vita$var                 # separate out selected predictive hashes
var <- gsub("X", "", var)           # remove the X from the beginning of hashes
write.table(var, snakemake@output[['vita_vars']],
            quote = F, col.names = F, row.names = F)
# write.table(var, "outputs/vita_rf/vita_vars.txt",
#             quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write training/testing set (novalidation) to files
fc_filt <- fc[ , colnames(fc) %in% var] # subset fc to hashes in fc_vita
write.csv(fc_filt, snakemake@output[['fc_filt']], quote = F)
# write.csv(fc_filt, "outputs/vita_rf/fc_filt.csv", quote = F)

