setwd("~/github/2020-calprotectin/")
set.seed(1)

# make fecal calprotectin metadata

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)

files <- list.files("inputs/ENA", full.names = T)
# filter to cohorts with calprotectin measurements
files <- files[files %in% c('inputs/ENA/PRJNA385949.txt', "inputs/ENA/PRJNA400072.txt", "inputs/ENA/SRP057027.txt")]
prj <- files %>%
  map(read_tsv) %>%
  reduce(rbind) %>%
  filter(library_strategy != "AMPLICON") %>%
  filter(library_strategy != "RNA-Seq") %>%
  filter(read_count >= 1000000) %>%
  filter(library_layout  == "PAIRED")

table(prj$study_accession)
# PRJNA385949 time series
# SRP057027 time series
# pediatric crohn's data set: https://raw.githubusercontent.com/louiejtaylor/sbx_lewis2015/master/metadata/metadata_with_SRR.csv
# PRJNA400072 matches library name

# join with study-specific metadata
# PRJNA400072 -- library_name
# SRP057027 -- run_accession, library_name
# PRJNA385949 -- library_name

prjna400072 <- read_csv("inputs/metadata/PRJNA400072_metadata.csv")
prj <- left_join(prj, prjna400072, by = "library_name") 
srp057027 <- read_csv("inputs/metadata/SRP057027_metadata.csv")
prj <- left_join(prj, srp057027, by = c("run_accession", "library_name"))
prjna385949 <- read_csv("inputs/metadata/PRJNA385949_metadata.csv")
prj <- left_join(prj, prjna385949, by = "library_name")


# coalesce colnames
prj <- prj %>%
  mutate(diagnosis = coalesce(diagnosis.x, diagnosis.y, diagnosis)) %>%
  select(-diagnosis.x, -diagnosis.y) %>%
  mutate(antibiotic = coalesce(antibiotic.x, antibiotic.y)) %>%
  select(-antibiotic.x, -antibiotic.y) %>%
  mutate(steroids = coalesce(steroids.x, steroids.y)) %>%
  select(-steroids.x, -steroids.y) %>%
  mutate(subject = coalesce(as.character(subject), as.character(patient_id))) %>%
  select(-patient_id) %>%
  mutate(fecal_calprotectin = coalesce(fecal_calprotectin.x, as.character(fecal_calprotectin.y), as.character(FCP))) %>%
  select(-fecal_calprotectin.x, -fecal_calprotectin.y, -FCP) %>%
  filter(diagnosis != "IC") # also removes NAs

# Filter NAs for fecal calprotectin
prj <- prj %>%
  filter(fecal_calprotectin != "#N/A") %>%
  mutate(fecal_calprotectin = as.numeric(fecal_calprotectin))

ggplot(prj, aes(x = fecal_calprotectin, fill = study_accession)) +
  geom_density(alpha = .3) +
  theme_minimal()

ggplot(prj, aes(x = fecal_calprotectin, fill = diagnosis)) +
  geom_density(alpha = .3) +
  theme_minimal()

table(prj$diagnosis)

prj %>% group_by(study_accession, diagnosis) %>% tally
wrk <- select(prj, study_accession, run_accession, library_name, read_count, 
              sample_alias, diagnosis, subject, fecal_calprotectin)

wrk$subject <- ifelse(is.na(wrk$subject), wrk$library_name, wrk$subject)

wrk <- wrk[order(wrk$library_name, decreasing = F), ]
wrk <- wrk[!duplicated(wrk$subject), ]

length(unique(wrk$library_name))
length(unique(wrk$subject))

# combine with hmp metadata -----------------------------------------------

hmp <- read_tsv("inputs/hmp2_mgx_metadata.tsv") %>%
  mutate(study_accession = "iHMP") %>%
  mutate(run_accession = NA) %>%
  mutate(library_name = External.ID) %>%
  mutate(read_count = reads_raw) %>%
  mutate(sample_alias = NA) %>%
  mutate(subject = Participant.ID) %>%
  mutate(fecal_calprotectin = fecalcal) %>%
  arrange(subject, week_num) %>%
  select(study_accession, run_accession, library_name, read_count, sample_alias,
         diagnosis, subject, fecal_calprotectin) %>%
  filter(!is.na(fecal_calprotectin))

# remove time series samples
hmp <- hmp[!duplicated(hmp$subject), ]

wrk <- rbind(wrk, hmp)

length(unique(wrk$subject))


# set train, test, and valid ----------------------------------------------

wrk %>%
  group_by(study_accession, diagnosis) %>%
  tally()

wrk %>% 
  filter(study_accession == "PRJNA385949") %>%
  distinct(subject) %>%
  nrow()

wrk %>% 
  filter(study_accession == "iHMP") %>%
  distinct(subject) %>%
  nrow()

# set validation cohort as PRJNA385949
wrk$set <- ifelse(wrk$study_accession == "PRJNA385949", "valid", NA)

# Determine testing and training sets. No subject should be in both "train" and
# "test". 
# set train/test based on subject identifier

subjects <- wrk %>%
  filter(study_accession != "PRJNA385949") %>%
  distinct(subject)

subjects
train <- sample(nrow(subjects), 0.7*nrow(subjects), replace = FALSE)
train_set <- subjects[train, ]
test_set <- subjects[-train, ]

for(i in 1:nrow(wrk)){
  if(wrk$subject[i] %in% train_set$subject){
    wrk$set[i] <- "train"
  } else if(wrk$subject[i] %in% test_set$subject){
    wrk$set[i]  <- "test"
  } else {
    wrk$set[i] <- "valid"
  }
}

write_tsv(x = wrk, path = "inputs/working_metadata_fecal_calprotectin.tsv")
