# STAR Counts - TCGA LUAD -------------------------

# August, 2025
# Source: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_counts.tsv.gz

raw_counts <- read.delim("data-raw/TCGA-LUAD.star_counts.tsv.gz", sep = "\t", header = TRUE)
sampledata <- read.csv("../TCGA-LUAD.samples.reduced.tsv", sep = "\t", header = TRUE)

raw_counts[, -1] <- as.matrix(raw_counts[, -1])

gids <- sub("\\.\\d+$", "", raw_counts$Ensembl_ID)
raw_counts <- raw_counts[,-1]
rownames(raw_counts) <- gids

raw_counts <- 2^raw_counts - 1

sampledata$patient_id <- gsub("-", ".", sampledata$patient_id)
raw_counts <- raw_counts[, sampledata$patient_id]

usethis::use_data(raw_counts, compress = "xz", overwrite = TRUE)
