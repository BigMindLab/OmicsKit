# Metadata - TCGA LUAD -------------------------

# August, 2025
# Source: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_counts.tsv.gz

sampledata <- read.csv("../TCGA-LUAD.samples.reduced.tsv", sep = "\t", header = TRUE)
sampledata$patient_id <- gsub("-", ".", sampledata$patient_id)

usethis::use_data(sampledata, compress = "xz", overwrite = TRUE)
