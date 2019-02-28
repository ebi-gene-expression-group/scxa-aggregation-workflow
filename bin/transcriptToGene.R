#!/usr/bin/env Rscript

library(rtracklayer)

cl <- commandArgs(trailingOnly = TRUE)

gtf_file <- cl[1]
transcript_attribute <- cl[2]
gene_attribute <- cl[3]
output_file <- cl[4]

# Read annotation to set up transcript/gene relationships

annotation <- elementMetadata(import(gtf_file))
tx2gene <- unique(annotation[ ! is.na(annotation[[transcript_attribute]]), c(transcript_attribute, gene_attribute)] )

# Check for inconsistencies

transcript_id_counts <- table(tx2gene[[transcript_attribute]])

if (any(any(transcript_id_counts > 1))){
  write(paste("Some transcripts have more than one gene mapping:", paste(names(transcript_id_counts[transcript_id_counts > 1]), collapse=', ')), stderr())
  q(status = 1)
}

# Write the outputs

write.csv(tx2gene, file=output_file, row.names=FALSE)
