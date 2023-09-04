#!/usr/bin/Rscript
#' OTU decontamination (Negative control normalization)
#'
#' @param OTU_filtered_table OTU filtered table file path with 
#' nonchimeric counts
#' 
#' @param neg A character pattern to identify negative controls in
#' OTU table rownames
#' 
#' @param  output_dir Decontaminated OTU table output directory
#' (usually the previous microbiome working directory employed)
#' 
#' @param metadata Metadata file path
#' 
#' @return decontam_OTU_table.csv QIIME2 artifact with the chimera-filtered 
#' and decontaminated OTU counts

args <- commandArgs(trailingOnly = TRUE)

# Variables definition
OTU_filtered_table <- NULL
neg <- NULL
metadata_fp <- NULL
output_dir <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--OTU_filtered_table") {
    OTU_filtered_table <- args[i + 1]
  } else if (args[i] == "--neg") {
    neg <- args[i + 1]
  } else if (args[i] == "--metadata_fp") {
    metadata_fp <- args[i + 1]
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
  } 
}

print(paste("OTU table file path:", OTU_filtered_table))
print(paste("Output directory:", output_dir))
print(paste("Metadata file path:", metadata_fp))
print(paste("Negative control identification pattern:", neg))

library(decontam)

## Importing OTU table
OTU_table <- t(read.delim(OTU_filtered_table, 
header=FALSE, 
comment.char=""))

rownames(OTU_table) <- OTU_table[,1]
OTU_table <- OTU_table[,-1]
colnames(OTU_table) <- OTU_table[1,]
OTU_table <- OTU_table[-1,]

rownames <- rownames(OTU_table)
OTU_table <- apply(OTU_table, 2, function(x) as.numeric(x))
rownames(OTU_table) <- rownames


## Negative logical vector for controls localization
neg <- grepl(neg, rownames(OTU_table))

## Logical vector indicating contaminants
Contaminant <- as.logical(isContaminant(
  seqtab = OTU_table,
  conc = NULL,
  neg = neg,
  method = "auto",
  batch = NULL,
  threshold = 0.1,
  normalize = FALSE,
  batch.combine = "minimum",
  detailed = F
))

## Removing frequencies from contaminated microorganisms
OTU_table[,Contaminant]=0

## Exporting decontaminated table without experimental controls
noncontam_export <- OTU_table[!neg ,]
write.csv(t(noncontam_export),
          row.names = T, 
          col.names = T, 
          file = paste(output_dir, "/decontam_OTU_table.csv", sep = ""))

