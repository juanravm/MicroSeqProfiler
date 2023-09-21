#!/usr/bin/Rscript
#' LEfSe analysis
#'
#' @param species_fp species.tsv file path
#' 
#' @param metadata_fp Metadata file path 
#' 
#' @param  sampling Minimum sampling depth for LEfSe analysis
#' 
#' @param col Metadata column name of compared groups
#' 
#' @param ref_group Reference group for LEfSe analysis
#' 
#' @param group2 Second group in LEfSe analysis
#'  
#' @param minimum_LDA Minimum LDA score for showing in .tsv output
#'  
#' @param output_dir LEfSe directory path for output location
#'  
#' @return LEfSe_output.tsv Output with significant microorganisms
#' LDA score 

library(lefser)


args <- commandArgs(trailingOnly = TRUE)

## Variables definition
species_fp <- NULL
metadata_fp <- NULL
sampling <- NULL
ref_group <-NULL
col <- NULL
group2 <- NULL
minimum_LDA <- NULL
output_dir <- NULL

## Argument assign
for (i in seq_along(args)) {
  if (args[i] == "--species_fp") {
    species_fp <- args[i + 1]
  } else if (args[i] == "--metadata_fp") {
    metadata_fp <- args[i + 1]
  } else if (args[i] == "--sampling") {
    sampling <- args[i + 1]
  } else if (args[i] == "--ref_group") {
    ref_group <- args[i + 1]
  } else if (args[i] == "--col") {
    col <- args[i + 1]
  } else if (args[i] == "--group2") {
    group2 <- args[i + 1]
  } else if (args[i] == "--minimum_LDA") {
    minimum_LDA <- args[i + 1]
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
}
}

## Importing data
comparison <- read.delim(species_fp, header=T, row.names = 1)
colnames(comparison) <- gsub("\\.", "-", colnames(comparison))
colnames(comparison) <- gsub("^X", "", colnames(comparison))

metadata <- read.delim(metadata_fp, comment.char="#")

# - Removing samples with low sampling depth such as in diversity analysis
sampling_depth <- colSums(comparison)
sampling_depth <- sampling_depth > sampling
comparison <- comparison[, sampling_depth]

# Metadata information
group <- metadata[, col][match(colnames(comparison), metadata[,1])]
na <- !is.na(group)
comparison <- comparison[,na]
group <- group[na]
group <- factor(group, levels = c(ref_group, group2))

filter <- !is.na(group)
comparison <- comparison[,filter]
group <- group[filter]

# Per sample normalization to 1M reads (as Galaxy does)
colsums<-colSums(comparison)

for (i in 1:length(colsums)) {
  comparison[, i]<- (comparison[, i]/colsums[i])*1000000
}

comparison <- as.matrix(comparison)

# Modifying input format
comparison <- SummarizedExperiment(assays = SimpleList(counts = comparison),
                               colData = DataFrame(grupo = group))

## LefSe calculation
set.seed(324514)

Lefse<-lefser(
  expr = comparison,
  kruskal.threshold = 0.05,
  wilcox.threshold = 0.05,
  lda.threshold = minimum_LDA,
  groupCol = "grupo",
  blockCol = NULL,
  assay = 1,
  trim.names = F,
  checkAbundances = F
)

write.table(Lefse, file = paste(output_dir, "/LEfSe_output.tsv", sep=""), 
            sep = "\t", 
            quote = FALSE, 
            row.names = F, 
            col.names = T)
            