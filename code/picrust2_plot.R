#!/bin/bash
#' OTU decontamination (Negative control normalization)
#'
#' @param ref_group Reference group for pathway differential
#' abundance analysis
#' 
#' @param col Metadata column where groups for differential analysis are
#' 
#' @param  metadata_fp Metadata file path
#' 
#' @param  method Select inference method for differential analysis
#' (KO, EC or MetaCyc)
#' 
#' @param  picrust_dir PICRUSt2 output directory
#' 
#' @param Visualizations PICRUSt2 Visualizations directory
#' 
#' @return file in .tsv format with pathway differential abundance analysis
#' 
#' @return ggpicrust2 visualization with differential abundance results

args <- commandArgs(trailingOnly = TRUE)

## Variables definition
ref_group <- NULL
col <- NULL
metadata_fp <- NULL
method <-NULL
picrust_dir <- NULL
Visualizations <- NULL

## Argument assign
for (i in seq_along(args)) {
  if (args[i] == "--ref_group") {
    ref_group <- args[i + 1]
  } else if (args[i] == "--col") {
    col <- args[i + 1]
  } else if (args[i] == "--metadata_fp") {
    metadata_fp <- args[i + 1]
  } else if (args[i] == "--method") {
    method <- args[i + 1]
  } else if (args[i] == "--picrust_dir") {
    picrust_dir <- args[i + 1]
  } else if (args[i] == "--Visualizations") {
    Visualizations <- args[i + 1]
  }
}

KO <- paste(picrust_dir, "/KO_metagenome_out/pred_metagenome_unstrat.tsv", sep = "")
MetaCyc <- paste(picrust_dir, "/pathways_out/path_abun_unstrat.tsv", sep = "")
EC <- paste(picrust_dir, "/EC_metagenome_out/pred_metagenome_unstrat.tsv", sep = "")


##······················ DIFFERENTIAL ABUNDANCE ANALYSIS
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggplot2)
library(MicrobiomeStat)
library(Maaslin2)
library(dplyr)
library(ggh4x)
library(GGally)

metadata <- as.data.frame(read_delim(metadata_fp, 
                                     comment="#", 
                                     escape_double = FALSE,
                                     trim_ws = TRUE))

colnames(metadata)[1] <- "sample_name"

metadata <- metadata[!is.na(metadata[,col]),]

if (method=="KO"){
  kegg_abundance <- ko2kegg_abundance(KO)
  intersect <- intersect(colnames(kegg_abundance),metadata$sample_name)
  metadata <- metadata[match(intersect,metadata$sample_name) ,]
  kegg_abundance <- kegg_abundance[, metadata$sample_name]
  
  group <- col
  kegg_daa_df <-
    pathway_daa(
      abundance = kegg_abundance,
      metadata = metadata,
      group = group,
      daa_method = "ALDEx2",
      select = NULL,
      p.adjust = "BH",
      reference = ref_group
    )
  
  if (length(unique(metadata[,col]))==2) {
    kegg_daa_sub_method_results_df <-
      kegg_daa_df[kegg_daa_df$method == "ALDEx2_Wilcoxon rank test", ]
  } else if (length(unique(metadata[,col]))>2) {
    kegg_daa_sub_method_results_df <-
      kegg_daa_df[kegg_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]
  }
  
  kegg_daa_annotated_sub_method_results_df <-
    pathway_annotation(pathway = "KO",
                       daa_results_df = kegg_daa_sub_method_results_df,
                       ko_to_kegg = TRUE)
  
  write.table(kegg_daa_annotated_sub_method_results_df, 
              file = paste(Visualizations, "/KO.tsv", sep = ""), 
              sep = "\t", 
              quote = FALSE, 
              row.names = F, 
              col.names = T)
  
  Group <- Group <-metadata [,col]
  
  # Removing NA in pathways name and rounding p_adjust to 5 decimals
  kegg_daa_annotated_sub_method_results_df <- kegg_daa_annotated_sub_method_results_df[!is.na(kegg_daa_annotated_sub_method_results_df$pathway_name),]
  kegg_daa_annotated_sub_method_results_df$p_adjust <- round(kegg_daa_annotated_sub_method_results_df$p_adjust,5)
  
  # Selecting the 20 more significant pathways in the errorbar for plotting
  low_p_feature <- kegg_daa_annotated_sub_method_results_df[order(kegg_daa_annotated_sub_method_results_df$p_adjust), ]$feature[1:20]
  
  p <- pathway_errorbar(abundance = kegg_abundance,
                        daa_results_df = kegg_daa_annotated_sub_method_results_df,
                        Group = Group,
                        ko_to_kegg = TRUE,
                        p_values_threshold = 0.05,
                        order = "pathway_class",
                        select = low_p_feature,
                        p_value_bar = TRUE,
                        colors = NULL,
                        x_lab = "pathway_name")
  
  ggsave("KO.png",
         path = Visualizations,
         plot = p, 
         dpi = 320, 
         width = 3,
         height = 2)
  
  a <- pathway_pca(abundance = kegg_abundance,
                   metadata = tibble(metadata),
                   group = col)
  
  ggsave("KO_PCA.png",
         path = Visualizations,
         plot = a, 
         dpi = 320)
  
} else if (method=="MetaCyc"){
  ## Data import
  MetaCyc_abundance <- read.table(MetaCyc, sep = "\t", header = F, row.names = NULL)
  rownames(MetaCyc_abundance) <- MetaCyc_abundance[,1]
  colnames(MetaCyc_abundance) <- MetaCyc_abundance[1,]
  MetaCyc_abundance <- MetaCyc_abundance[-1,-1]
  
  intersect <- intersect(colnames(MetaCyc_abundance),metadata$sample_name)
  metadata <- metadata[match(intersect,metadata$sample_name) ,]
  MetaCyc_abundance <- MetaCyc_abundance[, metadata$sample_name]
  
  for (i in 1:ncol(MetaCyc_abundance)) {
    MetaCyc_abundance[,i] <- as.numeric(MetaCyc_abundance[,i])
  }
  
  group <- col
  Metacyc_daa_df <-
    pathway_daa(
      abundance = MetaCyc_abundance,
      metadata = metadata,
      group = group,
      daa_method = "ALDEx2",
      select = NULL,
      p.adjust = "BH",
      reference = ref_group
    )
  
  if (length(unique(metadata[,col]))==2) {
    Metacyc_daa_sub_method_results_df <-
      Metacyc_daa_df[Metacyc_daa_df$method == "ALDEx2_Wilcoxon rank test", ]
  } else if (length(unique(metadata[,col]))>2) {
  Metacyc_daa_sub_method_results_df <-
    Metacyc_daa_df[Metacyc_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]
  }
  
  Metacyc_daa_annotated_sub_method_results_df <-
    pathway_annotation(pathway = "MetaCyc",
                       daa_results_df = Metacyc_daa_sub_method_results_df,
                       ko_to_kegg = F)
  
  write.table(Metacyc_daa_annotated_sub_method_results_df, 
              file = paste(Visualizations, "/MetaCyc.tsv", sep = ""), 
              sep = "\t", 
              quote = FALSE, 
              row.names = F, 
              col.names = T)
  
  Group <-metadata[,col]
  
  Metacyc_daa_annotated_sub_method_results_df$p_adjust <- round(Metacyc_daa_annotated_sub_method_results_df$p_adjust,5)
  
  low_p_feature <- Metacyc_daa_annotated_sub_method_results_df[order(Metacyc_daa_annotated_sub_method_results_df$p_adjust), ]$feature[1:20]
  
  p <- pathway_errorbar(abundance = MetaCyc_abundance,
                        daa_results_df = Metacyc_daa_annotated_sub_method_results_df,
                        Group = Group,
                        ko_to_kegg = F,
                        p_values_threshold = 0.05,
                        order = "group",
                        select = low_p_feature,
                        p_value_bar = TRUE,
                        colors = NULL,
                        x_lab = "description")
  
  ggsave("MetaCyc.png",
         path = Visualizations,
         plot = p, 
         dpi = 320,
         width = 3,
         height = 2)
  
  a <- pathway_pca(abundance = MetaCyc_abundance,
                   metadata = tibble(metadata),
                   group = col)
  
  ggsave("MetaCyc_PCA.png",
         path = Visualizations,
         plot = a, 
         dpi = 320)
  
} else if (method=="EC"){
  EC_abundance <-
    read.table(EC, sep = "\t", header = F, row.names = NULL)
  rownames(EC_abundance) <- EC_abundance[,1]
  colnames(EC_abundance) <- EC_abundance[1,]
  EC_abundance <- EC_abundance[-1,-1]
  
  intersect <- intersect(colnames(EC_abundance),metadata$sample_name)
  metadata <- metadata[match(intersect,metadata$sample_name) ,]
  EC_abundance <- EC_abundance[, metadata$sample_name]
  
  for (i in 1:ncol(EC_abundance)) {
    EC_abundance[,i] <- as.numeric(EC_abundance[,i])
  }
  
  group <- col
  EC_daa_df <-
    pathway_daa(
      abundance = EC_abundance,
      metadata = metadata,
      group = group,
      daa_method = "ALDEx2",
      select = NULL,
      p.adjust = "BH",
      reference = ref_group
    )
  
  if (length(unique(metadata[,col]))==2) {
    EC_daa_sub_method_results_df <-
      EC_daa_df[EC_daa_df$method == "ALDEx2_Wilcoxon rank test", ]
  } else if (length(unique(metadata[,col]))>2) {
    EC_daa_sub_method_results_df <-
      EC_daa_df[EC_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]
  }
  
  EC_daa_annotated_sub_method_results_df <-
    pathway_annotation(pathway = "EC",
                       daa_results_df = EC_daa_sub_method_results_df,
                       ko_to_kegg = F)
  
  write.table(EC_daa_annotated_sub_method_results_df, 
              file = paste(Visualizations, "/EC.tsv", sep = ""), 
              sep = "\t", 
              quote = FALSE, 
              row.names = F, 
              col.names = T)
  
    EC_daa_annotated_sub_method_results_df$p_adjust <- round(EC_daa_annotated_sub_method_results_df$p_adjust,5)
  
  # Selecting the 20 more significant pathways in the errorbar for plotting
  low_p_feature <- EC_daa_annotated_sub_method_results_df[order(EC_daa_annotated_sub_method_results_df$p_adjust), ]$feature[1:20]
  
  Group <-metadata[,col]
  p <- pathway_errorbar(abundance = EC_abundance,
                   daa_results_df = EC_daa_annotated_sub_method_results_df,
                   Group = Group,
                   ko_to_kegg = F,
                   p_values_threshold = 0.05,
                   order = "group",
                   select = low_p_feature,
                   p_value_bar = TRUE,
                   colors = NULL,
                   x_lab = "description")
  
  ggsave("EC.png",
         path = Visualizations,
         plot = p, 
         dpi = 320,
         width = 3,
         height = 2)
  
  a <- pathway_pca(abundance = EC_abundance,
                   metadata = tibble(metadata),
                   group = col)
  
  ggsave("EC_PCA.png",
         path = Visualizations,
         plot = a, 
         dpi = 320)
}