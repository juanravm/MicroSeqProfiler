#!/bin/bash
#' OTU decontamination (Negative control normalization)
#'
#' @param OTU_table_fp decontam_OTU_table.qza QIIME2 artifact file path 
#' with decontaminated OTU counts
#' 
#' @param OTU_seq_fp OTU_filtered_seqs.qza QIIME2 artifact file path 
#' with chimera filtered sequences
#' 
#' @param  cores Number of processor cores available to the program
#' 
#' @return output Directory with the inferred metabolic pathways matrix to
#' KEGG Orthologs, EC and MetaCyc for each sample

# Variables definition
ref_group=""
col=""
metadata_fp=""
method=""
picrust_dir=""


KO="\"$(realpath $picrust_dir/KO_metagenome_out/pred_metagenome_unstrat.tsv)\""
MetaCyc="\"$(realpath $picrust_dir/pathways_out/path_abun_unstrat.tsv)\""
EC="\"$(realpath $picrust_dir/EC_metagenome_out/pred_metagenome_unstrat.tsv)\""
mkdir $picrust_dir/Visualizations
Visualizations="\"$(realpath $picrust_dir/Visualizations)"
##······················ DIFFERENTIAL ABUNDANCE ANALYSIS
Rscript - <<RSCRIPT
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

method <- "$method"
metadata <- as.data.frame(read_delim("$metadata_fp", 
                                     comment="#", 
                                     escape_double = FALSE,
                                     trim_ws = TRUE))

colnames(metadata)[1] <- "sample_name"
path <- "$Visualizations"

if (method=="KO"){
  kegg_abundance <- ko2kegg_abundance("$KO")
  kegg_abundance <- kegg_abundance[, metadata$sample_name]
  
group <- "$col"
kegg_daa_df <-
    pathway_daa(
      abundance = kegg_abundance,
      metadata = metadata,
      group = group,
      daa_method = "ALDEx2",
      select = NULL,
      p.adjust = "BH",
      reference = "$ref_group"
    )
  
kegg_daa_sub_method_results_df <-
    kegg_daa_df[kegg_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]
  
kegg_daa_annotated_sub_method_results_df <-
    pathway_annotation(pathway = "KO",
                       daa_results_df = kegg_daa_sub_method_results_df,
                       ko_to_kegg = TRUE)

write.table(kegg_daa_annotated_sub_method_results_df, 
              file = paste(path, "/KO.tsv"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = F, 
              col.names = T)
  
Group <- Group <-metadata [,"$col"]
  
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
       path = path,
       plot = p, 
       dpi = 320)
       
} else if (method=="MetaCyc"){
## Data import
MetaCyc_abundance <- read.table("$MetaCyc", sep = "\t", header = F, row.names = NULL)
rownames(MetaCyc_abundance) <- MetaCyc_abundance[,1]
colnames(MetaCyc_abundance) <- MetaCyc_abundance[1,]
MetaCyc_abundance <- MetaCyc_abundance[-1,-1]

MetaCyc_abundance <- MetaCyc_abundance[, metadata$sample_name]

group <- "$col"
Metacyc_daa_df <-
  pathway_daa(
    abundance = MetaCyc_abundance,
    metadata = metadata,
    group = group,
    daa_method = "ALDEx2",
    select = NULL,
    p.adjust = "BH",
    reference = "$ref_group"
  )

Metacyc_daa_sub_method_results_df <-
  Metacyc_daa_df[Metacyc_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]
  
Metacyc_daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "MetaCyc",
                     daa_results_df = Metacyc_daa_sub_method_results_df,
                     ko_to_kegg = F)
                     
write.table(Metacyc_daa_annotated_sub_method_results_df, 
              file = paste(path, "/MetaCyc.tsv"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = F, 
              col.names = T)
                     
Group <-metadata[,"$col"]

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
       path = path,
       plot = p, 
       dpi = 320)
       
} else if (method=="EC"){
EC_abundance <-
  read.table("$EC", sep = "\t", header = F, row.names = NULL)
rownames(EC_abundance) <- EC_abundance[,1]
colnames(EC_abundance) <- EC_abundance[1,]
EC_abundance <- EC_abundance[-1,-1]

EC_abundance <- EC_abundance[, metadata$sample_name]

group <- "$col"
EC_daa_df <-
  pathway_daa(
    abundance = EC_abundance,
    metadata = metadata,
    group = group,
    daa_method = "ALDEx2",
    select = NULL,
    p.adjust = "BH",
    reference = "$ref_group"
  )

EC_daa_sub_method_results_df <-
  EC_daa_df[EC_daa_df$method == "ALDEx2_Kruskal-Wallace test", ]

EC_daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "EC",
                     daa_results_df = EC_daa_sub_method_results_df,
                     ko_to_kegg = F)

EC_daa_annotated_sub_method_results_df$p_adjust <- round(EC_daa_annotated_sub_method_results_df$p_adjust,5)

# Selecting the 20 more significant pathways in the errorbar for plotting
low_p_feature <- EC_daa_annotated_sub_method_results_df[order(EC_daa_annotated_sub_method_results_df$p_adjust), ]$feature[1:20]

Group <-metadata[,"$col"]
pathway_errorbar(abundance = EC_abundance_daa,
                 daa_results_df = EC_daa_annotated_sub_method_results_df,
                 Group = Group,
                 ko_to_kegg = F,
                 p_values_threshold = 0.05,
                 order = "group",
                 select = low_p_feature,
                 p_value_bar = TRUE,
                 colors = NULL,
                 x_lab = "description")
}