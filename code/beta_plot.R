#!/usr/bin/Rscript
#' Microbiome beta diversity plotting
#'
#' @param metadata_fp metadata.tsv file path
#' 
#' @param bray_fp bray_curtis.tsv file path
#' 
#' @param unweighted_fp unweighted.tsv file path
#' 
#' @param weigthed_fp weighted.tsv file path
#'  
#' @param jaccard_fp jaccard.tsv file path
#' 
#' @param col metadata column with the group of samples to compare
#' 
#' @param output_dir Output directory (Beta Visualizations directory)
#' 
#' @return Visualizations diversity visualizations and comparisons significance
#' by PERMANOVA statistical test

#···················· BETA DIVERSITY STATISTICS
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
require(devtools)
library(pairwiseAdonis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
library(ape)
library(glue)

## Argument vector definition
args <- commandArgs(trailingOnly = TRUE)

## Variables definition
metadata_fp <- NULL
bray_fp <- NULL
unweighted_fp <- NULL
weigthed_fp <-NULL
jaccard_fp <-NULL
col <- NULL
output_dir <- NULL

## Argument assign
for (i in seq_along(args)) {
  if (args[i] == "--metadata_fp") {
    metadata_fp <- args[i + 1]
  } else if (args[i] == "--bray_fp") {
    bray_fp <- args[i + 1]
  } else if (args[i] == "--unweighted_fp") {
    unweighted_fp <- args[i + 1]
  } else if (args[i] == "--weigthed_fp") {
    weighted_fp <- args[i + 1]
  } else if (args[i] == "--jaccard_fp") {
    jaccard_fp <- args[i + 1]
  } else if (args[i] == "--col") {
    col <- args[i + 1]
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
  }
}
output_dir <- gsub("/$", "", output_dir)


## Importing bray_curtis diversity matrix 
# Importing and ordering metadata and beta diversity matrix
metadata <- read.delim(metadata_fp, row.names = 1)
metadata <- metadata[order(rownames(metadata)) ,]

bray_curtis <- read.delim(bray_fp, row.names = 1)
colnames(bray_curtis) <- rownames(bray_curtis)

unweighted <- read.delim(unweighted_fp, row.names = 1)
colnames(unweighted) <- rownames(unweighted)

weighted <- read.delim(weighted_fp, row.names = 1)
colnames(weighted) <- rownames(weighted)

jaccard <- read.delim(jaccard_fp, row.names = 1)
colnames(jaccard) <- rownames(jaccard)

# Match metadata names to beta diversity
metadata <- metadata[(rownames(bray_curtis)) ,]
metadata <- metadata[match(rownames(bray_curtis), rownames(metadata)) ,]

# Filtering controls
bray_curtis <- bray_curtis[!is.na(metadata[,col]),!is.na(metadata[,col])]
unweighted <- unweighted[!is.na(metadata[,col]),!is.na(metadata[,col])]
weighted <- weighted[!is.na(metadata[,col]),!is.na(metadata[,col])]
jaccard <- jaccard[!is.na(metadata[,col]),!is.na(metadata[,col])]

## Comparisons vectors
metadata <- metadata[!is.na(metadata[,col]),]
group <- metadata[,col]

# Dist object creation for pairwise.adonis() analysis
dist_bray <- as.dist(bray_curtis, diag = F, upper = F)
dist_unweighted <- as.dist(unweighted, diag = F, upper = F)
dist_weighted <- as.dist(weighted, diag = F, upper = F)
dist_jaccard <- as.dist(jaccard, diag = F, upper = F)

## Permanova statistical analysis
permanova_bray <- pairwise.adonis(x = dist_bray,
                                     factors = group,
                                     p.adjust.m = "fdr",
                                     reduce = NULL,
                                     perm = 999)

write.csv(permanova_bray, 
          file = paste(output_dir, "/bray_pval.csv", sep = ""), 
          row.names = F, 
          col.names = T)

permanova_unweighted <- pairwise.adonis(x = dist_unweighted,
                                  factors = group,
                                  p.adjust.m = "fdr",
                                  reduce = NULL,
                                  perm = 999)

write.csv(permanova_unweighted, 
          file = paste(output_dir, "/unweighted_pval.csv", sep = ""), 
          row.names = F, 
          col.names = T)

permanova_weighted <- pairwise.adonis(x = dist_weighted,
                                  factors = group,
                                  p.adjust.m = "fdr",
                                  reduce = NULL,
                                  perm = 999)

write.csv(permanova_weighted, 
          file = paste(output_dir, "/weighted_pval.csv", sep = ""), 
          row.names = F, 
          col.names = T)

permanova_jaccard <- pairwise.adonis(x = dist_jaccard,
                                  factors = group,
                                  p.adjust.m = "fdr",
                                  reduce = NULL,
                                  perm = 999)

write.csv(permanova_jaccard, 
          file = paste(output_dir, "/jaccard_pval.csv", sep = ""), 
          row.names = F, 
          col.names = T)

#···················· BETA DIVERSITY PLOTTING

color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
color <- color[1:length(unique(metadata[,col]))]
names(color) <- unique(metadata[,col])

#### BRAY CURTIS
## PCoA statistics calculation
PCoA <- cmdscale(bray_curtis, k = 3, eig = TRUE, add = TRUE)
positions <- PCoA$points
colnames(positions) <- c("pcoa1", "pcoa2", "pcoa3")

positions <- as.data.frame(cbind(positions,
                                 "Group"=metadata[,col][match(rownames(positions), rownames(metadata))]))

positions[,"Group"] <- factor(positions[,"Group"])

## Percentage explained
percent_explained <- format(round((100*PCoA$eig/sum(PCoA$eig)), digits = 1),
                            nsmall = 1, 
                            trim = T)
labs <- c(glue("PCoA 1 ({percent_explained[1]} %)"), 
          glue("PCoA 2 ({percent_explained[2]} %)"))

positions[, 1] <- as.numeric(positions[, 1]) 
positions[, 2] <- as.numeric(positions[, 2]) 
positions[, 3] <- as.numeric(positions[, 3]) 

# Cluster comparison
a <- ggplot(positions, aes(x = pcoa1, y = pcoa2, color = Group))+
  geom_point()+
  labs(x = labs[1],
       y = labs[2],
       title = "Bray Curtis")+
  theme_linedraw() +
  theme(panel.grid = element_blank())+
  scale_color_manual(values = color[levels(positions[,"Group"])])

ggsave("bray_curtis.png",
       path = output_dir,
       plot = a, 
       dpi = 320)

#### UNWEIGHTED
## PCoA statistics calculation
PCoA <- cmdscale(unweighted, k = 3, eig = TRUE, add = TRUE)
positions <- PCoA$points
colnames(positions) <- c("pcoa1", "pcoa2", "pcoa3")

positions <- as.data.frame(cbind(positions,
                                 "Group"=metadata[,col][match(rownames(positions), rownames(metadata))]))

positions[,"Group"] <- factor(positions[,"Group"])

## Percentage explained
percent_explained <- format(round((100*PCoA$eig/sum(PCoA$eig)), digits = 1),
                            nsmall = 1, 
                            trim = T)
labs <- c(glue("PCoA 1 ({percent_explained[1]} %)"), 
          glue("PCoA 2 ({percent_explained[2]} %)"))

positions[, 1] <- as.numeric(positions[, 1]) 
positions[, 2] <- as.numeric(positions[, 2]) 
positions[, 3] <- as.numeric(positions[, 3]) 

# Cluster comparison
b <- ggplot(positions, aes(x = pcoa1, y = pcoa2, color = Group))+
  geom_point()+
  labs(x = labs[1],
       y = labs[2],
       title = "Unweighted UniFrac")+
  theme_linedraw() +
  theme(panel.grid = element_blank())+
  scale_color_manual(values = color[levels(positions[,"Group"])])

ggsave("unweighted.png",
       path = output_dir,
       plot = b, 
       dpi = 320)

#### WEIGHTED
## PCoA statistics calculation
PCoA <- cmdscale(weighted, k = 3, eig = TRUE, add = TRUE)
positions <- PCoA$points
colnames(positions) <- c("pcoa1", "pcoa2", "pcoa3")

positions <- as.data.frame(cbind(positions,
                                 "Group"=metadata[,col][match(rownames(positions), rownames(metadata))]))

positions[,"Group"] <- factor(positions[,"Group"])

## Percentage explained
percent_explained <- format(round((100*PCoA$eig/sum(PCoA$eig)), digits = 1),
                            nsmall = 1, 
                            trim = T)
labs <- c(glue("PCoA 1 ({percent_explained[1]} %)"), 
          glue("PCoA 2 ({percent_explained[2]} %)"))

positions[, 1] <- as.numeric(positions[, 1]) 
positions[, 2] <- as.numeric(positions[, 2]) 
positions[, 3] <- as.numeric(positions[, 3]) 

# Cluster comparison
c <- ggplot(positions, aes(x = pcoa1, y = pcoa2, color = Group))+
  geom_point()+
  labs(x = labs[1],
       y = labs[2],
       title = "Weighted UniFrac")+
  theme_linedraw() +
  theme(panel.grid = element_blank())+
  scale_color_manual(values = color[levels(positions[,"Group"])])

ggsave("weighted.png",
       path = output_dir,
       plot = c, 
       dpi = 320)

#### WEIGHTED
## PCoA statistics calculation
PCoA <- cmdscale(jaccard, k = 3, eig = TRUE, add = TRUE)
positions <- PCoA$points
colnames(positions) <- c("pcoa1", "pcoa2", "pcoa3")

positions <- as.data.frame(cbind(positions,
                                 "Group"=metadata[,col][match(rownames(positions), rownames(metadata))]))

positions[,"Group"] <- factor(positions[,"Group"])

## Percentage explained
percent_explained <- format(round((100*PCoA$eig/sum(PCoA$eig)), digits = 1),
                            nsmall = 1, 
                            trim = T)
labs <- c(glue("PCoA 1 ({percent_explained[1]} %)"), 
          glue("PCoA 2 ({percent_explained[2]} %)"))

positions[, 1] <- as.numeric(positions[, 1]) 
positions[, 2] <- as.numeric(positions[, 2]) 
positions[, 3] <- as.numeric(positions[, 3]) 

# Cluster comparison
d <- ggplot(positions, aes(x = pcoa1, y = pcoa2, color = Group))+
  geom_point()+
  labs(x = labs[1],
       y = labs[2],
       title = "Jaccard Index")+
  theme_linedraw() +
  theme(panel.grid = element_blank())+
  scale_color_manual(values = color[levels(positions[,"Group"])])

ggsave("jaccard.png",
       path = output_dir,
       plot = d, 
       dpi = 320)
