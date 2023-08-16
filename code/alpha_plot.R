#!/usr/bin/Rscript
#' Microbiome alpha diversity plotting
#'
#' @param shannon_fp shannon.tsv file path
#' 
#' @param evenness_fp evenness_pd.tsv file path
#'  
#' @param faith_fp faith_pd.tsv file path
#' 
#' @param metadata_fp metadata.tsv file pach
#' 
#' @param col metadata column with the group of samples to compare
#' 
#' @param output_dir Output directory (Visualizations directory)
#' 
#' @return alpha diversity visualizations and comparisons significance
#' by Kruskal-Wallis statistical test

#···················· ALPHA DIVERSITY PLOTTING
library(ggplot2)
library(dplyr)

## Argument vector definition
args <- commandArgs(trailingOnly = TRUE)

## Variables definition
shannon_fp <- NULL
evenness_fp <- NULL
faith_fp <- NULL
metadata_fp <-NULL
col <- NULL
output_dir <- NULL

## Argument assign
for (i in seq_along(args)) {
  if (args[i] == "--shannon_fp") {
    shannon_fp <- args[i + 1]
  } else if (args[i] == "--evenness_fp") {
    evenness_fp <- args[i + 1]
  } else if (args[i] == "--faith_fp") {
    faith_fp <- args[i + 1]
  } else if (args[i] == "--metadata") {
    metadata_fp <- args[i + 1]
  } else if (args[i] == "--col") {
    col <- args[i + 1]
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
  }
}

output_dir <- gsub("/$", "", output_dir)

## Metadata input
metadata <- read.delim(metadata_fp)

# Faith PD input
faith.pd <- read.delim(faith_fp)
faith.pd <- cbind(faith.pd, "Group"=metadata[, col][match(faith.pd$X, metadata$Sample.id)])
faith.pd[,Group] <- factor(faith.pd[,Group])

# Pielou evenness input
evenness.pd <- read.delim(evenness_fp)
evenness.pd <- cbind(evenness.pd, Group=metadata[, Group][match(evenness.pd$X, metadata$Sample.id)])
evenness.pd[,Group] <- factor(evenness.pd[,Group])

# Shannon input
shannon <- read.delim(shannon_fp)
shannon <- cbind(shannon, Group=metadata[, Group][match(shannon$X, metadata$Sample.id)])
shannon[,Group] <- factor(shannon[,Group])

#### Kruskal-Wallis statistical test (pairwise)
## Comparisons dataframe creation for each metric
# Faith PD:
p_values1 <- c()
for (i in 1:(length(unique(faith.pd[,Group]))-1)) {
  for (j in (i+1):length(unique(faith.pd[,Group]))) {
    kruskal1 <- kruskal.test(faith_pd ~ Group, 
                             data = faith.pd, 
                             subset = faith.pd$Cluster %in% c(levels(faith.pd[,Group])[i],levels(faith.pd[,Group])[j]))
    p_values1 <- c(p_values1, kruskal1$p.value)}}

group1 <- c()
group2 <- c()

for (i in 1:(length(unique(faith.pd[,Group]))-1)) {
  for (j in (i+1):length(unique(faith.pd[,Group]))) {
    group1 <- c(group1, i)
    group2 <- c(group2, j)}}

faith_pval <- data.frame(group1=levels(faith.pd[,Group])[group1],
                            group2=levels(faith.pd[,Group])[group2],
                            p_values1=p_values1)

write.csv(faith_pval, 
          file = paste(visualization_fp, "/faith_pval.csv"), 
          row.names = T, 
          col.names = T)

## Pielou evenness
p_values1 <- c()
for (i in 1:(length(unique(evenness.pd[,Group]))-1)) {
  for (j in (i+1):length(unique(evenness.pd[,Group]))) {
    kruskal1 <- kruskal.test(pielou_evenness ~ Group, 
                             data = evenness.pd, 
                             subset = evenness.pd$Cluster %in% c(levels(evenness.pd[,Group])[i],levels(evenness.pd[,Group])[j]))
    p_values1 <- c(p_values1, kruskal1$p.value)}}

group1 <- c()
group2 <- c()

for (i in 1:(length(unique(evenness.pd[,Group]))-1)) {
  for (j in (i+1):length(unique(evenness.pd[,Group]))) {
    group1 <- c(group1, i)
    group2 <- c(group2, j)}}

evenness_pval <- data.frame(group1=levels(evenness.pd[,Group])[group1],
                         group2=levels(evenness.pd[,Group])[group2],
                         p_values1=p_values1)

write.csv(evenness_pval, 
          file = paste(visualization_fp, "/evenness_pval.csv"), 
          row.names = T, 
          col.names = T)

# Shannon:
p_values1 <- c()
for (i in 1:(length(unique(shannon[,Group]))-1)) {
  for (j in (i+1):length(unique(shannon[,Group]))) {
    kruskal1 <- kruskal.test(shannon_entropy ~ Group, 
                             data = shannon, 
                             subset = shannon$Cluster %in% c(levels(shannon[,Group])[i],levels(shannon[,Group])[j]))
    p_values1 <- c(p_values1, kruskal1$p.value)}}

group1 <- c()
group2 <- c()

for (i in 1:(length(unique(shannon[,Group]))-1)) {
  for (j in (i+1):length(unique(shannon[,Group]))) {
    group1 <- c(group1, i)
    group2 <- c(group2, j)}}

shannon_pval <- data.frame(group1=levels(shannon[,Group])[group1],
                            group2=levels(shannon[,Group])[group2],
                            p_values1=p_values1)

write.csv(shannon_pval, 
          file = paste(visualization_fp, "/shannon_pval.csv"), 
          row.names = T, 
          col.names = T)

## Alpha diversity boxplots visualizations
color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
color <- color[1:length(unique(metadata[,col]))]
names(color) <- unique(metadata[,col])

#### Evenness boxplot
# - Evenness Cluster comparison
a <-ggplot(evenness.pd, aes(x=Group, y=pielou_evenness)) + 
  geom_boxplot(fill=color[levels(evenness.pd[,Group])], color="black", width = 0.6) + 
  theme_linedraw() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) +
  labs(title = "Pielou evenness", 
       x = levels(evenness.pd[,Group]), 
       y = "Pielou evenness") +
  theme(axis.title = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("pielou_evenness.png",
       path = output_dir,
       plot = a, 
       dpi = 320)


#### Faith PD boxplot
# - Faith PD Cluster Comparison
b <- ggplot(faith.pd, aes(x=Cluster, y=faith_pd)) + 
  geom_boxplot(fill=color[levels(faith.pd[,Group])], color="black", width = 0.6) + 
  theme_linedraw() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) +
  labs(title = "Faith PD", 
       x = levels(faith.pd[,Group]), 
       y = "Faith PD") +
  theme(axis.title = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("faith_pd.png",
       path = output_dir,
       plot = b, 
       dpi = 320)

#### Shannon boxplot
# - Shannon Cluster Comparison
c <- ggplot(shannon, aes(x=Cluster, y=shannon_entropy)) + 
  geom_boxplot(fill=color[levels(shannon[,Group])], color="black", width = 0.6) + 
  theme_linedraw() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) +
  labs(title = "Shannon index", 
       x = levels(shannon[,Group]), 
       y = "Shannon Entropy") +
  theme(axis.title = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("shannon.png", 
       path = output_dir,
       plot = c,
       dpi = 320)
