#!/bin/bash
#' Microbiome alpha and beta diversity metrics analysis
#'
#' @param OTU_filtered_seqs OTU_rep_seqs.qza QIIME2 artifact with the 
#' representative chimera-filtered sequence of each OTU
#' 
#' @param n_trees Number of phylogenetic trees to be constructed
#' for diversity analysis
#'  
#' @param cores Number of processor cores to use during the analysis
#' 
#' @param OTU_table OTU_table.qza QIIME2 artifact with the OTU counts
#' 
#' @param metadata Metadata file path
#' 
#' @param controls Metadata column character indicating experimental
#' controls (e.g. "[Coltrols]!='Yes'", removing samples "yes" from 
#' metadata column "Controls", retaining non controls)
#' 
#' @param sampling Minimum sampling depth to calculate diversity metrics
#' 
#' @return diversity A diversity directory divided in other 2 
#' subdirectories with alpha and beta diversity metrics

# Variables definition
OTU_filtered_seqs=""
n_trees=""
cores=""
OTU_table=""
metadata=""
controls=""
sampling=""

## Arguments assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --OTU_filtered_seqs)
            OTU_filtered_seqs="$2"
            shift 2
            ;;
        --n_trees)
            n_trees="$2"
            shift 2
            ;;
        --cores)
            cores="$2"
            shift 2
            ;;
        --OTU_table)
            OTU_table="$2"
            shift 2
            ;;
        --metadata)
            metadata="$2"
            shift 2
            ;;
        --controls)
            controls="$2"
            shift 2
            ;;
        --sampling)
            sampling="$2"
            shift 2
            ;;
        *)
            echo "Error in option: $1"
            exit 1
            ;;
    esac
done


#···················· PHYLOGENETIC TREE CONSTRUCTION
## Reads alignment
qiime alignment mafft \
  --i-sequences $OTU_filtered_seqs \
  --o-alignment intermediate/aligned_rep_seqs.qza

## Alignment mask
qiime alignment mask \
  --i-alignment intermediate/aligned-rep-seqs.qza \
  --o-masked-alignment intermediate/masked_aligned_rep_seqs.qza

## Phylogenetic tree construction
qiime phylogeny raxml \
  --i-alignment masked_aligned_rep_seqs.qza \
  --p-substitution-model GTRCAT \
  --p-raxml-version SSE3 \
  --p-seed 1723 \
  --p-n-searches $n_trees \
  --o-tree intermediate/unrooted_tree.qza \
  --p-n-threads $cores \
  --verbose

## Setting a midpoint root
qiime phylogeny midpoint-root \
   --i-tree intermediate/unrooted_tree.qza \
   --o-rooted-tree intermediate/rooted_tree.qza

## Filter OTU_table to aligned sequences
qiime phylogeny filter-table \
   --i-table $OTU_table \
   --i-tree intermediate/rooted_tree.qza \
   --o-filtered-table intermediate/aligned_OTU_table.qza

## Filter OTU_table to remove experimental controls
qiime feature-table filter-samples \
  --i-table intermediate/aligned_OTU_table.qza \
  --m-metadata-file $metadata \
  --p-where "$controls" \
  --o-filtered-table intermediate/aligned_OTU_table.qza

#···················· DIVERSITY METRICS CALCULATION
## Generating diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny intermediate/rooted_tree.qza \
  --i-table intermediate/aligned_OTU_table.qza \
  --p-sampling-depth $sampling \
  --m-metadata-file $metadata \
  --output-dir diversity \
  --p-n-jobs-or-threads $cores

# Moving diversity metrics
cd diversity
mkdir alpha
mkdir beta
mv  *faith* *evenness* *shannon* ./alpha
mv *jaccard* *bray* *unweighted* *weighted* ./beta

## Exporting alpha diversity metrics to .tsv
cd alpha

qiime tools export \
--input-path faith_pd_vector.qza \
--output-path .

mv alpha-diversity.tsv faith_pd.tsv 

qiime tools export \
--input-path evenness_vector.qza \
--output-path .

mv alpha-diversity.tsv evenness_pd.tsv 

qiime tools export \
--input-path shannon_vector.qza \
--output-path .

mv alpha-diversity.tsv shannon.tsv 

mkdir Visualizations

## Exporting beta diversity metrics to .tsv
cd ../beta/
  
qiime tools export \
--input-path bray_curtis_distance_matrix.qza \
--output-path .

mv distance-matrix.tsv bray_curtis.tsv 

qiime tools export \
--input-path jaccard_distance_matrix.qza \
--output-path .

mv distance-matrix.tsv jaccard.tsv 

qiime tools export \
--input-path unweighted_unifrac_distance_matrix.qza \
--output-path .

mv distance-matrix.tsv unweighted.tsv 

qiime tools export \
--input-path weighted_unifrac_distance_matrix.qza \
--output-path .

mv distance-matrix.tsv weighted.tsv 

mkdir Visualizations
