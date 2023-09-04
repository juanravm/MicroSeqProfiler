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
#' @param decontam_OTU_table decontam_OTU_table.csv file path obtained 
#' from OTU_decontam.R
#' 
#' @param taxonomy taxonomy.qza QIIME2 artifact file path
#' 
#' @param metadata_fp Metadata file path
#' 
#' @param column Metadata column to indicate experimental controls 
#' 
#' @param pattern Character string to identify controls in metadata column
#' 
#' @param sampling Minimum sampling depth to calculate diversity metrics
#' 
#' @return diversity A diversity directory divided in other 2 
#' subdirectories with alpha and beta diversity metrics

# Variables definition
OTU_filtered_seqs=""
n_trees=""
cores=""
decontam_OTU_table=""
taxonomy=""
metadata_fp=""
column=""
pattern=""
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
        --decontam_OTU_table)
            decontam_OTU_table="$2"
            shift 2
            ;;
        --metadata_fp)
            metadata_fp="$2"
            shift 2
            ;;
        --column)
            column="$2"
            shift 2
            ;;
        --pattern)
            pattern="$2"
            shift 2
            ;;
        --sampling)
            sampling="$2"
            shift 2
            ;;
        --taxonomy)
            taxonomy="$2"
            shift 2
            ;;
        *)
            echo "Error in option: $1"
            exit 1
            ;;
    esac
done

output_dir="${decontam_OTU_table%/*}"

# Modifying .csv to .tsv for importing
sed -e 's/,/\t/g' $output_dir/decontam_OTU_table.csv > $output_dir/decontam_OTU_table.tsv
sed -i '1s/^""\t/"#OTU ID"\t/' $output_dir/decontam_OTU_table.tsv
sed -i 's/"//g' $output_dir/decontam_OTU_table.tsv

rm $output_dir/decontam_OTU_table.csv

# .tsv to BIOM
biom convert \
-i $output_dir/decontam_OTU_table.tsv \
-o $output_dir/decontam_OTU_table.biom \
-m $metadata_fp \
--table-type="OTU table" \
--to-hdf5

## Importing decontaminated FeatureTable[Frequency] to QIIME2
qiime tools import \
--input-path $output_dir/decontam_OTU_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path $output_dir/decontam_OTU_table.qza

OTU_table="$(realpath $output_dir/decontam_OTU_table.qza)"
mv $output_dir/decontam_OTU_table.biom ./picrust2/input/
mv ./picrust2/input/decontam_OTU_table.biom ./picrust2/input/feature-table.biom


# Species taxa collapse
qiime taxa collapse \
    --i-table $output_dir/decontam_OTU_table.qza \
    --i-taxonomy $taxonomy \
    --p-level 7 \
    --o-collapsed-table $output_dir/species.qza

qiime composition add-pseudocount \
  --i-table $output_dir/species.qza \
  --o-composition-table $output_dir/species.qza
  
qiime tools export \
   --input-path $output_dir/species.qza \
   --output-path $output_dir


# Convert biom format to tsv
biom convert \
-i $output_dir/feature-table.biom \
-m $metadata_fp \
-o $output_dir/species.txt \
--to-tsv

rm $output_dir/feature-table.biom
mv $output_dir/species.txt $output_dir/species.tsv

sed -i '1d' $output_dir/species.tsv
sed -i '1s/^#OTU ID/id/' $output_dir/species.tsv
mkdir $output_dir/LEfSe

#···················· PHYLOGENETIC TREE CONSTRUCTION
## Reads alignment
qiime alignment mafft \
  --i-sequences $OTU_filtered_seqs \
  --o-alignment ./intermediate/aligned_rep_seqs.qza

## Alignment mask
qiime alignment mask \
  --i-alignment ./intermediate/aligned_rep_seqs.qza \
  --o-masked-alignment ./intermediate/masked_aligned_rep_seqs.qza

## Phylogenetic tree construction
qiime phylogeny raxml \
  --i-alignment ./intermediate/masked_aligned_rep_seqs.qza \
  --p-substitution-model GTRCAT \
  --p-raxml-version SSE3 \
  --p-seed 1723 \
  --p-n-searches $n_trees \
  --o-tree ./intermediate/unrooted_tree.qza \
  --p-n-threads $cores \
  --verbose

## Setting a midpoint root
qiime phylogeny midpoint-root \
   --i-tree ./intermediate/unrooted_tree.qza \
   --o-rooted-tree ./intermediate/rooted_tree.qza

## Filter decontam_OTU_table to aligned sequences
qiime phylogeny filter-table \
   --i-table $OTU_table \
   --i-tree ./intermediate/rooted_tree.qza \
   --o-filtered-table ./intermediate/aligned_OTU_table.qza

## Filter decontam_OTU_table to remove experimental controls
qiime feature-table filter-samples \
  --i-table ./intermediate/aligned_OTU_table.qza \
  --m-metadata-file $metadata_fp \
  --p-where "[$column]!='$pattern'" \
  --o-filtered-table ./intermediate/aligned_OTU_table.qza

#···················· DIVERSITY METRICS CALCULATION
## Generating diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ./intermediate/rooted_tree.qza \
  --i-table ./intermediate/aligned_OTU_table.qza \
  --p-sampling-depth $sampling \
  --m-metadata-file $metadata_fp \
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
