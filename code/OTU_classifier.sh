#!/bin/bash

#' Import fastq.gz files and quality control (QC)
#'
#' @param QC_table QIIME2 FeatureTable[Frequency] artifact with 
#' QC read counts
#' 
#' #' @param QC_seqs QIIME2 FeatureData[Sequence] artifact with QC rep-seqs
#' 
#' @param perc_identity Identity percentage to cluster OTUs: From 0 to 1
#' 
#' @param f_primer Forward primer sequence in 5' -> 3' 
#' 
#' @param r_primer Reverse primer sequence in 5' -> 3' 
#' 
#' @param class_seq Classifier reference sequences for microorganisms, usually
#' from GreenGenes or Silva databases
#' 
#' @param class_tax Taxonomy associated to classifier references
#' 
#' @param metadata Metadata file path
#' 
#' @return taxa-barplot QIIME2 visualization with taxonomic composition of
#' samples
#' 
#' @return OTU_filtered_seqs.qza QIIME2 artifact with the representative
#' chimera-filtered sequence of each OTU
#' 
#' @return OTU_filtered_table.qza QIIME2 artifact with the OTU filtered counts
#' 
#' 

# Variables definition
QC_table=""
QC_seqs=""
perc_identity=""
f_primer=""
r_primer=""
class_seq=""
class_tax=""
metadata=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --QC_table)
            QC_table="$2"
            shift 2
            ;;
        --QC_seqs)
            QC_seqs="$2"
            shift 2
            ;;
        --perc_identity)
            perc_identity="$2"
            shift 2
            ;;
        --f_primer)
            f_primer="$2"
            shift 2
            ;;
        --r_primer)
            r_primer="$2"
            shift 2
            ;;
        --class_seq)
            class_seq="$2"
            shift 2
            ;;
        --class_tax)
            class_tax="$2"
            shift 2
            ;;
        --metadata)
            metadata="$2"
            shift 2
            ;;
        *)
            echo "Opción inválida: $1"
            exit 1
            ;;
    esac
done

echo "QC_table file path: $QC_table"
echo "QC_seqs file path: $QC_seqs"
echo "Identity for clustering: $perc_identity"
echo "Forward primer sequence: $f_primer"
echo "Reverse primer sequence: $r_primer"
echo "Reference classifier seqs file path: $class_seq"
echo "Reference classifier taxa file path: $class_taxa"
echo "Metadata file path: $class_seq"

#···················· CLUSTERING OTUS
## 99% Identity clustering sequences: 
qiime vsearch cluster-features-de-novo \
--i-table $QC_table \
--i-sequences $QC_seqs \
--p-perc-identity $perc_identity \
--o-clustered-table intermediate/OTU_table.qza \
--o-clustered-sequences intermediate/OTU_seqs.qza

## Chimera filtering 
qiime vsearch uchime-denovo \
--i-table intermediate/OTU_table.qza \
--i-sequences intermediate/OTU_seqs.qza \
--output-dir intermediate/chimera_filter

qiime feature-table filter-features \
  --i-table intermediate/OTU_table.qza \
  --m-metadata-file intermediate/chimera_filter/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-table OTU_filtered_table.qza
  
qiime feature-table filter-seqs \
  --i-data intermediate/OTU_seqs.qza \
  --m-metadata-file intermediate/chimera_filter/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-data OTU_filtered_seqs.qza

#···················· TAXONOMIC CLASSIFICATION

# Sequences import
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $class_seq \
--output-path intermediate/ref-seq.qza

# Taxonomy import
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $class_tax \
--output-path intermediate/ref-taxonomy.qza

## Extracting classifier reads for our primers
qiime feature-classifier extract-reads \
--i-sequences intermediate/ref-seq.qza \
--p-f-primer $f_primer \
--p-r-primer $r_primer \
--o-reads intermediate/classifier-seqs.qza

## Creating classifier from reads and taxonomy
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads intermediate/classifier-seqs.qza \
--i-reference-taxonomy intermediate/ref-taxonomy.qza \
--o-classifier intermediate/classifier.qza

## Taxonomic classification of our reads
qiime feature-classifier classify-sklearn \
--i-classifier intermediate/classifier.qza \
--i-reads OTU_filtered_seqs.qza \
--o-classification intermediate/taxonomy.qza

qiime metadata tabulate \
--m-input-file intermediate/taxonomy.qza \
--o-visualization intermediate/taxonomy.qzv

# Taxonomic visualization
qiime taxa barplot \
--i-table OTU_filtered_table.qza \
--i-taxonomy intermediate/taxonomy.qza \
--m-metadata-file $metadata \
--o-visualization taxa-barplots.qzv
