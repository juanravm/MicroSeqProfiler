#!/bin/bash

#' Import fastq.gz files and quality control (QC)
#'
#' @param input_path Folder name with fastq files
#' @param type Type of the input sequences suitable with QIIME2 importing
#' types
#' @param format Format of the imput sequences suitable with qiime2 importing
#' formats
#' @param trim_left Number of nucleotides removed from the start of the read
#' @param trunc_length Reads truncation length (maximum length of reads)
#' @param cores Number of available cores to the program
#' @param metadata_fp Metadata file path
#' 
#' @return QC-seq.qza QIIME2 artifact with QC unique sequences
#' @return QC-table.qza QIIME2 artifact with the counts of each unique 
#' sequence
#' @return QC-stats.qza QIIME2 artifact with the QC statistics
#' 
#' 

# Variables definition
input_path=""
type=""
format=""
trim_left=""
trunc_length=""
metadata_fp=""
cores=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input_path)
            input_path="$2"
            shift 2
            ;;
        --type)
            type="$2"
            shift 2
            ;;
        --format)
            format="$2"
            shift 2
            ;;
        --trim_left)
            trim_left="$2"
            shift 2
            ;;
        --trunc_length)
            trunc_length="$2"
            shift 2
            ;;
        --metadata_fp)
            metadata_fp="$2"
            shift 2
            ;;
        --cores)
            cores="$2"
            shift 2
            ;;
        *)
            echo "Error in option: $1"
            exit 1
            ;;
    esac
done

echo "Raw sequences folder: $input_path"
echo "QIIME2 type: $type"
echo "QIIME2 format: $format"
echo "$trim_left nt removed left"
echo "Truncation length: $trunc_length nt "
echo "$cores cores employed"


## Import fastq files to QIIME2: In a FastqFiles folder
mkdir intermediate

qiime tools import \
  --type $type \
  --input-path $input_path \
  --input-format $format \
  --output-path ./intermediate/Raw.qza

## Summary of the imported fastq
qiime demux summarize \
  --i-data ./intermediate/Raw.qza \
  --o-visualization ./intermediate/Raw.qzv


## QIIME2 Quality Control:
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./intermediate/Raw.qza \
  --p-trim-left $trim_left \
  --p-trunc-len $trunc_length \
  --p-n-threads $cores \
  --o-representative-sequences QC-seq.qza \
  --o-table QC-table.qza \
  --o-denoising-stats QC-stats.qza

qiime metadata tabulate \
  --m-input-file QC-stats.qza \
  --o-visualization QC-stats.qzv

qiime feature-table tabulate-seqs \
  --i-data QC-seq.qza \
  --o-visualization QC-seq.qzv
  
qiime feature-table summarize \
--i-table QC-table.qza \
--o-visualization QC-table.qzv \
--m-sample-metadata-file $metadata_fp