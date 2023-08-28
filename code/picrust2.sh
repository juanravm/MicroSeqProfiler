#!/bin/bash
#' OTU decontamination (Negative control normalization)
#'
#' @param biom feature-table.biom with decontam OTU table 
#' located in de picrust/input directory
#' 
#' @param fasta dna-sequences.fasta with OTU filtered sequences
#' 
#' @param picrust2_fp picrust2_pipeline.py file path downloaded 
#' in scripts folder from PICRUSt2 GitHub
#' 
#' @param  cores Number of processor cores available to the program
#' 
#' @return output Directory with the inferred metabolic pathways matrix to
#' KEGG Orthologs, EC and MetaCyc for each sample


# Variables definition
cores=""
fasta=""
biom=""
output_dir=""
picrust2_fp=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores)
            cores="$2"
            shift 2
            ;;
        --fasta)
            fasta="$2"
            shift 2
            ;;
        --biom)
            biom="$2"
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --picrust2_fp)
            picrust2_fp="$2"
            shift 2
            ;;
        *)
            echo "Error in option: $1"
            exit 1
            ;;
    esac
done


#QIIME2 exporting sequences and feature-table :
# input for Picrust2
python $picrust2_fp -s  $fasta \
-i $biom \
-o $output_dir \
-t sepp \
-p $cores \
-e 0 \
--verbose


## Unzip compressed outputs
# Recursively all .gz files in the picrust2 output
gunzip -r *.gz $output_dir