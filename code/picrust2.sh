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
OTU_table_fp=""
OTU_seq_fp=""
cores=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --OTU_table_fp)
            OTU_table_fp="$2"
            shift 2
            ;;
        --OTU_seq_fp)
            OTU_seq_fp="$2"
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

## PICRUSt2 code download
wget https://github.com/picrust/picrust2/archive/v2.5.2.tar.gz
tar xvzf  v2.5.2.tar.gz

rm v2.5.2.tar.gz
mv picrust2-2.5.2/ picrust2
cd picrust2/

#QIIME2 exporting sequences and feature-table :
# input for Picrust2
qiime tools export \
--input-path $OTU_seq_fp \
--output-path ./input
fasta="'$(realpath input/dna-sequences.fasta)'"

qiime tools export \
--input-path $OTU_table_fp \
--output-path ./input
biom="'$(realpath input/feature-table.biom)'"

## Enviroment creation and picrust installation with pip
conda env create -n picrust2 -f picrust2-env.yaml
conda activate picrust2
pip install --editable .

## Run picrust2_pipeline.py downloaded in scripts folder
cd scripts
python picrust2_pipeline.py -s  $fasta \
-i $biom \
-o ../output \
-t sepp \
-p $cores \
-e 0 \
--verbose

cd ..

## Unzip compressed outputs
# Recursively all .gz files in the picrust2 output
gunzip -r *.gz ./output
