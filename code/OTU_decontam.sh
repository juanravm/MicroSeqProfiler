#!/bin/bash
#' OTU decontamination (Negative control normalization)
#'
#' @param OTU_filtered_table OTU filtered table file path with 
#' nonchimeric counts
#' 
#' @param neg A character pattern to identify negative controls in
#' OTU table rownames
#' 
#' @param  output_dir Decontaminated OTU table output directory
#' (usually the previous microbiome working directory employed)
#' 
#' @param metadata Metadata file path
#' 
#' @param taxonomy taxonomy.qza file path
#' 
#' @return decontam_OTU_table.qza QIIME2 artifact with the chimera-filtered 
#' and decontaminated OTU counts
#'  
#' @return species.tsv File with the chimera-filtered and decontaminated 
#' OTU counts in .tsv format

# Variables definition
OTU_filtered_table=""
neg=""
output_dir=""
metadata_fp=""
taxonomy=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --OTU_filtered_table)
            OTU_filtered_table="$2"
            shift 2
            ;;
        --neg)
            neg="\"$2\""
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --metadata_fp)
            metadata_fp="$2"
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

echo "OTU table file path: $OTU_filtered_table"
echo "Output directory: $output_dir"
echo "Metadata file path: $metadata_fp"
echo "Negative control identification pattern: $neg"
echo "taxonomy.qza file path: $taxonomy"

## Changing to R scripting
Rscript - <<RSCRIPT

library(decontam)

## Importing OTU table
OTU_table <- t(read.delim("$OTU_filtered_table", 
header=FALSE, 
comment.char=""))

rownames(OTU_table) <- OTU_table[,1]
OTU_table <- OTU_table[,-1]
colnames(OTU_table) <- OTU_table[1,]
OTU_table <- OTU_table[-1,]

rownames <- rownames(OTU_table)
OTU_table <- apply(OTU_table, 2, function(x) as.numeric(x))
rownames(OTU_table) <- rownames


## Negative logical vector for controls localization
neg <- grepl($neg, rownames(OTU_table))

## Logical vector indicating contaminants
Contaminant <- as.logical(isContaminant(
  seqtab = OTU_table,
  conc = NULL,
  neg = neg,
  method = "auto",
  batch = NULL,
  threshold = 0.1,
  normalize = FALSE,
  batch.combine = "minimum",
  detailed = F
))

## Removing frequencies from contaminated microorganisms
OTU_table[,Contaminant]=0

## Exporting decontaminated table without experimental controls
noncontam_export <- OTU_table[!neg ,]
file <- paste("$output_dir", "/decontam_OTU_table.csv", sep = "")
write.csv(t(noncontam_export),
          row.names = T, 
          col.names = T, 
          file = file)

RSCRIPT

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