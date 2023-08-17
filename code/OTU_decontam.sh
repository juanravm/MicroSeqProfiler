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
metadata=""
taxonomy=""

# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --OTU_filtered_table)
            OTU_filtered_table="$2"
            shift 2
            ;;
        --neg)
            neg="$2"
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --metadata)
            metadata="$2"
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
echo "Metadata file path: $metadata"
echo "Negative control identification pattern: $neg"
echo "taxonomy.qza file path: $taxonomy"

## Changing to R scripting
Rscript - <<RSCRIPT

library(decontam)

## Importing OTU table
OTU_table <- as.matrix(t(read.delim("$OTU_filtered_table", row.names=1)))

rownames(OTU_table) <- gsub("^X", "", rownames(OTU_table))
rownames(OTU_table) <- gsub("\\.", "-", rownames(OTU_table))

## Negative logical vector for controls localization
neg <- grepl("$neg", rownames(OTU_table))

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
write.csv(t(noncontam_export),
          row.names = T, 
          col.names = T, 
          file = "$output_dir/decontam_OTU_table.csv")

RSCRIPT

# Modifying .csv to .tsv for importing
sed -e 's/,/\t/g' $output_dir/decontam_OTU_table.csv > $output_dir/decontam_OTU_table.tsv
sed -i '1s/^\t/#OTU ID\t/' $output_dir/intermediate/decontam_OTU_table.tsv
rm $output_dir/decontam_OTU_table.csv

# .tsv to BIOM
biom convert \
-i $output_dir/intermediate/decontam_OTU_table.tsv \
-o $output_dir/intermediate/decontam_OTU_table.biom \
-m $metadata \
--table-type="OTU table" \
--to-hdf5

## Importing decontaminated FeatureTable[Frequency] to QIIME2
qiime tools import \
--input-path $output_dir/decontam_OTU_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path $output_dir/decontam_OTU_table.qza

rm $output_dir/intermediate/decontam_OTU_table.biom

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
-m $metadata \
-o $output_dir/species.txt \
--to-tsv

rm $output_dir/feature-table.biom
mv $output_dir/species.txt $output_dir/species.tsv
