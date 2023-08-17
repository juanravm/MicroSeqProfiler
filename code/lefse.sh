#!/bin/bash
#' LEfSe analysis
#'
#' @param species_fp species.qza QIIME2 artifact file path
#' 
#' @param neg A character pattern to identify negative controls in
#' OTU table rownames
#' 
#' @param  output_dir Decontaminated OTU table output directory
#' (usually the previous microbiome working directory employed)
#' 
#' @param metadata_fp Metadata file path
#' 
#' @param taxonomy taxonomy.qza file path
#' 
#' @return decontam_OTU_table.qza QIIME2 artifact with the chimera-filtered 
#' and decontaminated OTU counts
#'  
#' @return species.tsv File with the chimera-filtered and decontaminated 
#' OTU counts in .tsv format

# Variables definition
species_fp=""
neg=""
output_dir=""
metadata_fp=""
taxonomy=""
sampling=""
col=""
ref_group=""
group2=""


# Argument assign to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --species_fp)
            species_fp="$2"
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

# Exporting species composition
qiime tools export \
   --input-path $species_fp \
   --output-path LEfSe

biom convert \
-i LEfSe/feature-table.biom \
-m $metadata_fp \
-o LEfSe/species.txt \
--to-tsv

rm feature-table.biom
mv LEfSe/species.txt LEfSe/species.tsv

# First row elimination for .tsv files 
# (elimination of #Constructed from BIOM file) 
sed -i '1d' LEfSe/species.tsv

## We must eliminate "#OTU ID" from 1st cell for LEfSe input
# First cell substitution to "id"
sed -i '1s/^[^\t]*/id/' LEfSe/species.tsv
species=realpath LEfSe/species.tsv
LEfSe=realpath LEfSe/
Rscript - <<RSCRIPT
library(lefser)

## Importing data
comparison <- read.delim("$species", header=T, row.names = 1)
colnames(comparison) <- gsub("\\.", "_", colnames(comparison))
colnames(comparison) <- gsub("^X", "", colnames(comparison))

metadata <- read.delim($metadata_fp, comment.char="#")

# - Removing samples with low sampling depth such as in diversity analysis
sampling_depth <- colSums(comparison)
sampling_depth <- sampling_depth > $sampling
comparison <- comparison[, sampling_depth]

# Metadata information
group <- metadata[, col][match(colnames(comparison), metadata[,1])]
group <- factor(group, levels = c("$ref_group", "$group2"))

# Per sample normalization to 1M reads (as Galaxy does)
colsums<-colSums(comparison)

for (i in 1:length(colsums)) {
  comparison[, i]<- (comparison[, i]/colsums[i])*1000000
}

# Modifying input format
comparison <- SummarizedExperiment(assays = SimpleList(counts = comparison),
                               colData = DataFrame(grupo = group))

## LefSe calculation
Lefse<-lefser(
  expr = comparison,
  kruskal.threshold = 0.05,
  wilcox.threshold = 0.05,
  lda.threshold = 2,
  groupCol = "grupo",
  blockCol = NULL,
  assay = 1,
  trim.names = F
)

write.table(Lefse, file = paste("$LEfSe", "/LEfSe_output.tsv"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = F, 
            col.names = T)
            
RSCRIPT
