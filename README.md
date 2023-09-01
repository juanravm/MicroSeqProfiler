# MicroSeqProfiler

## Description

*MicroSeqProfiler* package is a powerful tool for in-depth analysis and exploration of microbiome data based on 16S ribosomal RNA marker genes. This package provides a complete and efficient solution for researchers to understand the composition and diversity of microbial communities present in different biological groups, by compiling a variety of previously developed tools.

With *MicroSeqProfiler*, you can perform comprehensive analyses of 16S RNA sequencing data, describing and quantifying the microbial species present in your samples. You can obtain diversity metrics to determine the divergences in the species in the microbial community and their relative abundance between sample groups by using [QIIME2 software](https://qiime2.org/).

Also, it performs LEfSe analysis to identify differentially abundant microorganisms between groups that can be associated with sample condition. Finally, it uses PICRUSt2 software to infer in the metagenomic composition of the microbial community and performs advanced differential abundance methods to reveal microbial metabolic pathways disparities between samples.

## Installation

As this package depends on other packages and softwares, you must visit their official website if you want acquire their latest versions:

-   [QIIME2 software](https://docs.qiime2.org/): QIIME2 is a python based bioinformatic tool that compiles different metagenomic programs that lets you analyze raw 16S reads of your sample. It is really useful for starting your metagenomic analysis and the diversity metrics.

-   [decontam package](https://github.com/benjjneb/decontam): *decontam* is an R package that performs a statistical OTU decontamination based on negative controls read counts.
  
-   [Lefser package](https://github.com/waldronlab/lefser): *Lefser* is an R package that performs taxonomic differential abundance analysis to identify differentially abundant microorganisms and infers on the effect of the abundance change in the group of samples. Bioconductor: (https://www.bioconductor.org/packages/release/bioc/html/lefser.html)

-   [PICRUSt2 software](https://github.com/picrust/picrust2): PICRUSt2 is also a python based tool that infers the rest of the metagenome from the 16S marker gene sequences. This inference is performed for KEGG Orthologs and MetaCyc databases, allowing us to carry out differential abundance methods.

-   [*ggpicrust2* package](https://github.com/cafferychen777/ggpicrust2): *ggpicrust2* is an R package used for the these differential abundance methods and the graphical illustration of the differentially active microbial metabolic pathways.

Each package has his own requirements and the installation method is explained in more detail in the [GitHub wiki web tab](https://github.com/juanravm/MicroSeqProfiler/wiki).

## Example data

For a better understanding of these functions, we prepared example data inputs for each script (a directory with example inputs for each function). You can obtain this example data through the following [link to Google Drive.](https://drive.google.com/drive/folders/10HhEDj57BGYD3aDi5E_5vNvG5jb54V4Q?usp=drive_link "ExampleData")
