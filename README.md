# RNAseq PostProcesing
## Leticia G Leon

### Introduction

Those rmarkdown code were created to perform post-processing analysis of the RNAseq data, from Frozen or FFPE samples.

There is NOT a pipeline or workFlow but we recomend an order to use those scripts

The suggested order would be:

- Filtering script, to remove the Noise in the counts files removing the low counts rows
- Normalization script, to calculate rlog, rpkm, fpkm or tpm 
- Correlation, and  Zscore script, those one can be use before and after normalization, to check teh correlation between samples and unsupervised clustering (before and after Norm)
- Overview script, in this one we used DESeq2 tool and some edgeR to take a overview (clustering with different distance, PCA and other techniques)
- DE_analysis with DESeq2
- DE_analysis with edgeR
- CMSclasification/CMScaller
- immunedeconv(MCP_counter, EPIC)
- Select a set of genes and perform some visualization



The input for those script is always a count matrix. This matrix could come form featureounts or HTseq or any other tool for gene expression cuantification 


### Installation

There is no instalation needed, you need to have R installed to run those scrpts
Immunedeconv needs a bit more attention (https://icbi-lab.github.io/immunedeconv/)

**Required R-version etc.**

Those script were make and tested in R 3.6, using BiocManager to install the neccesary packages from Bioconductor or Cran. Be aware that you need admin permision to install packages


### Usage

The best way to use those script is creating a copy in your working directory.

Those stepts are not interctive, yo need to take care of few features

- path to your working directory, where you want to save the outputs
- Path to the raw counts files 
- Indicate the path to your metadata


#### To Run those scripts: 

__R scripts__

This is a R markdown script, the msot friendly way to modify is in Rstudio. Once you modifly the path you can use the option Knit to run the entire code and generate an html(pdf etc) to keep a record of the settings that you ran for each study. 

If you prefer not to depend on Rstudio, you can modify the paths and bin size using any editor that you like and run it in CLI with this option:
R -e "rmarkdown::render('/path to the Rmarkdown script /testCLI_Knit/test_knitCLI.Rmd')"


#### Input files

- Counts
Most of those script just need the raw counts files or the filter count files, for the Overview or visualization use teh Norm files that you generate with teh filter and Norm scripts

- Metadata
You will also need a metadata sheet, csv, txt (No excel), this metadata have the samples as row names in the exact way as they are in columns of the count matrix. This metadata contain as much phenotipic information as possible and it is as complete as possible

#### How to run 

Run in Rstudio using the Knit option or cell by cell
Run the the terminal:
R -e "rmarkdown::render('/path to the Rmarkdown script /testCLI_Knit/test_knitCLI.Rmd')"


#### Output Files
___From the Filter, You will get:___
A count matrix without the very low counts rows, the thershold can be change in the script

___From Norm, you will get:___
You can obtain a Normalized counts matrix, the normalization that you can calculate are: 
rlog (using DESeq2)
rpkm (using edgeR function)
fpkm and tpm, manual calculations

___Overview and Correlations, Zscore___
Several pdf witht he corelation plots, clustering analysis, using DESeq2 functions or normla hclus in R, PCA plots ans analysis

___DE___
A table with the most differenciate genes, order acording to the FDR (adj.pval)
several other plots with heatmap with the most variable genes, heatmap (pdf format)
Bed file 

___CMScaller___
Csv file with the probabilities of that sample belonging to an specific CMS group, the first column os the prediction column
cluster plot using those probabilities. As a visualization you will get a heatmap, genset enrichment plot and PCA.

___CMSclasifier___
Csv file with the probabilities of that sample belonging to an specific CMS group, the first column os the prediction column
cluster plot using those probabilities

___immunedeconv(MCP_counter, EPIC)___
Bar and scater plot with the different cell clasifications
