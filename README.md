[//]: # (Title: Description Deconvolution method - Biogem team)  
[//]: # (Author: Gianni Monaco)  
[//]: # (Date: March 18, 2022) 



# DREAM Challenge on Deconvolution -- Biogem team

This repository contains the code used to participate to the Dream Challenge on RNA-seq deconvolution. 


---
## Pre-requisites

Install the following dependencies for your R (tested with v4.2.0) instance:

* optparse (tested with v1.7.1)
* MASS (tested with v7.3-57)
* dplyr (tested with v1.0.9)
* readr (tested with v2.1.2)
* tidyr (tested with v1.2)
* tibble (tested with v3.1.7)

## Installation
The installation time on a should be around 10min in a normal laptop or desktop computer.  

## Usage

To use the method, download this repository to your local computer. 

From the command line, navigate to the folder containing the file "RLMdeconvolution.R" and run the R script in the following way:

```
Rscript ./RLMdeconvolution.R --ExprData <path/to/input/CPMexpr.csv> [--SubChallenge <coarse,fine>] [--Cancer <CRC,BRCA>]
[-- OutFormat <wide,long>] [--PrefixOutput <path/to/output/Prefix>] [--Submission <first,second,third>]
```  
  

The only required argument is the gene expression data. The other arguments are optional. In the table below, you can find a descriptions of all arguments.  


argument|description
---|---
--ExprData <path/to/input/CPMexpr.csv> | This is the only required argument. It is where you define your input gene expression data file. The files should be a csv file with the first column containing the gene symbols and the remaining columns containing the gene expression of your samples normalized as Counts Per Million (CPM).
--SubChallenge <coarse,fine> | Choose Coarse to deconvolute the major cell types. It has less granularity but produces more robust results (B cells, CD4 T cells, CD8 T cells, NK cells, Neutrophils, Monocytic-lineage cells, Fibroblasts, Endothelial cells). Choose Fine to deconvolution more cell types (Naive B cells, Memory B cells, Naive CD4 T cells, Memory CD4 T cells, Regulatory T cells, Naive CD8 T cells, Memory CD8 T cells, NK cells, Neutrophils, Monocytes, Myeloid Dendritic cells, Macrophages, Fibroblasts and Endothelial cells).  
--Cancer <CRC,BRCA> | With this option, you can choose to deconvolute also cancer cells from Colon (CRC) or Breast (BRCA) cancer.  
--OutFormat <wide,long> | You can choose if the results shoudl be formatted as wide or long table. The default parameter is wide.  
--OutFilename <path/to/output/Filename.csv> | The filename to use for your output file. Default name is "Predictions.csv". You could also add the path where you want to save your results.  
--Submission <first,second,third> |Choose which signature matrix you want to use between the ones developed for the three submissions to the DREAM challenge. The default parameter is third, as better results were obtained for the last submission for the fine subchallenge.  

## Output
The output is a CSV table that by default is named "Predictions.csv", where columns are the samples, rows are the cell types and the values are cell type frequences.

## Usage with a Demo dataset
This repository includes a demo dataset in the folder Demo. You can run the demo dataset in the following way:
```
Rscript ./RLMdeconvolution.R --ExprData Demo/TESTdataset.csv  --OutFilename Demo/Demo_Predictions.csv
```
The output file "Demo_Prediction.csv" will be generated in the Demo directory. The run time on a normal laptop or desktop computer is only a few seconds (about 2 second on a Mac Book Pro with a 2.9 Ghz CPU and 16 GB of RAM).

## License
This software is released under the MIT License. 

©GianniMonaco
  




