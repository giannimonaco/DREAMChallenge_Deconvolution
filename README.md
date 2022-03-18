[//]: # (Title: Description Deconvolution method - Biogem team)  
[//]: # (Author: Gianni Monaco)  
[//]: # (Date: March 18, 2022) 



# DREAM Challenge on Deconvolution -- Biogem team

This repository contains the code used to participate to the Dream Challenge on RNA-seq deconvolution. 


---
## Pre-requisites

Install the following dependencies for your R instance:

* optparse
* MASS
* dplyr
* readr
* tidyr
* tibble 


## Usage

To use the method, download this repository to your local computer. 

Move inside the folder containing the file "RLMdeconvolution.R" and run the R script in the following way:

```
Rscript ./RLMdeconvolution.R --ExprData <path/to/input/CPMexpr.csv> [--SubChallenge <coarse,fine>] [--Cancer <CRC,BRCA>]
[-- OutFormat <wide,long>] [--PrefixOutput <path/to/output/Prefix>] [--Submission <first,second,third>]
```

required argument|description
---|---
--ExprData <path/to/input/CPMexpr.csv> | This is the only required argument. It is where you define your input gene expression data file. The files should be a csv file with the first column containing the gene symbols and the remaining column containing the gene expression of your samples as Counts Per Million (CPM) values.


optional argument|description
---|---
--SubChallenge <coarse,fine> | Choose Coarse to deconvolute the major cell types. It has less granularity but produces more robust results (B cells, CD4 T cells, CD8 T cells, NK cells, Neutrophils, Monocytic-lineage cells, Fibroblasts, Endothelial cells). Choose Fine to deconvolution more cell types (Naive B cells, Memory B cells, Naive CD4 T cells, Memory CD4 T cells, Regulatory T cells, Naive CD8 T cells, Memory CD8 T cells, NK cells, Neutrophils, Monocytes, Myeloid Dendritic cells, Macrophages, Fibroblasts and Endothelial cells).  
--Cancer <CRC,BRCA> | With this option, you can choose to deconvolute also cancer cells from Colon (CRC) or Breast (BRCA) cancer.  
-- OutFormat <wide,long> | You can choose if the results shoudl be formatted as wide or long table. The default parameter is wide.  
-- PrefixOutput <path/to/output/Prefix> | If you want to add a prefix to your result file add it through this argument. You could also add the path where you want to save your results.  
--Submission <first,second,third> |Choose which signature matrix you want to use between the ones developed for the three submissions to the DREAM challenge. The default parameter is third, as better results were obtained for the last submission.  


  




