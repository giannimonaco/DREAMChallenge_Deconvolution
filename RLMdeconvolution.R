
###############################################################################
## This code uses Robust Linear Modelling to perform deconvolution
## It works for both the fine-grained and coarse-grained sub-Challenges defined 
## by the DREAM challenge on deconvolution opened in 2019.
###############################################################################
## Author: Gianni Monaco
## Date: 18/03/2022
###########


####################################################
################  LOAD REQUIRED LIBRARIES   


require(optparse, quietly = TRUE, warn.conflicts=FALSE)
require(MASS, quietly = TRUE, warn.conflicts=FALSE)
require(dplyr, quietly = TRUE, warn.conflicts=FALSE)
require(readr, quietly = TRUE, warn.conflicts=FALSE)
require(tidyr, quietly = TRUE, warn.conflicts=FALSE)
require(tibble, quietly = TRUE, warn.conflicts=FALSE)



option_list = list(
  make_option(c("--ExprData"), type="character", default=NULL, 
              help="Count data as CPM values", metavar="character"),
  make_option(c("--Submission"), type="character", default="third", 
              help="Choose between first, second or third", metavar="character"),   
  make_option(c("--SubChallenge"), type="character", default="coarse",
              help="Define granularity of deconvolution as defined by the DREAM challenge. Choose between coarse or fine.", metavar="character"),
  make_option(c("--Cancer"), type="character", default=NULL, 
              help="Choose between CRC (Colon Cancer) or BRCA (Breast Cancer)", metavar="character"),
  make_option(c("--OutFormat"), type="character", default="wide", 
              help="Output results can be saved in either wide or long format. Choose between wide or long", metavar="character"),
  make_option(c("--PrefixOutput"), type="character", default="", 
              help="Prefix to add to the output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$ExprData)){
  print_help(opt_parser)
  stop("You need to supply the count data as CPM values (ExpressionData).n", call.=FALSE)
}


######## load signature matrices
if(opt$SubChallenge == "coarse"){
  if(       opt$Submission =="first"){
    load("SignatureMatrix/SigMatCoarseCPM_1st.RData")
  }else if( opt$Submission =="second"){
    load("SignatureMatrix/SigMatCoarseCPM_2nd.RData")
  }else if( opt$Submission =="third"){
    load("SignatureMatrix/SigMatCoarseCPM_3rd.RData")
  }
  SigMat0 <- SigMatCoarse
}

if(opt$SubChallenge == "fine"){
  if(      opt$Submission =="first"){
    load("SignatureMatrix/SigMatFineCPM_1st.RData")
  }else if( opt$Submission =="second"){
    load("SignatureMatrix/SigMatFineCPM_2nd.RData")
  }else if( opt$Submission =="third"){
    load("SignatureMatrix/SigMatFineCPM_3rd.RData")
  }
  SigMat0 <- SigMatFine
}

####### Obtain cancer type
cancer_type <- opt$Cancer


if(      !is.null(cancer_type) && cancer_type == "CRC"){
  SigMat <- SigMat0[, grep("BRCA", colnames(SigMat0), invert = T)]
}else if(!is.null(cancer_type) && cancer_type == "BRCA") {
  SigMat <- SigMat0[, grep("CRC", colnames(SigMat0), invert = T)]
}else{
  SigMat <- SigMat0[, grep("CRC|BRCA", colnames(SigMat0), invert = T)]
}
SigMat <- as.matrix(SigMat)  



####### READ Input dataset
# This reads in the input file and converts to a matrix which will be
# input to the RLM deconvolution
expression_df <- opt$ExprData %>% 
  readr::read_csv( progress = FALSE, show_col_types=FALSE) %>% 
  as.data.frame() 

# check expression data
if(any(colSums(expression_df[,2:ncol(expression_df), drop=F]) !=10^6)){ 
  warning("Not all column values sum up to 1e+06. Make sure your expression table is in Counts Per Million (CPM)")
  #  expression_matrix <- apply(expression_matrix, 2, function(x) x/sum(x) * 10^6)
}

genesComm <- intersect(expression_df[,1], rownames(SigMat))

expression_df <- expression_df[expression_df[,1] %in% genesComm, ]
colnames(expression_df)[1] <- "Gene"
rownames(expression_df) <- NULL

expression_matrix <- expression_df %>%
  tibble::column_to_rownames("Gene") %>% 
  as.matrix() 


##### Perform Deconvolution

#check for NA data
if(sum(is.na(expression_matrix) >0 )) { expression_matrix[is.na(expression_matrix)] = 0 }

##### Robust Linear modelling
# We are using the HUGO version of the expression file
result_matrix <- apply(expression_matrix[genesComm,], 2, function(x) coef(rlm( SigMat[genesComm, ], x, maxit =100  ))) *100
result_matrix[result_matrix < -2] <- 0
result_matrix <- result_matrix + abs(min(result_matrix))
# result_matrix <- apply(result_matrix, 2, function(x) x/sum(x))   # sum to 1
result_matrix <- apply(result_matrix, 2, function(x) x/100)   # NO sum to 1
# result_matrix <- result_matrix[ grep("CRC|BRCA", rownames(result_matrix), invert = T), ]

# Convert the result matrix back to a dataframe
result_df <- result_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell.type") %>% 
  dplyr::as_tibble()

# Stack the predictions into one column
result_df_stacked <- tidyr::gather(
  result_df,
  key = "sample.id", 
  value = "prediction", 
  -cell.type) 



## Write result into output directory

if(      opt$OutFormat == "wide"){
  readr::write_csv(result_df, paste0(opt$PrefixOutput, "predictions.csv"))
}else if(opt$OutFormat == "long"){
  readr::write_csv(result_df_stacked, paste0(opt$PrefixOutput, "predictions.csv"))
}



