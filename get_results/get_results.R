################ Get Copy-number signatures results and figures ###################################################
# This script is compiled from BRITROC CNsignatures which contains all code necessary to reproduce the analyses #
# for the accompanying manuscript ["Copy-number signatures and mutational processes in ovarian carcinoma"]        #                                                                                                                #                                                                                                            
###################################################################################################################

## First steps
# 1st Select your working directory. It has to include all data necessary:
# main_functions.R and helper_functions.R scripts
# .rds file generated with the Absolute Copy Number pipeline

## Select your path in main_functions.R script. 

  
##Set up 
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(message = FALSE)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(QDNAseq))
suppressMessages(library(flexmix))
suppressMessages(library(NMF))
num_cores<-16

## Run main_functions.R

source(paste(this_path,"main_functions.R",sep="/"))

##Loading the Absolute Copy Number .rds file into a variable
a <- readline("Enter .rds file: ")
abscn <- readRDS(a)
## Redefine abscn only with the information of the 18Kp and 18K2. The rest are discarded
abscn_18Kp_18K2 <- abscn[,2:3]

## Extracting CNfeatures
CN_features_18Kp_18K2 <- extractCopynumberFeatures(abscn_18Kp_18K2)

## Generate Sample By Component Matrix
sample_by_component <- generateSampleByComponentMatrix(CN_features_18Kp_18K2)

## Quantify signatures
signature_quantification <- quantifySignatures(sample_by_component)

## Plot Sample by Component Matrix
NMF::aheatmap(sample_by_component,fontsize = 7,Rowv=FALSE,Colv=FALSE,legend = T,breaks=c(seq(0,199,2),500),main="Component-by-sample matrix")

## Plot Patient x Signature matrix
nmfalg<-"brunet"
seed<-77777
nsig<-7
sigs<-NMF::nmf(t(sample_by_component),nsig,seed=seed,nrun=1000,method=nmfalg,.opt = paste0("p",num_cores))
coefmap(sigs,Colv="consensus",tracks=c("basis:"), main="Sample-by-signature",)

## Plot Signature x component matrix
basismap(sigs,Rowv=NA,main="Signature x Component matrix")




