###### SCRIPT FOR OBTAINING THE CN DATA OF THE CELL LINES FROM THE DATABASE #####
#######################################  INFO ############################################
#QDNASeqCopyNumber generated using https://bitbucket.org/britroc/cnsignatures/src/master/#
#Copy number profiles from 1421 obtained with ASCAT                                      #                     
##########################################################################################

####### LOAD PACKAGES & ENVIRONMENT ######
#Directorios
inputDir <- "C:/Users/bchavesu/Desktop/datos_sWGS/CDKs/GRCH37/downstream_analyses/britroc-cnsignatures-bfb69cd72c50"


options(java.parameters = "-Xmx8000m")

#Cargar librerias
library(QDNAseqmod)
library(dplyr)
library(splitstackshape)
library(reshape2)
library(stats)


##### LOAD DATA  #####
##### Cell lines copy number obtained with ASCAT #####
cells_data <- readRDS(paste0(inputDir, "/cellLine_segment_ascat_sc_fixed_purity_tCN.rds"))
cells_data <- cells_data[cells_data$chromosome != "X",] #remove X chromosome
cells_data$chromosome <- as.integer(cells_data$chromosome) 

cells_map <- read.delim(paste0(inputDir, "/cellLine_CEL_file_mapping.tsv"), sep="\t", head=T)

##### Blas cell lines #####
blas_data <- readRDS(paste0(inputDir, "/cdks_30_kb_14_30kb_ds_absCopyNumber.rds"))
blas_cn   <- blas_data@assayData$copynumber #data from 6 different cell lines

#add genomic coordenates as columns
name_rows <- as.data.frame(row.names(blas_cn))
df <- cSplit(name_rows,"row.names(blas_cn)",":")
df <- cSplit(df,"row.names(blas_cn)_2","-")
colnames(df) <- c("chromosome", "start", "end") 
blas_cn <- cbind(blas_cn,df)
blas_cn$binID <- 1:nrow(blas_cn) #add ID for each 30kb bin


##### GET CN VALUES PER 30KB BIN IN ALL CELL LINES #####
samples <- unique(cells_data$sample)
nsamp <- length(unique(cells_data$sample))
nbins <- nrow(blas_cn)

cn_matrix <- matrix(nrow = nbins, ncol = nsamp) 
colnames(cn_matrix) <- samples

#for each sample, get data from cell lines within this genomic range
cn <- c()
i=1
for (i in 1:nrow(blas_cn)){
  print(paste0("Starting bin #", i))
  chrom <- blas_cn$chromosome[i]
  start <- blas_cn$start[i]
  end   <- blas_cn$end[i]
  
  cn <- cells_data[(cells_data[,1] %in% chrom & cells_data[,2] <= start & cells_data[,3] >= end), ]
  s=1
  for (s in 1:nsamp){
    #select rows in cells_data that overlaps positions of row[i]
    samp <- samples[s]
    segVal <- cn[cn$sample==samp, "segVal"]
    cn_matrix[i,s] <- ifelse(length(segVal)!=0, round(segVal), NA) ##ROUND SEG VAL !!!!!
  }
}
# Save .rds file with the copy number information of all cell lines
setwd(inputDir)
saveRDS(cn_matrix, "CNMatrix_CellLines_round.rds")