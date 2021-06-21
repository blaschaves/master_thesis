###### SCRIPT FOR COMPARING CN PROFILES BETWEEN QDNAseqCopyNumber DATA AND 1421 CELL LINES #####

#######################################  INFO ############################################
#QDNASeqCopyNumber generated using https://bitbucket.org/britroc/cnsignatures/src/master/#
#Copy number profiles from 1421 obtained with ASCAT                                      #                          
##########################################################################################

####### LOAD PACKAGES & ENVIRONMENT ######

#Choose your input directory
inputDir <-  "C:/Users/bchavesu/Desktop/datos_sWGS/CDKs/GRCH37/downstream_analyses/britroc-cnsignatures-bfb69cd72c50"
options(java.parameters = "-Xmx8000m")

#Loading of the diferent libraries
library(QDNAseqmod)
library(dplyr)
library(splitstackshape)
library(reshape2)
library(stats)
library(corrplot)
library(RColorBrewer)
library(ggplot2)


##### LOAD DATA  #####
##### Cell lines #####
cn_matrix <- readRDS(paste0(inputDir,"/CNMatrix_CellLines.rds"))

##### QDNAseqCopyNumber cells data #####
blas_data <- readRDS(paste0(inputDir, "/cdks_30_kb_14_30kb_ds_absCopyNumber.rds"))
blas_cn   <- round(blas_data@assayData$copynumber) #Extract the copy number data from the different cell lines

#add genomic coordenates as columns
name_rows <- as.data.frame(row.names(blas_cn))
df <- cSplit(name_rows,"row.names(blas_cn)",":")
df <- cSplit(df,"row.names(blas_cn)_2","-")
colnames(df) <- c("chromosome", "start", "end") 

blas_cn <- cbind(blas_cn,df)
blas_cn$binID <- 1:nrow(blas_cn) #add ID for each 30kb bin
cn18k1016kp <- blas_cn[,c(1,7:10)]
cn18K2      <- blas_cn[,c(2,7:10)]
cn18kp      <- blas_cn[,c(3,7:10)]
cnEVA1      <- blas_cn[,c(4,7:10)]
cnEVB      <- blas_cn[,c(5,7:10)]
cn18K10     <- blas_cn[,c(6,7:10)]

#Prepare dataframe from  QDNAseqCopyNumber cell lines for correlation analyses
cn18k1016kp$sample <- "18k1016kp"
colnames(cn18k1016kp)[1] <- "segVal"
cn18K2$sample <- "18K2"
colnames(cn18K2)[1] <- "segVal"
cn18kp$sample <- "18kp"
colnames(cn18kp)[1] <- "segVal"
cnEVA1$sample <- "EVA1"
colnames(cnEVA1)[1] <- "segVal"
cnEVB$sample <- "EVAB"
colnames(cnEVB)[1] <- "segVal"
cn18K10$sample <- "18K10"
colnames(cn18K10)[1] <- "segVal"


#### CORRELATION QDNAseqCopyNumber CELL LINES VERSUS THE 1421 CELL LINES #####

#FOR 18K2 CELL LINE
res_cor <- c() #Define the vector where the correlation test results will be recorded
cells   <- ncol(cn_matrix) 


i=1
for(i in 1:cells){
  cor <- cor.test(cn18K2$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_cor18K2 <- as.data.frame(res_cor)
res_cor18K2$q.value <- p.adjust(res_cor18K2$p.value,method="BH") #adjust p.values by multiple comparisons

#FOR 18kp CELL LINE
res_cor <- c()
cells   <- ncol(cn_matrix)

i=1
for(i in 1:cells){
  cor <- cor.test(cn18kp$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_cor18kp <- as.data.frame(res_cor)
res_cor18kp$q.value <- p.adjust(res_cor18kp$p.value,method="BH") #adjust p.values by multiple comparisons


#FOR EVA1 CELL LINE
res_cor <- c()
cells   <- ncol(cn_matrix)

i=1
for(i in 1:cells){
  cor <- cor.test(cnEVA1$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_corEVA1 <- as.data.frame(res_cor)
res_corEVA1$q.value <- p.adjust(res_corEVA1$p.value,method="BH") #adjust p.values by multiple comparisons

#FOR EVB CELL LINE
res_cor <- c()
cells   <- ncol(cn_matrix)

i=1
for(i in 1:cells){
  cor <- cor.test(cnEVB$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_corEVB <- as.data.frame(res_cor)
res_corEVB$q.value <- p.adjust(res_corEVB$p.value,method="BH") #adjust p.values by multiple comparisons

#FOR 18K10 CELL LINE
res_cor <- c()
cells   <- ncol(cn_matrix)

i=1
for(i in 1:cells){
  cor <- cor.test(cn18K10$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_cor18K10 <- as.data.frame(res_cor)
res_cor18K10$q.value <- p.adjust(res_cor18K10$p.value,method="BH") #adjust p.values by multiple comparisons

#FOR 18K1016kp CELL LINE
res_cor <- c()
cells   <- ncol(cn_matrix)

i=1
for(i in 1:cells){
  cor <- cor.test(cn18k1016kp$segVal, cn_matrix[,i], method = "pearson")
  cor <- cbind(cellLine=colnames(cn_matrix)[i],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  res_cor <- rbind(res_cor, cor)
}

res_cor18k1016kp <- as.data.frame(res_cor)
res_cor18k1016kp$q.value <- p.adjust(res_cor18k1016kp$p.value,method="BH") #adjust p.values by multiple comparisons


#Matrix with all correlations
cor_matrix <- cbind(as.numeric(res_cor18kp[,2]), as.numeric(res_cor18K2[,2]), 
                    as.numeric(res_cor18K10[,2]), as.numeric(res_cor18k1016kp[,2]), 
                    as.numeric(res_corEVA1[,2]), as.numeric(res_corEVB[,2]))
colnames(cor_matrix)  <- c("18kp", "18K2","18K10", "18k1016k", "EVA1", "EVAB")
row.names(cor_matrix) <- res_cor18K10$cellLine

#Donwload correlation matrix data as a .csv file

write.csv(cor_matrix,"correlation_results.csv")

###CORRELATION TEST BETWEEN ALL CELL LINES DATA ###

cor_all <- c()
ndata <- ncol(cn_matrix)

i=1
j=1

for(i in 1:ndata){
     for (j in 1:ndata){
    cor <- cor.test(cn_matrix[,i], cn_matrix[,j], method = "pearson")
    cor <- rbind(cellLine_1=colnames(cn_matrix)[i],cellLine_2=colnames(cn_matrix)[j],cor=cor$estimate,p.value=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
    cor_all <- cbind(cor_all, cor)
    print(paste0("Comparing n#  ", i))
   }
  }

cor_all <- as.data.frame(cor_all)
cor_all$q.value <- p.adjust(cor_all$p.value,method="BH") #adjust p.values by multiple comparisons
setwd(inputDir)
saveRDS(cor_all, "correlations_matrix.rds")

##Represent the correlation between all the ASCAT cell lines
corr_coeficients <- correlations_matrix[3,]

setwd(inputDir)
saveRDS(corr_coeficients, "corr_coeficients.rds")

  ggplot(aes(x=corr_coeficients))+
  geom_density(fill='firebrick4')+
  xlab("Pearson correlation coefficient")+
  ylab("Frequency")+
  ggtitle("Distribution of the different correlation coeficients")+
  theme_minimal()
  



###### PLOTS #######
#Plot heatmap of correlations
#Plot data using heatmap 
library(nlme)
library(ComplexHeatmap)

#heatmap with cluster

setwd(inputDir)
tiff("CorrelationHeatmap.tiff", width = 8, height = 16, res=600, units = "in")
Heatmap(cor_matrix, row_names_gp = gpar(fontsize = 3), 
        column_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_names = T, 
        show_column_names = T, cluster_rows = T, cluster_columns = T,
        clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
        #width = unit(150, "mm"), 
        heatmap_legend_param = list(title = "Correlation", title_gp = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 8)))
dev.off()

## Example of heatmap for showing the correlation results of one cell lines vs all QDNAseqCopyNumber data
# cor.matrix.HuH7 <- cor_matrix["HuH-7",]
# tiff("CorrelationHeatmap_HuH7.tiff", width = 10, height = 10, res = 600, units = "in")
# Heatmap(cor.matrix.HuH7, row_names_gp = gpar(fontsize = 12), 
#         column_names_gp = gpar(fontsize = 0), row_names_side = "left", show_row_names = T, 
#         show_column_names = T, cluster_rows = T, cluster_columns = T,
#         clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
#         width = unit(150, "mm"), height = unit(150, "mm"),
#         heatmap_legend_param = list(title = "Correlation", title_gp = gpar(fontsize = 10),
#                                     labels_gp = gpar(fontsize = 8)))
# dev.off()



#Plot a copy number profile  from a selected cell line
cells_data$segVal <- round(cells_data$segVal)
ncih23 <-cells_data[cells_data$sample=="NCI-H23",]
colnames(ncih23)[4] <- "cn" # plot_from_segTab.R requieres SegVal data named as "cn"

setwd(inputDir)
source(paste0(inputDir,"/CNPlot_from_segTab.R"))
pdf("ncih23_CNprofile.pdf",width = 20, height = 15)
CNPlot_from_segTab.R(ncih23) 
dev.off()

#Plot a copy number profile from a selected QDNAseqCopyNUmber cell line
segkp <- getSegTable(blas_data [,3]) # plot_fromsegTab.R requieres segments data
colnames(segkp)[4] <- "cn" 
segkp$cn <- round(as.numeric(segkp$cn))
segkp$cn <- segkp$cn 

segk2 <- getSegTable(blas_data [,2]) # plot_fromsegTab.R requieres segments data
colnames(segk2)[4] <- "cn" 
segk2$cn <- round(as.numeric(segk2$cn))
segk2$cn <- segk2$cn + 0.1 ## Values lift up 0.1 in order to see diferences between the 2 CN profiles



setwd(inputDir)
source(paste0(inputDir,"/CNPlot_from_segTab.R"))
pdf("segeva1_CNprofile.pdf",width = 20, height = 15)
CNPlot_from_segTab.R(segeva1) 
dev.off()


#Plot both copy number profiles at the same time
setwd(inputDir)
source(paste0(inputDir,"/CNPlot_from_segTab.R"))
pdf("kp_k2.pdf",width = 20, height = 15)
CNPlot_from_segTab.R(segkp,segk2) ## OMIT IN plot_from_segTab.R the data of chromosomes X and Y !!!!!!!!!
dev.off()

## Plot both copy number profiles of a selected chromosome
setwd(inputDir)
source(paste0(inputDir,"/CNPlot_by_chrom.R"))
pdf("cromosoma.pdf",width = 20, height = 15)
chr_number <- 22
CNPlot_by_chrom.R(segkp,segk2,chr_number) ## OMIT IN plot_from_segTab.R the data of chromosomes X and Y !!!!!!!!!
dev.off()
  

