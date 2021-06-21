###### SCRIPT FOR FIXING SEGMENT BOUNDARIES ACROSS CELLS #########
## This script is used for for fixing segment boundaries across### 
## cells after obtaining the unique copy number changes###########

###Calculate the unique copy number events

inputDir<- "C:/Users/bchavesu/Desktop/datos_sWGS/CDKs/GRCH37/downstream_analyses/britroc-cnsignatures-bfb69cd72c50"
source(paste0(inputDir,"/unify_segments_boundaries.R"))
fix <- unify_segment_boundaries(abscn_18Kp_18K2,500000,30000)

### Plot the unique copy number events

## Manipulate the data for plotting
unique_K2<-fix[,1:4]
colnames(unique_K2)[4]<-"cn"
unique_K2$cn <- as.numeric(unique_K2$cn)
unique_K2$cn <- round(unique_K2$cn)

unique_Kp<-fix[,c(1,2,3,5)]
colnames(unique_Kp)[4]<-"cn"
unique_Kp$cn <- round(as.numeric(unique_Kp$cn))


##Plot unique unique copy number from both data

setwd(inputDir)
source(paste0(inputDir,"/CNPlot_from_segTab.R"))
pdf("18K_unique.pdf",width = 20, height = 15)
CNPlot_from_segTab.R(unique_K2,unique_Kp) ## OMIT IN plot_from_segTab.R the data of chromosomes X and Y !!!!!!!!!
dev.off()

## PLot both data selecting which is different between the pool and the single ko
setwd(inputDir)
source(paste0(inputDir,"/plot_dif_cn.R"))
pdf("18K_unique_diferences_round.pdf",width = 20, height = 15)
plot_dif_cn(unique_K2,unique_Kp) ## OMIT IN plot_from_segTab.R the data of chromosomes X and Y !!!!!!!!!
dev.off()

## PLot both data selecting which is different between the pool and the single ko for all the chromosomes
for (i in 1:22){
setwd(inputDir)
source(paste0(inputDir,"/plot_dif_by_chrom.R"))
pdf(paste0("18K_unique_diferences_chr_",i,"_round.pdf"),width = 20, height = 15)
chr_number<-i
p<-plot_dif_by_chrom(unique_K2,unique_Kp,chr_number) ## OMIT IN plot_from_segTab.R the data of chromosomes X and Y !!!!!!!!!
print(p)
dev.off()
}






