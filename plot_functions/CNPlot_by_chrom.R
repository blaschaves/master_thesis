##### FUNCTION FOR PLOTTING A CN PROFILE FROM  1 or 2 SEGTAB DATA SET CHROMOSOME BY CHROMOSOME ##############
# Function based on plot_from_segTab.R                                                                      #
# To call this function: CNPlot_from_segTab.R(<segTabb1>,<segTab>)                                          #
#############################################################################################################

# To plot a CN profile chromosome by chromosome: CNPlot_by_chrom.R
### Ploting function
CNPlot_chrom <- function(events,events_2,color=F){
  #Determine the Y limit for not to loose any of the cn 
  if (max(events$cn)==max(events_2$cn)){
    y_axis<-round(max(events$cn))
  }
  else if (max(events$cn)>max(events_2$cn)){
    y_axis<-round(max(events$cn))
  }
  else{
    y_axis<-round(max(events_2$cn))
  }
  par(pin=c(14,10),las=1)
  plot.new()
  plot.window(xlim=c(0,sum(events$length)),ylim=c(0,y_axis))
  abline(h=0,lty=1)
  abline(h=1:y_axis,lty=2,col="purple")
  title(main="Copy Number Profile",xlab=paste0("Chromosome   ",chr_number),ylab="Copy Number",
        cex.lab=2,cex.main=2,font.lab=2,family="Palatino")
  ##Set the axis
  axis(1,labels=FALSE, tick=F,
       at = chr_sizes$offset+chr_sizes$length/2,cex.axis=1)
  # labels = chr_number
  # lablist<-as.vector(c(0:6))
  # text(par("usr")[4], seq(0, 6, by=1), labels = lablist, srt = 0, pos = 2, xpd = TRUE, offset = 1 ,
  #      cex=1.3)
  
  axis(2,labels=as.vector(c(0:y_axis)),
       at=seq(0, y_axis, by=1))
  
  cols = colorRampPalette(c("Green","Red"))(10)
  for (i in 1:nrow(events)){
    e <- events[i,]
    if (color == T){
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,
               col=cols[[1+round(9*(abs(round(e$cn) - e$cn)/0.5))]])
    }
    else{
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,col="darkgreen")
    }
  }
  
  # Plot the second data set 
  for (i in 1:nrow(events_2)){
    e <- events_2[i,]
    if (color == T){
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,
               col=cols[[1+round(9*(abs(round(e$cn) - e$cn)/0.5))]])
    }
    else{
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,col="darkblue")
    }
  }
}

# Call functions and set chromosomes sizes

CNPlot_by_chrom.R<-function(cn_profile,cn_profile_2,chr_number)
{
  chr_number <- chr_number
## Choose chromosome and set stablish CNprofile data for that chromosome
  chr_sizes <- read.csv(file="hg19.chrom.sizes.txt",sep="\t",stringsAsFactors = F,
                        header=FALSE,col.names = c("chr","length"))[chr_number,] 
  chr_sizes$length <- 1e-6 * chr_sizes$length # Convert sizes to Mb
  chr_sizes$offset <- cumsum(chr_sizes$length) - chr_sizes$length # Calculate the interval range of each chromosome
## Obtain the values for only the selcted chromosome
  chrSize_cn1 <- cn_profile[cn_profile[,1]==chr_number,]
  chrSize_cn2 <- cn_profile_2[cn_profile_2[,1]==chr_number,]
  CNPlot_chrom(CNconvert_chrom(chrSize_cn1,chr_sizes),CNconvert_chrom(chrSize_cn2,chr_sizes))
}

### Convert function -> Prepare data to plot it
CNconvert_chrom <- function(cs,base){
  cs$start <- as.numeric(cs$start)
  cs$end <- as.numeric(cs$end)
  cs$cn <- as.numeric(cs$cn)
  cs$start <- 1e-6 *cs$start
  cs$end <- 1e-6 * cs$end
  for (c in unique(cs$chr)){
    cs$start[cs$chr == c] <- cs$start[cs$chr == c]+ 
      base$offset[base$chr == c]
    cs$end[cs$chr == c]   <- cs$end[cs$chr == c]  + 
      base$offset[base$chr == c]
  }
  cs$length <- cs$end - cs$start
  return(cs)
}


