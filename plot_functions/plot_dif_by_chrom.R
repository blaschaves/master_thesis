##### FUNCTION FOR PLOTTING A CN PROFILECN FROM  1 or 2 SEGTAB DATA and show differences by chromosome###########
# Function based on plot_from_segTab.R                                                                          #
# To call this function: plot_dif_by_chrom(<segTab1>,<segTab2>,chr_number)                                     #
#################################################################################################################
### Ploting function
CNPlot_chrom <- function(events,events_2,color=F){
  #Determine the Y limit for not to loose any of the cn 
  if (max(events$cn)==max(events_2$cn)){
    y_axis<-max(events$cn)
  }
  else if (max(events$cn)>max(events_2$cn)){
    y_axis<-max(events$cn)
  }
  else{
    y_axis<-max(events_2$cn)
  }
  ##Plot
  par(pin=c(14,10),las=1)
  plot.new()
  plot.window(xlim=c(0,sum(events$length)),ylim=c(0,y_axis))
  abline(h=0,lty=1)
  abline(h=1:y_axis,lty=2,col="purple")
  title(main="Copy Number Profile",xlab=paste0("Chromosome   ",chr_number),ylab="Copy Number",
        cex.lab=3,cex.main=3,font.lab=2,family="Palatino")
  ##Set the axis
  axis(1,labels=FALSE, tick=F,
       at = chr_sizes$offset+chr_sizes$length/2,cex.axis=3)
  # labels = chr_number
  # lablist<-as.vector(c(0:6))
  # text(par("usr")[4], seq(0, 6, by=1), labels = lablist, srt = 0, pos = 2, xpd = TRUE, offset = 1 ,
  #      cex=1.3)
  
  axis(2,labels=as.vector(c(0:y_axis)),
       at=seq(0, y_axis, by=1),cex.axis=2)
  
  cols = colorRampPalette(c("Green","Red"))(10)
  for (i in 1:nrow(events)){
    e <- events[i,]
    e2<- events_2[i,]
    if (color == T){
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,
               col=cols[[1+round(9*(abs(round(e$cn) - e$cn)/0.5))]])
    }
    else if (e$cn!=e2$cn){
      segments(e$start,e$cn+0.1,e$end,e$cn+0.1,lwd=5.,col="red")
      segments(e2$start,e2$cn,e2$end,e2$cn,lwd=5.,col="blue")
    }
    else {
      segments(e$start,e$cn,e$end,e$cn,lwd=5.,col="darkgreen")
      segments(e2$start,e2$cn,e2$end,e2$cn,lwd=5.,col="darkgreen")
    }
  }
}

# Call functions and set chromosomes sizes

plot_dif_by_chrom<-function(cn_profile,cn_profile_2,chr_number)
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


