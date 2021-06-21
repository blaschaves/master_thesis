##### FUNCTION FOR PLOTTING A CN PROFILECN FROM  1 or 2 SEGTAB DATA ###########
# Function based on plot_from_segTab.R                                        #
# To call this function: CNPlot_from_segTab.R(<segTabb1>,<segTab>)            #
###############################################################################

### To plot a CN profile call: CNplot_from_segTab.R###

### Ploting function
CNPlot_events <- function(events,events_2,color=F){
  par(pin=c(14,10),las=1)
  plot.new()
  plot.window(xlim=c(0,sum(events$length)),ylim=c(0,10),yaxs="i",xaxs="r")
  abline(v=chr_sizes$offset,lty=2,lwd=1)
  abline(v=0,lty=1,lwd=1)
  title(main="Copy Number Profile",xlab="Chromosome",ylab="Copy Number",
        cex.lab=1.5,cex.main=2,font.lab=2,family="Palatino",outer=F)
  ##Set the axis
  axis(1,labels = c(1:22),tick=T,
       at = chr_sizes$offset+chr_sizes$length/2,cex.axis=1)
  
  lablist<-as.vector(c(0:10))
  text(par("usr")[4], seq(0, 10, by=1), labels = lablist, srt = 0, pos = 2, xpd = TRUE, offset = 1 ,
       cex=1.3)
  
  # axis(2,labels=as.vector(c(0:6)),
  #      at=seq(0, 6, by=1))
  
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

### Call functions and set chromosomes sizes ###
CNPlot_from_segTab.R<-function(cn_profile,cn_profile_2)
{
  chr_sizes <- read.csv(file="hg19.chrom.sizes.txt",sep="\t",stringsAsFactors = F,
                        header=FALSE,col.names = c("chr","length"))[1:22,] ## [1:22,]
  chr_sizes$length <- 1e-6 * chr_sizes$length # Convert sizes to Mb
  chr_sizes$offset <- cumsum(chr_sizes$length) - chr_sizes$length # Calculate the interval range of each chromosome
  CNPlot_events(CNconvert(cn_profile,chr_sizes),CNconvert(cn_profile_2,chr_sizes))
}

### Convert function -> Prepare data to plot it
CNconvert <- function(e,base){
  e$start <- as.numeric(e$start)
  e$end <- as.numeric(e$end)
  e$cn <- as.numeric(e$cn)
  e$start <- 1e-6 * e$start
  e$end <- 1e-6 * e$end
  for (c in unique(e$chr)){
    e$start[e$chr == c] <- e$start[e$chr == c]+ 
      base$offset[base$chr == c]
    e$end[e$chr == c]   <- e$end[e$chr == c]  + 
      base$offset[base$chr == c]
  }
  e$length <- e$end - e$start
  return(e)
}
