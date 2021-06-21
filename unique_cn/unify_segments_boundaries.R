##### FUNCTION FOR UNIFYING SEGMENT BOUNDARIES ################################                                      
# To call this function: unify_segment_boundaries(<CN>,<window>,<bin_width>)###
# This function is used to determine the unique copy number events          ###
###############################################################################

unify_segment_boundaries<-function(cn,window,bin_width)
{
  segTabs<-list()
  for(i in colnames(cn))
  {
    segTabs[[i]]<-getSegTable(cn[,i])
  }
  
  starts<-cbind(as.numeric(unlist(lapply(segTabs,"[[",1))),as.numeric(unlist(lapply(segTabs,"[[",2))))
  ends<-cbind(as.numeric(unlist(lapply(segTabs,"[[",1))),as.numeric(unlist(lapply(segTabs,"[[",3))))
  
  starts<-unique(starts[order(starts[,1],starts[,2]),])
  ends<-unique(ends[order(ends[,1],ends[,2]),])
  
  #collapse start segments that are within window distance
  starts_collapsed<-c()
  for(i in 1:(nrow(starts)-1))
  {
    if(abs(starts[i+1,2]-starts[i,2])>(window+1))
    {
      starts_collapsed<-rbind(starts_collapsed,starts[i,])
    }
  }
  
  #create segment ends vector
  ends_collapsed<-c()
  for(i in unique(starts_collapsed[,1]))
  {
    ends_collapsed<-c(ends_collapsed,starts_collapsed[starts_collapsed[,1]==i,2][-1]-1)
    ends_collapsed<-c(ends_collapsed,max(ends[ends[,1]==i,2]))
  }
  
  new_seg<-cbind(starts_collapsed,ends_collapsed)
  
  #generate output matrix of unified copy number
  cn_out<-c()
  for(i in colnames(cn))
  {
    curr_seg<-c()
    for(j in 1:nrow(new_seg))
    {
      chr<-new_seg[j,1]
      start<-new_seg[j,2]
      end<-new_seg[j,3]
      curr_cn<-assayDataElement(cn[,i],"copynumber")
      start_ind<-which(rownames(curr_cn)==paste0(chr,":",start,"-",as.integer(start+bin_width-1)))
      end_ind<-which(rownames(curr_cn)==paste0(chr,":",as.integer(end)-bin_width+1,"-",as.integer(end)))
      if(length(end_ind)==0)
      {
        end_ind<-which(grepl(paste0(chr,":[0-9]+-",as.integer(end)),rownames(curr_cn)))
        if(length(end_ind)==0)
        {
          end_ind<-which(rownames(curr_cn)==paste0(chr,":",end+1,"-",as.integer(end+bin_width)))
        }
      }
      curr_seg<-c(curr_seg,median(curr_cn[start_ind:end_ind,],na.rm=T))
    }
    cn_out<-cbind(cn_out,curr_seg)
  }
  colnames(cn_out)<-colnames(cn)
  colnames(new_seg)<-c("chromosome","start","end")
  out<-cbind(new_seg,cn_out)
  out<-data.frame(out,stringsAsFactors = F)
  out
}