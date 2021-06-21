### Function for extracting the unique event from a selected chr ####
# To call this function: unique_events_from_chr(segtable,chr_number)#
#####################################################################

unique_events_from_chr <- function(fix_segtable,chr_number){
  unique_chr_cn <-c()
  for (i in 1:nrow(fix_segtable)){
    
    if (fix_segtable$chromosome[i]==chr_number){
      
      unique_chr_cn <- rbind (fix_segtable[i,],unique_chr_cn)
      
    }
    else{}
  }
  return(unique_chr_cn)
}