revise_bed=function(bed, min_span=30, merge=TRUE){
  
  bed=read.table('inputdata/cnv.anno.bed')
  bed=bed[bed[,3]-bed[,2]>=min_span,]
  nu_chr=substr(bed[,1], 4, 7)
  nu_chr=gsub( 'X', 23, nu_chr)
  nu_chr=gsub( 'Y', 23, nu_chr)
  nu_chr=as.numeric(nu_chr)*1e10
  
  ###this sorts the starting position
  I=order(nu_chr+bed[,2])
  bed=bed[I,]
  nu_chr=nu_chr[I]
  
  ###make sure bed regions do not overmap, merge overlaps if any
  bed_list=list()
  tt=cbind(nu_chr+bed[,2], nu_chr+bed[,3])
  
  unique_range=tt[1,]
  for (i in 1:(nrow(tt)-1)){
    if (unique_range[2]>=tt[i+1, 1]){
      ###in case the previous record totally overlaps the next record
      unique_range=c(unique_range[1], max(unique_range[2], tt[i+1, 2]))
    } else {
      bed_list[[i]]=unique_range
      unique_range=tt[i+1,]
    }
  } 
  bed_list[[i+1]]=unique_range
  return(do.call(rbind, bed_list))
  
}



