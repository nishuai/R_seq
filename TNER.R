###############################################################
# File name: TNER.R
# detect ctDNA mutation from background of mutational signatures  
# Author: Ni Shuai nishuai@yahoo.com
# Inspired by: shibing Deng, Tao Xie (xietao2000@Pfizer.com)
# YuceBio R&D department
# Input: test_sample file in  atgc format, All_training file,
# output directory (gwd by default)
#################################################################
# args=c('training/18B99984XJ_1.umi.pileup.txt', 'trained_BMER.rds',  0.1)
args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input
# processing the input
#####################
patient_pileup_file=trained_BMER=input_alpha=output_dir=NULL
if (length(args)<2) {
  stop("Usage: Rscript TNER.R patient_pileup_file trained_model input_alpha output_dir", call.=FALSE)
} else if (length(args)>1) {
  # default output file
  patient_pileup_file=args[1]
  trained_BMER=readRDS(args[2])  # Average bkgrd rate file
  if (length(args)>2) input_alpha=as.numeric(args[3])  # 3nt bkgrd rate file
}
if (is.null(input_alpha)) input_alpha=0.001
if (is.null(output_dir)) output_dir='.'

if (dir.exists(patient_pileup_file)){
  patient_pileup_file=file.path(patient_pileup_file, list.files(patient_pileup_file, pattern = 'pileups.freq$'))
}
if(length(patient_pileup_file)==0) {
  stop("Usage: Input file should be one file in pileup format, or a dirctory contains at least one pileup file", call.=FALSE)
}

print(patient_pileup_file)
 

bg_rate=trained_BMER[[1]]
bg_depth=trained_BMER[[2]]
shrinkage_rate=trained_BMER[[3]]

# function to define BMER and call mutation for input file
polish.pred.3nt=function(shrinkage_rate,  #basian posterior mutation rate and depth information
                         bg_depth,
                         input_sample, #input pileup file
                         input_alpha,     #alpha = significance level, control polishing strength
                         output_dir='.'){ 
  cat(paste("Polishing",input_sample,"\n"))
  
  input_pileup=read.table(input_sample, header = TRUE, stringsAsFactors = FALSE) 
  input_pileup=input_pileup[,-dim(input_pileup)[2]]
  colnames(input_pileup)=c('CHR','POSITION','REF','DEPTH', 'A', 'C', 'G','T')
  input_pileup$REF=toupper(input_pileup$REF)
  rownames(input_pileup)=paste(input_pileup$CHR,input_pileup$POSITION,sep="_")
  input_pileup=input_pileup[order(rownames(input_pileup)),]
  
  ###remove counts of the original base
  input_rate=input_pileup[,5:8]/input_pileup$DEPTH
  input_rate[input_pileup$REF=='A','A']=0
  input_rate[input_pileup$REF=='C','C']=0
  input_rate[input_pileup$REF=='G','G']=0
  input_rate[input_pileup$REF=='T','T']=0
  
  ##function to adjust input_rate using base-coverage information in patients
  p_adj=function(p, depth, limit_n=300){
    depth_adj=depth/limit_n
    depth_adj[depth_adj>1]=1
    return(p*depth_adj**0.5)
  }
  
  input_rate_adj=apply(input_rate, 2, function(x) p_adj(x, depth = input_pileup$DEPTH, limit_n = 400))
  
  #calculate the average count using input rate and depth 
  avg_count=shrinkage_rate*bg_depth 
  
  cat('Computing site specific frequency\n')
  #define background for binomial using upper CI limit
  cat(paste0('Computing site specific upper CI limit with alpha = ', input_alpha, '\n'))
  CI_limit=qbeta(1-as.numeric(input_alpha)/2,as.matrix(avg_count+1),as.matrix(bg_depth-avg_count))
  rownames(CI_limit)=rownames(bg_depth)
  
  ##just to make sure the entries are in the same order
  I=match(rownames(input_rate), rownames(CI_limit))
  CI_limit=CI_limit[I,]
  # call mutation if observed data > background
  # also the rate should not higher than 0.3 to remove SNPs to remove SNP calls
  pred=as.logical(rowSums(input_rate_adj>CI_limit) & apply(input_rate,1,max)<0.3)
  
  
  # plot(CI_limit[,4],bg_rate[I,][,7], main='mutation rate upper CI limit
       # gets higher when background rate inflates',  ylab='background',xlab='CI')
  # plot(CI_limit[,4],bg_depth[I,][,1], main='mutation rate upper CI limit
       # gets higher when background depth is low',  ylab='background',xlab='CI')
  
  
  # polished noise data output (with mutation but not pass the filter)
  # polish.noise=cbind(input_pileup, noise.num=rowSums(input_pileup[,5:8])>0,
  #                    MajorMutFreq=apply(d1.2,1,max)/pileup.data$DEPTH,
  #                    threshold=d3[cbind(1:nrow(d3),apply(d1.2,1,which.max))])
  # ###only select mutations that did not pass threshold filter
  # polish.noise=polish.noise[rowSums(d1.2)>0 & pred==FALSE,]
  # 
  # #output the polished data
  # nt.column.loc=match(c(rbind(s,tolower(s))),colnames(pileup.data)) #NT columns, used to be hard coded  7:14
  # out.pileup=pileup.data
  # out.pileup[pred==0,nt.column.loc]=0 # polish all positions without detected mutation
  # 
  # #handling output file name and directory
  out.file.name=input_sample
  # out.file.name1=paste(out.file.name,".clean.txt",sep="")
  # out.file.name2=paste(out.file.name,".noise.txt",sep="")
  # write.table(out.pileup,file=out.file.name1,sep="\t",quote=F,row.names=F)
  # write.table(polish.noise,file=out.file.name2,sep="\t",quote=F,row.names=F)
  write.table(input_pileup[pred,],file=file.path(
    output_dir, paste0(out.file.name,'.real_mutations.txt')),sep="\t",quote=F,row.names=F)
  # return(pred)
  # 
}

# call the function to run each input sample
##############################################
n_sample=length(patient_pileup_file)
for (i in 1:n_sample) {
  polish.pred.3nt(shrinkage_rate = shrinkage_rate, 
                  bg_depth= bg_depth,   
                  input_sample = patient_pileup_file[i],  
                  input_alpha=input_alpha,  
                  output_dir='.')
}
