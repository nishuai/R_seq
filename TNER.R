###############################################################
# File name: TNER.R
# detect ctDNA mutation from background of mutational signatures  
# Author: Ni Shuai nishuai@yahoo.com
# Inspired by: shibing Deng, Tao Xie (xietao2000@Pfizer.com)
# YuceBio R&D department
# Input: test_sample file in  atgc format, All_training file,
# output directory (gwd by default)
#################################################################

args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input
source('R/depth_p_adj.R')
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
  patient_path=patient_pileup_file
  patient_pileup_file=file.path(patient_path, list.files(patient_path, pattern = 'pileups.freq$'))
    if (length(patient_pileup_file)==0){
      patient_pileup_file=file.path(patient_path, list.files(patient_path, pattern = 'umi.pileup.txt$'))
	}
}
if(length(patient_pileup_file)==0) {
  stop("Usage: Input file should be one file in pileup format, or a dirctory contains at least one pileup file", call.=FALSE)
}

# function to define BMER and call mutation for input file
polish.pred.3nt=function(trained_BMER,  #basian posterior mutation rate and depth information and more
                         input_sample, #input pileup file
                         input_alpha,     #alpha = significance level, control polishing strength
                         output_dir='.'){
  cat(paste("Polishing",input_sample,"\n"))
  ##get information form model
  bg_rate=trained_BMER[[1]]
  bg_depth=trained_BMER[[2]]
  shrinkage_rate=trained_BMER[[3]]
  shrinkage_rate[is.na(shrinkage_rate)]= median(as.matrix(shrinkage_rate), na.rm = TRUE)
  
  ##get information form input pileup
  input_pileup=read.table(input_sample, header = TRUE, stringsAsFactors = FALSE) 
  input_pileup=input_pileup[,-dim(input_pileup)[2]]
  colnames(input_pileup)=c('CHR','POSITION','REF','DEPTH', 'A', 'C', 'G','T')
  input_pileup$REF=toupper(input_pileup$REF)
  chr_id=substr(input_pileup$CHR, 4, 7)
  chr_id[chr_id=='X']=23
  chr_id=as.integer(chr_id)
  rownames(input_pileup)=input_pileup$POSITION+chr_id*1e10
  
  ###remove counts of the original base
  input_rate=input_pileup[,5:8]/input_pileup$DEPTH
  input_rate[input_pileup$REF=='A','A']=0
  input_rate[input_pileup$REF=='C','C']=0
  input_rate[input_pileup$REF=='G','G']=0
  input_rate[input_pileup$REF=='T','T']=0
  Major_freq=apply(input_rate, 1, max)
  input_pileup=cbind(input_pileup, Major_freq)
  
  
  ###adjusting base frequency with depth lower than 500
  input_rate_adj=apply(input_rate, 2, function(x) depth_p_adj(x, depth = input_pileup$DEPTH, limit_n = 20))
  
  #calculate the average count using input rate and depth 
  ###test
  # I=match(rownames(input_pileup), rownames(shrinkage_rate))
  # shrinkage_rate=shrinkage_rate[I,]
  # bg_depth=bg_depth[I,]
  ###
  
  if (input_alpha==0.001){
    cat('Read model upper and lower CI limit\n')
    CI_UpperLimit=trained_BMER[[4]]
    CI_LowerLimit=trained_BMER[[5]]
  } else if (input_alpha==0.01){
    cat('Read model upper and lower CI limit\n')
    CI_UpperLimit=trained_BMER[[6]]
    CI_LowerLimit=trained_BMER[[7]]
  } else{
    cat('Computing site specific frequency\n')
    #define background for binomial using upper CI limit
    cat(paste0('Computing site specific upper CI limit with alpha = ', input_alpha, '\n'))
    # final_depth=ifelse(bg_depth>input_pileup$DEPTH, bg_depth, input_pileup$DEPTH)
    # CI_UpperLimit=qbeta(1-as.numeric(input_alpha)/2,as.matrix(avg_count+1),as.matrix(bg_depth-avg_count))
    # CI_UpperLimit2=qbeta(1-as.numeric(input_alpha)/2,as.matrix(avg_count+1),as.matrix(final_depth-avg_count))
    avg_count=shrinkage_rate*bg_depth 
    CI_UpperLimit=qbeta(1-as.numeric(input_alpha)/2,as.matrix(avg_count+1),as.matrix(bg_depth-avg_count))
    cat(paste0('Computing site specific lower CI limit with alpha = ', input_alpha, '\n'))
    CI_LowerLimit=qbeta(as.numeric(input_alpha)/2,as.matrix(avg_count),as.matrix(bg_depth-avg_count)+1)
    ###The lower CI is set to 0 if lower than 0.000001
    CI_LowerLimit[CI_LowerLimit<0.00001]=0
    rownames(CI_UpperLimit)=rownames(bg_depth)
    rownames(CI_LowerLimit)=rownames(bg_depth)
  }
  
  ##just to make sure the entries are in the same order
  I=match(rownames(input_rate), rownames(CI_UpperLimit))
  CI_UpperLimit=CI_UpperLimit[I,]
  CI_LowerLimit=CI_LowerLimit[I,]
  input_pileup=cbind(input_pileup, CI_UpperLimit)
  input_pileup=cbind(input_pileup, CI_LowerLimit)
  input_pileup=cbind(input_pileup, bg_rate[I,4:7])
  colnames(input_pileup)[10:21]=c(paste0('UpperCI_', c('A','C','G','T')),
                                  paste0('LowerCI_', c('A','C','G','T')),
                                  paste0('bg_', c('A','C','G','T')))
  
  ##add hotspot information
  if ('Hotspot' %in% colnames(bg_rate)){
    Hotspot_VAR=bg_rate$Hotspot[I]
    input_pileup=cbind(input_pileup, Hotspot_VAR)
  }
  row_which_max=as.vector(apply(input_rate_adj, 1, which.max))
  VAR=c('A','C','G','T')[row_which_max]
  input_pileup=cbind(input_pileup[,c(1, 2, 2, 3)], VAR, input_pileup[,4:ncol(input_pileup)])
  
  # call mutation if observed data > background
  # also the rate should not higher than the trained cutoff value to remove SNPs to remove SNP calls
  # the trained cutoff value is the lower CI confidence for haploid SNP in training patients
  pred_CI_upper=as.logical(rowSums(input_rate_adj>CI_UpperLimit)) 
  row_max=input_rate_adj[cbind(1:nrow(input_rate_adj), row_which_max)]
  row_max_CI_lower=CI_LowerLimit[cbind(1:nrow(input_rate_adj), row_which_max)]
  
  min_pred_CI_upper=as.vector(apply(CI_UpperLimit, 1, min))
  ###this is to avoid cases when SNP mutations meet real mutations
  ###but to keep a normal standard level of AF threshhold
  
   pred_CI_lower=row_max_CI_lower-row_max>0.1 & row_max>min_pred_CI_upper 
   pred=pred_CI_upper | pred_CI_lower 
   #pred=pred_CI_upper 
   pred=pred & row_max>0.05 
  # pred=pred_CI_upper & row_max	> 0.05
  # pred=pred_CI_upper 
  ###the minimium detection rate is to be discussed
  # pred=pred & apply(input_rate,1,max)<0.271
  
  # plot(CI_UpperLimit[,4],bg_rate[I,][,7], main='mutation rate upper CI limit
       # gets higher when background rate inflates',  ylab='background',xlab='CI')
  # plot(CI_UpperLimit[,4],bg_depth[I,][,1], main='mutation rate upper CI limit
       # gets higher when background depth is low',  ylab='background',xlab='CI')
  
  # polished noise data output (with mutation but not pass the filter)
  ###only select mutations that did not pass threshold filter
  polish_noise=input_pileup[input_pileup$Major_freq>0 & pred==FALSE,]
  # #output the polished data
  # nt.column.loc=match(c(rbind(s,tolower(s))),colnames(pileup.data)) #NT columns, used to be hard coded  7:14
  # out.pileup=pileup.data
  # out.pileup[pred==0,nt.column.loc]=0 # polish all positions without detected mutation
  # 
  # #handling output file name and directory
  output_name=gsub(pattern = "\\.txt$", "",input_sample)
  #out.file.name1=paste(out.file.name,".clean.txt",sep="")
  output_name2=paste(output_name,".noise.txt",sep="")
  write.table(polish_noise,file=output_name2, sep="\t", quote=F, row.names=F)
  #write.table(polish.noise,file=out.file.name2,sep="\t",quote=F,row.names=F)
  real_mu_file=input_pileup[pred,]
  write.table(real_mu_file,file=file.path(
    output_dir, paste0(output_name,'.real_mutations.txt')),sep="\t",quote=F,row.names=F)
  
  ##write a table in vcf file format
  # if (dim(real_mu_file)[1]>0){
  #   mu_vcf=cbind(real_mu_file[,c(1,2)], ID=rownames(real_mu_file),real_mu_file[,c(3,4)], QUAL=30, FILTER='PASS', INFO='0')
  # } else{
  #   mu_vcf=as.data.frame(matrix(ncol = 8, nrow = 0))
  # }
  # colnames(mu_vcf)=c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER', 'INFO')
  # write.table(mu_vcf,file=file.path(output_dir, paste0(output_name,'.real_mutations.vcf')),sep="\t",quote=F,row.names=F)
  # return(pred)
  # 
  #write a table in vcf file format
  if (dim(real_mu_file)[1]>0){
    mu_tmp_cal=real_mu_file[,c(1,2, 11, 4,5)]
  } else{
    mu_tmp_cal=as.data.frame(matrix(ncol = 5, nrow = 0))
  }
  colnames(mu_tmp_cal)=c('Chr','Pos','Vaf','Ref','Alt')
  write.table(mu_tmp_cal,file=file.path(output_dir, paste0(output_name,'.real_mutations.vaf.xls')),sep="\t",quote=F,row.names=F)
  return(pred)

}

# call the function to run each input sample
##############################################
n_sample=length(patient_pileup_file)
for (i in 1:n_sample) {
  polish.pred.3nt(trained_BMER = trained_BMER, 
                  input_sample = patient_pileup_file[i],  
                  input_alpha=input_alpha,  
                  output_dir='.')
}
