#################################################################
# File name: TNER_training.R
# Render training sets to have the same format for easier bg error rate estimation
# Author: Ni Shuai nishuai@yahoo.com
# Ri&De YuceBio
# Input: directory of all training pileup files, reference sequnce used in the panel, 
# bed file defining coordinates of the panel, output directory (cwd by default)
##training sample should be in one directory with file extension 'pileups.freq'
# impvoded over original TNER.training to release RAM after each file is processed
#################################################################
####get input parameters
cat('Start computing\n')
args = commandArgs(trailingOnly=TRUE)
 # args=c('training/','ref/Yuce1plus.fa','ref/YConePlus_v1.3.bed', 0.001,'.')
if (length(args)<3) {
  stop("Usage: Rscript TNER_training.R training_dir ref_seq bed_file input_alpha output_dir", call.=FALSE)
} else if (length(args)>5) {
  stop("Usage: Rscript TNER_training.R training_dir ref_seq bed_file input_alpha output_dir", call.=FALSE)
}
training_dir = output_dir= input_alpha=NULL
training_dir = args[1]  
input_alpha= as.numeric(args[4])
output_dir= args[5]
if (is.null(output_dir)) output_dir='.'  
if (is.na(output_dir)) output_dir='.'  
training_files=list.files(path=training_dir,pattern="pileups.freq$") # list all healthy subject pileup data files; require at least two or more normal samples

if (is.null(input_alpha)) input_alpha=0.001

###generate a master atgc file (each training may have missing positions)
ref_seq=scan(args[2], what='character')
ref_seq=ref_seq[seq(2, length(ref_seq), 2)]
ref_seq=paste0(ref_seq, collapse = '')
ref_seq=toupper(ref_seq)

bed_file=read.table(args[3], stringsAsFactors = FALSE, blank.lines.skip = TRUE)
all_positions=unlist(apply(bed_file,1, function(x) {(as.numeric(x[2])+1):x[3]}))
chr=rep(bed_file$V1, c(bed_file$V3 - bed_file$V2))

###encode chromosome name and position into one integer to avoid string operations
chr_id=substr(chr, 4, 7)
chr_id[chr_id=='X']=23
chr_id=as.integer(chr_id)

master_acgt=data.frame(chr=chr, position=all_positions, ref_base=unlist(strsplit(ref_seq, split = '')))
master_acgt$id=master_acgt$position+chr_id*1e10
colnames(master_acgt)=c('CHR','POSITON','REF','id')
rm('ref_seq','all_positions','bed_file')
####read into each training acgt data and render for missing bases
depth=matrix(0, nrow=0,ncol=dim(master_acgt)[1])
acgt_bg=matrix(0, nrow=dim(master_acgt)[1], ncol = 4)
acgt_bg_mean=matrix(0, nrow=dim(master_acgt)[1], ncol = 4)
acgt_bg_times=matrix(0, nrow=dim(master_acgt)[1], ncol = 4)
acgt_tnt=matrix(0, nrow=dim(master_acgt)[1], ncol = 4)
rownames(acgt_bg)=master_acgt$id
rownames(acgt_bg_mean)=master_acgt$id
rownames(acgt_tnt)=master_acgt$id

for (i in 1:length(training_files)){
  cat(paste('processing ', training_files[i], '\n'))
  training=read.table(file.path(training_dir, training_files[i]), as.is = TRUE, header = TRUE)
  training=training[,-dim(training)[2]]
  chr=substr(training$chr, 4, 7)
  chr[chr=='X']=23
  chr=as.integer(chr)
  ###create a unified identifier
  training$id=training$n_base+chr*1e10
  bb=merge(master_acgt, training, by ='id', all=TRUE)
  ##do not remove NA when calculating depth
  depth=rbind(depth,bb$read.depth)
  acgt=bb[,c('A','C','G','T')]
  acgt[master_acgt$REF=='A','A']=0
  acgt[master_acgt$REF=='C','C']=0
  acgt[master_acgt$REF=='G','G']=0
  acgt[master_acgt$REF=='T','T']=0
  acgt=acgt/bb$read.depth
  acgt[is.na(acgt)]=0
  acgt_bg=pmax(acgt_bg,as.matrix(acgt), na.rm = TRUE)
  
  acgt_bg_mean=acgt_bg_mean+as.matrix(acgt)
  acgt_bg_times=acgt_bg_times+as.matrix(acgt>0.05)
  tnt_rate=acgt
  tnt_rate[tnt_rate>0.271]=0
  acgt_tnt=acgt_tnt+tnt_rate
  # All_training[[filename]]=bb
  rm(bb)
  rm(training)
}
cat('removing catch')
acgt_bg_mean=acgt_bg_mean/i
print(head(acgt_bg))
print(head(acgt_bg_mean))
acgt_bg[acgt_bg_times<10]=acgt_bg_mean[acgt_bg_times<10]
acgt_bg=acgt_bg*5
print(head(acgt_bg))
acgt_bg[acgt_bg>1]=1
print(head(acgt_bg))
acgt_tnt=acgt_tnt/i

##for those pmax is large but only present in less than 2 samples, 
##change the rate to mean rate rather than pmax*5

depth=colMeans(depth, na.rm = TRUE)
depth[is.na(depth)]=0

colnames(acgt_tnt)=c('A','C','G','T')
colnames(acgt_bg)=c('A','C','G','T')
acgt_bg=cbind(master_acgt[,1:3],acgt_bg)
acgt_tnt=cbind(master_acgt[,1:3],acgt_tnt)

bg_depth=depth
cutoff_rate=acgt_tnt
bg_rate=acgt_bg


###############################################################
#create a Tri-nucleotide (tnt) code for each position
ref=bg_rate$REF
ref1=ref[1:(length(ref)-2)]
ref2=ref[2:(length(ref)-1)]
ref3=ref[3:length(ref)]
tnt=paste0(ref1, ref2, ref3)
tnt=c(NA, tnt, NA)##the first and last nucletide do not have a corresponding 3-mer
cat('Generating 3 mers \n')
# remove the tnt code for the position at sequence gaps
x=bg_rate$POSITION
###remove the 3-mer end and start with the 2 GAP-nucleotides
I=which(diff(x)!=1)
tnt[c(I, I+1)]=NA
## Empirical Bayesian - a shrinkage approach
##################################################
#calculate the tri-nucleotide average and SD
cat('Generating 3-mer global frequencies\n')

## mutation rate estimation
###TNT rate can be estimated sample-by-sample if sample size gets larger than ~50
###to count for inflated MR by a single SNP in one sample
tnt_rate=aggregate(cutoff_rate[,4:7],list(TNT=tnt),function(x) mean(x,na.rm=T))
row.names(tnt_rate)=tnt_rate$TNT
## quality check, helps to define the cutoff mutation rate for SNP exclusion
# library(pheatmap)
# pheatmap(tnt_rate[,2:5], main='3-mer specific mutation rate in 10 training samples')
# plot(density(as.matrix(tnt_rate[,2:5])))
# match the average tnt rate to the input data
# enumerate Background rate for each base in panel to match with local
tnt_rate_site=tnt_rate[match(tnt,tnt_rate$TNT), -1]## To remove the first column TNT
##deal with unavaliable rates due to disconnected 3-mers at breakpoints
tnt_rate_site[is.na(tnt_rate_site)] = median(as.matrix(tnt_rate_site), na.rm = TRUE)

#shrinkage weight parameter 
###The cause of high site specific backgroupnd rate of any reason will 
###reduce w, thus weights down the 3-mer mutation rate in calculating
###posterior mutation rate
w=tnt_rate_site/(tnt_rate_site+bg_rate[,4:7]) 
##Only when 3-mer and site specific muation rates are both zero will make w NA
###actually doesn't matter what w is
w[is.na(w)]=1

# average mutation error rate
cat('Calibrating global BMER with local BMER using shrinkage estimator\n')
shrinkage_rate=w*tnt_rate_site+(1-w)*bg_rate[,4:7]
rownames(shrinkage_rate)=master_acgt$id

avg_count=shrinkage_rate*bg_depth
  cat('Computing site specific frequency\n')
  #define background for binomial using upper CI limit
  cat(paste0('Computing site specific upper CI limit with alpha = ', 0.01, '\n'))
  CI_UpperLimit=qbeta(1-as.numeric(0.01)/2,as.matrix(avg_count+1),as.matrix(bg_depth-avg_count))
  cat(paste0('Computing site specific lower CI limit with alpha = ',0.01, '\n'))
  CI_LowerLimit=qbeta(as.numeric(0.01)/2,as.matrix(avg_count),as.matrix(bg_depth-avg_count)+1)
  ###The lower CI is set to 0 if lower than 0.000001
  CI_LowerLimit[CI_LowerLimit<0.00001]=0

  cat(paste0('Computing site specific upper CI limit with alpha = ', 0.001, '\n'))
  CI_UpperLimit1=qbeta(1-as.numeric(0.001)/2,as.matrix(avg_count+1),as.matrix(bg_depth-avg_count))
  cat(paste0('Computing site specific lower CI limit with alpha = ',0.001, '\n'))
  CI_LowerLimit1=qbeta(as.numeric(0.001)/2,as.matrix(avg_count),as.matrix(bg_depth-avg_count)+1)
  CI_LowerLimit[CI_LowerLimit<0.00001]=0
  
saveRDS(list(bg_rate, bg_depth, shrinkage_rate, CI_UpperLimit, CI_LowerLimit,
             CI_UpperLimit1, CI_LowerLimit1), file.path(output_dir, 'trained_BMER_pmax.rds'))

# end of processing input
