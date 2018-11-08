#################################################################
# File name: Training.R
# Render training sets to have the same format for easier bg error rate estimation
# Author: Ni Shuai nishuai@yahoo.com
# YuceBio R&D department
# Input: directory of training pileup files, reference sequnce used in the panel, 
# bed file defining coordinates of the panel, output directory (cwd by default)
#################################################################

####get input parameters
print('Start computing')
args = commandArgs(trailingOnly=TRUE)
 # args=c('training/','ref_seqs_panel.txt','Target_v2.call.bed','.')
if (length(args)<1) {
  stop("Usage: Rscript Training.R training_dir ref_seq bed_file output_dir", call.=FALSE)
} else if (length(args)>4) {
  stop("Usage: Rscript Training.R training_dir ref_seq bed_file output_dir", call.=FALSE)
}

training_dir =NULL; output_dir=NULL
training_dir = args[1]  
output_dir= args[4]
if (is.null(output_dir)) output_dir='.'  
training_files=list.files(path=training_dir,pattern="pileup.txt$") # list all healthy subject pileup data files; require at least two or more normal samples

###generate a master atgc file (each training may have missing positions)

ref_seq=scan(args[2], what='character')
ref_seq=ref_seq[seq(2, length(ref_seq), 2)]
ref_seq=paste0(ref_seq, collapse = '')


bed_file=read.table(args[3], stringsAsFactors = FALSE)
all_positions=unlist(apply(bed_file,1, function(x) {(as.numeric(x[2])+1):x[3]}))
chr=rep(bed_file$V1, c(bed_file$V3 - bed_file$V2))

master_acgt=data.frame(chr=chr, position=all_positions, ref_base=unlist(strsplit(ref_seq, split = '')))
master_acgt$id=paste(master_acgt$chr, master_acgt$position, sep='_')
# saveRDS( master_acgt, file.path(output_dir,'master_acgt.rds'))

####read into each training acgt data and render for missing bases
All_training=list()
for (filename in training_files){
  print(paste('processing ', filename))
  
  training=read.table(file.path(training_dir, filename), stringsAsFactors = FALSE, header = TRUE)
  training=training[,-dim(training)[2]]
  ###create a unified identifier
  training$id=paste(training$chr,training$n_base, sep='_')
  bb=merge(master_acgt, training, by='id', all=TRUE)
  row.names(bb)=bb$id
  bb=bb[,c(2,3,4,9:12)]
  colnames(bb)=c('CHR','POSITION','REF','A','C','G','T')
  All_training[[filename]]=bb
}

# saveRDS(All_training, 'All_training.rds')

####define funtions to process the input training set
get_nt_rate=function(All_training_list){
  nt_matrix=lapply(All_training_list, function(x) {a=x[,4:7]; a[is.na(a)]=0; return (a)})
  All_nt_matrix=Reduce('+', nt_matrix)
  All_nt_rate=All_nt_matrix/rowSums(All_nt_matrix)
  hs_bg_ave=cbind(All_training[[1]][,1:3],All_nt_rate)
  hs_bg_ave$REF=toupper(hs_bg_ave$REF)
  hs_bg_ave=hs_bg_ave[order(rownames(hs_bg_ave)),]
  
  ###the original base is not counted as mutation
  hs_bg_ave[hs_bg_ave$REF=='A','A']=0
  hs_bg_ave[hs_bg_ave$REF=='C','C']=0
  hs_bg_ave[hs_bg_ave$REF=='G','G']=0
  hs_bg_ave[hs_bg_ave$REF=='T','T']=0
  hs_bg_ave[is.na(hs_bg_ave)]=0
  return(hs_bg_ave)
}

get_depth=function(All_training_list){
  nt_matrix=lapply(All_training_list, function(x) {a=x[,4:7]; a[is.na(a)]=0; return (a)})
  All_nt_matrix=Reduce('+', nt_matrix)
  All_nt_matrix=All_nt_matrix/length(All_training_list)
  depth_count=matrix(rep(rowSums(All_nt_matrix), 4), ncol=4)
  rownames(depth_count)= rownames(All_nt_matrix)
  colnames(depth_count)=colnames(All_nt_matrix)
  depth_count=depth_count[order(rownames(depth_count)),]
  return(depth_count)
}

bg_depth=get_depth(All_training_list = All_training)
bg_rate=get_nt_rate(All_training_list = All_training)

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
cat('Generating 3-mer global frequenties\n')

## Sites with mutation rate larger than 0.001 (or even 0.0001) do not follow the 
## assumed pattern of distribution, therefor is excluded from the global 3-mer
## mutation rate estimation
###TNT rate should be estimated in another way if sample size gets larger than ~50
###to count for inflated MR by a single SNP in one sample
tnt_rate=aggregate(bg_rate[,4:7],list(TNT=tnt),function(x) mean(x[x<0.01],na.rm=T))
# row.names(tnt_rate)=tnt_rate$TNT
## quality check, helps to define the cutoff mutation rate for SNP exclusion
# library(pheatmap)
# pheatmap(tnt_rate[2:5], main='3-mer specific mutation rate in 10 training samples')

# match the average tnt rate to the input data
# enumerate Background rate for each base in panel to match with local
tnt_rate_site=tnt_rate[match(tnt,tnt_rate$TNT), -1]## To remove the first column TNT
rownames(tnt_rate_site)=rownames(bg_rate)
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

saveRDS(list(bg_rate, bg_depth, shrinkage_rate), file.path(output_dir, 'trained_BMER.rds'))


# end of processing input
