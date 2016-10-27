
codon_usage_bias=function(genome_seq, weighted=TRUE, human=TRUE){
    sst <- strsplit(genome_seq, "")[[1]]
    codons=paste0(sst[c(TRUE, FALSE,FALSE)],sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
    if (weighted==TRUE){
        codons=data.frame(codons, counts=1)
        cov_genomic=readRDS('inputdata/coverage_genomic.RDS')
        tail(cov_genomic)
        for ( i in 1:nrow(codons)){
            codons$counts[i]=cov_genomic$hour18[i*3]
        }
    } else(
        codons=data.frame(codons, counts=1)
    )
    
    codons=aggregate(codons$counts, by=list(codons$codons), sum)
    colnames(codons)=c('codon','counts')
    codons2aminoacid=read.table('inputdata/codons2aminoacid.txt', sep=',', stringsAsFactors = FALSE)
    colnames(codons2aminoacid)=c('name','Acro','codon')
    codons=merge(condons2aminoacid, codons, by='codon')
    codons=codons[order(codons$Acro),]
    sum_codon_usage=aggregate(codons$count, by=list(codons$Acro), sum)
    sum_codon_usage=rep(sum_codon_usage$x, rle(as.character(codons$Acro))$lengths)
    codons$Bias=codons$counts/sum_codon_usage
    codons$Bias=round(codons$Bias, 3)
    colnames(codons)=c('Codons','Name','Symbol','Counts','Bias')
    codons=codons[order(codons$Symbol, codons$Bias, decreasing = TRUE),]
    if (human==TRUE){
        human_codon_bias=read.csv('inputdata/codon_bias_human.txt', stringsAsFactors = FALSE)
        colnames(human_codon_bias)=c('Codons','Bias_H')
        human_codon_bias$Codons=gsub('U','T',human_codon_bias$Codons)
        codons=merge(human_codon_bias, codons, by='Codons')
        codons$Bias_H=codons$Bias_H/100
        codons$Bias_H=round(codons$Bias_H, 3)
        codons$diff_bias=codons$Bias-codons$Bias_H
        codons=codons[order(codons$Counts, decreasing =TRUE),]
        codons=codons[,c(1,3,4,5,6,2,7)]
        colnames(codons)[5]='Bias_V'
        row.names(codons)=NULL
    }
    codons
}
