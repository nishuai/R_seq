library(mixtools)

find_cutoff=function(rates){
  if(sum(is.na(rates))>0){stop('rates must not have NA values')}
  bb=density(log10(as.matrix(rates)))
  ###the local minimium of mutation rate 
  slide_window=50
  I=which.min(sapply(1:(length(bb$y)-slide_window), function(x){
    sum(bb$y[x:(x+slide_window)])}))+round(slide_window/2)
  soft_cutoff=10**bb$x[I]
  
  ###fit a mixture of 2 normal distributions (2 SNP muration rates,
  ###haploid and diploid SNPs)
  ###defind the goodness of fit, the sigmal of haploid SNP should 
  ###roughly double the value of deploid SNPs
  rates=as.vector(as.matrix(rates))
  rates=rates[rates>0.2]
  if (length(rates)>50) {
    ##add 20 dummies around 1 to avoid 0 variance at mean 1
    rates=c(rates, rnorm(20, mean = 0, sd = 0.01)+1)
    mixture<-normalmixEM(rates,k=2, mu = c(0.5, 1), verb = FALSE, maxrestarts = 10000) 
    K=which.min(mixture$mu)
    # qnorm(0.0001)=3.7
    
    lower_bound=mixture$mu[K]-3.7*mixture$sigma[K]
    if (lower_bound<0) return (soft_cutoff) else  lower_bound
  } else {
    ##if no enough training data
    return(soft_cutoff)
  }
  
}
