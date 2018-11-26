###ctDNA_training uses multipe ML algorithms trying to 
###predict the esiitance of ctDNA for cancer early diagnosis
###and to selecte the most promissing features for prediction
setwd('D:/wrk/TNER-master/')
library(randomForest)
library(MASS)
source('R/get_features.R')##this function also loads find_cutoff function
source('R/log_ratio_normal.R')
###Calculate features for patient data 
###don't be confused by the varialbe name, this is also part of training data
# test_features=list()
# for (filename in list.files('test/', pattern = 'pileups.freq$', full.names = TRUE)){
#   pileup=read.table(filename, header=TRUE, as.is=TRUE)
#   if(mean(pileup$read.depth)>800)
#   test_features[[filename]]=get_features(pileup)
# }
# saveRDS(test_features,'inputdata/test_feature_list.rds')
test_features=readRDS('inputdata/test_feature_list.rds')


# tr_features=list()
# for (filename in list.files('training/', pattern = 'pileup.txt$', full.names = TRUE)){
#   pileup=read.table(filename, header=TRUE, as.is=TRUE)
#   if(mean(pileup$read.depth)>800)
#   tr_features[[filename]]=get_features(pileup)
# }
# saveRDS(tr_features,'inputdata/tr_feature_list.rds')

tr_features=readRDS('inputdata/tr_feature_list.rds')
##training random forest 
tr_features=do.call(rbind, tr_features)
test_features=do.call(rbind, test_features)
all_features=rbind(tr_features, test_features)
all_features$y=c(rep('0', nrow(tr_features)), rep('1', nrow(test_features)))
all_features$y= as.factor(all_features$y)
# model1=randomForest(y ~ . , data =all_features, ntree =10000, mtry=2)
# importance(model1, type = 2, scale = TRUE)

##some visualizations
# for (i in colnames(all_features)){
#  plot(all_features$y, all_features[,i], main=paste(i))
# }
###linear model forward selection
all_features$y=as.numeric(as.character(all_features$y))
pval=c()
for (i in 1:(length(colnames(all_features))-1)){
  f <- as.formula(paste('y ~', colnames(all_features)[i]))
  modelAllHexSubscales <- lm(f, all_features)
  pvalue=summary(modelAllHexSubscales)$coefficients[2,4]
  pval=c(pval, pvalue)
}
names(pval)=colnames(all_features)[-dim(all_features)[2]]
pval=sort(pval)
linear_features=names(pval)[1:15]

# ###pca (doesn't do the job at all)
# pr_test=t(t(all_features[,1:7])/colMeans(all_features[,1:7]))
# pr_test=data.frame(pr_test)
# pr = prcomp(pr_test)
# summary(pr,loadings=TRUE)
# pr$sdev
# features=data.frame(cbind(pr$x, y=all_features$y))
# features$y=factor(features$y)
# 
# features=data.frame(cbind(pr_test, y=all_features$y))
# features$y=factor(features$y)
# model2=randomForest(y ~ ., data =features, ntree = 20000, mtry=4)
# model3=randomForest(y ~ ., data =all_features, ntree = 20000, mtry=4)
# ####logistic regression feature forward selection 
# colnames(all_features)

#logistic regression 
# logstic_pval=c()
# for (i in 1:(length(colnames(all_features))-1)){
#   f <- as.formula(paste('y ~', colnames(all_features)[i]))
#   modelAllHexSubscales <- glm(f, all_features, family = "binomial")
#   pvalue=summary(modelAllHexSubscales)$coefficients[2,4]
#   logstic_pval=c(logstic_pval, pvalue)
# }
# names(logstic_pval)=colnames(all_features)[-length(colnames(all_features))]
# sort(logstic_pval)
# 
# plot(all_features[names(sort(logstic_pval))[1:15]])
# 

###forward selection with glm
#summary(glm(y ~ mad_mean_depth+mut_count_2+panel_allbase_mr, data =all_features, family = "binomial"))


#####################
###Naive bayesian#### 
#####################
##fitting a Gaussian distribution for each feature in 2 conditions
### and save the parameters

###name feature_fit matrix for training and test group
feature_fit_tr=sapply(1:(ncol(all_features)-1), function(x) {
  fitdistr(all_features[all_features$y==0,x],densfun = 'normal')$estimate
})
colnames(feature_fit_tr)=colnames(all_features)[1:ncol(feature_fit_tr)]

feature_fit_test=sapply(1:(ncol(all_features)-1), function(x) {
  fitdistr(all_features[all_features$y==1,x],densfun = 'normal')$estimate
})
colnames(feature_fit_test)=colnames(all_features)[1:ncol(feature_fit_test)]

###defind clusters of features using hierarchical clustering
###for each feature get the log ratio

# performance=list()
# performance_feature=list()
# for (try in 1:10000){
#   patient=c(); control=c();
#   if (try %% 100 ==0) print(try)
#   try_features=sample(1:48, 10)
#   test_features_now=test_features[,try_features]
# 
#   ###for each person get features(for),
#   ###for feature get the log ratio(sapply)
#   for (i in 1:nrow(test_features_now)){
#   dummy=test_features_now[i,]
#   patient=c(patient, mean(sapply(1:length(dummy), function(x)
#     log_ratio(dummy[x], feature_fit_test[,try_features][,x], feature_fit_tr[,try_features][,x]))))
#   }
# 
#   ##the same try_features
#   tr_features_now=tr_features[,try_features]
#   for (i in 1:nrow(tr_features_now)){
#   dummy=tr_features_now[i,]
#   control=c(control, mean(sapply(1:length(dummy), function(x)
#     log_ratio(dummy[x], feature_fit_test[,try_features][,x], feature_fit_tr[,try_features][,x]))))
#   }
# 
#   performance[[try]]=c(patient, control)
#   performance_feature[[try]]=try_features
# }
# saveRDS(performance,'inputdata/performance.rds')
# saveRDS(performance_feature,'inputdata/performance_feature.rds')

###performance contains 10000 simulations of random selection of 10 features
###and their estimated performance, performance_features contains the 
###corresponding features randomly selected
performance=readRDS('inputdata/performance.rds')
performance_features=readRDS('inputdata/performance_feature.rds')
##the log values can be added up to select the best features
##the value is confined between -4 and 4 to avoid outliers
##for determin ctDNA abundance 
aa=sapply(performance, function(x){ 
  len1=nrow(test_features)
  len2=nrow(tr_features)
  sum(ifelse(x[1:len1]>4, 4, x[1:len1]))/len1-
    sum(ifelse(x[(len1+1):(len1+len2)]<-4, -4, x[(len1+1):(len1+len2)]))/2
  })

I=aa>3.6
best_f_index=performance_features[I]
bb=unlist(best_f_index)
best_f=as.numeric(names(sort(table(bb), decreasing = TRUE))[1:15])
best_features=colnames(all_features[,best_f])

#####test the best features
###remove correlated features in linear model 
feature_matrix=all_features[,linear_features]
rm_1cor=function(feature_matrix){
  cor_matrix=abs(cor(feature_matrix))
  limit=rowSums(as.matrix(cor_matrix)>0.8)
  if(all(limit<2)) return(feature_matrix)
  del_feature=which(limit>1)[1]
  feature_matrix=feature_matrix[,-del_feature]
  rm_1cor(feature_matrix)
}
linear_features=colnames(rm_1cor(feature_matrix))
selected_features=unique(c(rev(linear_features[!linear_features %in% best_features]), rev(best_features)))
selected_features=colnames(rm_1cor(all_features[,selected_features]))

###for each person get features(for), 
###for feature get the log ratio(sapply)
test_features_now=test_features[,selected_features]
patient=list()
for (i in 1:nrow(test_features_now)){
dummy=test_features_now[i,]
patient[[i]]=sapply(1:length(dummy), function(x)
  log_ratio(dummy[x], feature_fit_test[,selected_features][,x], feature_fit_tr[,selected_features][,x]))
}
##the same try_features

tr_features_now=tr_features[,selected_features]
control=list()
for (i in 1:nrow(tr_features_now)){
dummy=tr_features_now[i,]
control[[i]]=sapply(1:length(dummy), function(x)
  log_ratio(dummy[x], feature_fit_test[,selected_features][,x], feature_fit_tr[,selected_features][,x]))
}

##some visulalization
require(ggplot2)

aa=data.frame(value=c(sapply(patient,function(x) sum(x)),
                      sapply(control,function(x) sum(x))),
              condition=c(rep('patient',length(patient)), rep('control', length(control))))
ggplot(aa, aes(x = condition, y = value))+geom_boxplot(aes(colour =condition))+
  geom_jitter(aes(color=condition), position=position_jitterdodge(dodge.width=0.001))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
        geom_hline(yintercept = 0)

##post training selection
Reduce('+', lapply(patient, function(x) x>0))
Reduce('+', lapply(control, function(x) x<0))

Reduce('+',patient)
Reduce('+',control)

sum(sapply(patient, function(x) sum(x)>1))

