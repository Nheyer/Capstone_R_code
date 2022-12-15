##import data##
require(arm)
require(DescTools)
require(Zelig)
require(glmnet)
require(matrixcalc)
raw_data <- read.csv("/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/Original_22_People.csv", sep=",")
###   View(raw_data)


##Processing##
invariants                <- c(26,34,60,70,72)
RN_varinats               <- c(16,25,33,64)
indexes_of_best_varyables <- c(18,19,23,24,30,38,49,50,54,55,56,61,62,63,65)


Big       <-raw_data[,indexes_of_best_varyables]

### Priors Determination ###
#priors represented as beta hyper peramiters, hyper peramiters were selectedx by looking at distrabutions
Theo_people <- read.csv("/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/PriorsF.csv",sep = ",")
### View(Theo_people)

M = as.matrix(cbind(rep(1,15),Theo_people[,-c(15,16)]))
Minv = matrix.inverse(M)

#generate distrabutions of sample size 100,000
n = c()
for(i in 1:15){
  
  pi  = rbeta(100000,Theo_people$a[i],Theo_people$b[i])
  
  n <- c(n,as.numeric(log(pi/(1-pi))))
}
# convert to independent normal distrabutions
ln_beta = Matrix(n,nrow = 15, ncol = 100000,sparse = F)
mtrx = Minv %*% ln_beta

#aproximate mu and sigma^2 by the law of large numbers 
info <- function(x){ return(c(mean(x),sqrt(var(x))))}
indipendent_norm_priors = data.frame()



for(i in 1:15){
  indipendent_norm_priors = rbind(indipendent_norm_priors,info(mtrx[i,]))
}
names(indipendent_norm_priors) <- c("mu","Sd")

#define number of positive and negitive in sample
num_positive = sum(Big$ND.score == 1)
num_negitive = sum(Big$ND.score == 0)


#selection of lambda

x_vals <- Big[,names(Big) != rep("ND.score",length(names(Big)))]
  
lambda_poss = seq(0,1,0.1)
AUC_lambda = rep(NA,length(lambda_poss))

  
# AUC is not implemented in CV so do it manually 

lambda_poss = seq(0,0.1,0.1)
AUC_lambda = rep(NA,length(lambda_poss))
for (i in 1:length(lambda_poss)){
  # set the thresholds and space to hold 
  thresholds = seq(0,1,0.01)
  TPR = FPR = rep(NA,length(thresholds))
  for (k in 1:length(thresholds)) {
    num_TP = 0
    num_FP = 0
    for(c1 in 1:length(Big$ND.score)){
      ###c1 = 1
      working_model = glmnet(x_vals[-c(c1),]
                               ,Big$ND.score[-c(c1)]
                               ,lambda = lambda_poss[i]
                               ,alpha = 0.1
                               ,family = "binomial")
      
      
      beta0_to_beta_14 = c(working_model$a0,as.vector(working_model$beta))
      # get individual values for gene dosage
      X0_to_X14 = c(1,x_vals[c1,])
      # find liniar componat
      N = as.numeric(X0_to_X14) * as.numeric(as.vector(beta0_to_beta_14))
      # ecponentiate be couse logistic 
      E = exp(sum(N))
      prob = E/(1+E)
      # this is prediction we prob should have used predict.glmnet.... 
      if(prob > thresholds[k]){
        #we think they have an ND (1)
      num_TP = num_TP + Big$ND.score[c1]
      num_FP = num_FP + (1-Big$ND.score[c1])
      }# else {
      #we think they are nurotypical (0) so do nothing
      #}
    }
    TPR[k] = num_TP / num_positive
    FPR[k] = num_FP / num_negitive
  }
  AUC_lambda[i] = AUC(TPR,FPR)
}
# get value from graph 
best_lambda = 0.151
x = glmnet(x_vals,Big$ND.score,lambda = 0.151,alpha = 1,family = "binomial")
x$beta
lasso_suggests <- Big[,c(7,9,14,15)]
lasso_suggests_priors <- indipendent_norm_priors[c(7,9,14)+1,]


for(j in 1:length(thresholds)){
    num_TP = 0
    num_FP = 0
  
    for(i in 1:length(Big$ND.score)){
       model_working = bayesglm(ND.score ~ .,
                               data = Big[-c(i),],
                               family = binomial,
                               prior.mean.for.intercept = indipendent_norm_priors$mu[1],
                               prior.scale.for.intercept = indipendent_norm_priors$Sd[1],
                               prior.mean = indipendent_norm_priors$mu[2:15],
                               prior.scale = indipendent_norm_priors$Sd[2:15],
                               prior.df.for.intercept = Inf,
                               prior.df = Inf
                               )
      ND_likelyhood = predict.glm(model_working
                                  ,newdata = Big[i,]
                                  ,type="response")
      #print(ND_likelyhood)
      if(ND_likelyhood > thresholds[j]){
        if(Big$ND.score[i] == 1 ){
          num_TP = num_TP + 1
        }else{
          num_FP = num_FP + 1
        }
      }
    }
    TPR[j] = num_TP/num_positive
    FPR[j] = num_FP/num_negitive
  }
plot(FPR,TPR,main = paste("ROC curve of a model including all likely biologicaly informative allelies with an AUR of",AUC(FPR,TPR)) ,xlab = "False Positive Rate", ylab = "True Positive Rate",type = "l")


### now make it for lasso suggested peramiters ###
for(j in 1:length(thresholds)){
  num_TP = 0
  num_FP = 0
  
  for(i in 1:length(Big$ND.score)){
    model_working = bayesglm(ND.score ~ .,
                             data = lasso_suggests[-c(i),],
                             family = binomial,
                             prior.mean.for.intercept = indipendent_norm_priors$mu[1],
                             prior.scale.for.intercept = indipendent_norm_priors$Sd[1],
                             prior.mean = lasso_suggests_priors$mu,
                             prior.scale = lasso_suggests_priors$Sd,
                             prior.df.for.intercept = Inf,
                             prior.df = Inf
    )
    ND_likelyhood = predict.glm(model_working
                                ,newdata = lasso_suggests[i,]
                                ,type="response")
    #print(ND_likelyhood)
    if(ND_likelyhood > thresholds[j]){
      if(Big$ND.score[i] == 1 ){
        num_TP = num_TP + 1
      }else{
        num_FP = num_FP + 1
      }
    }
  }
  TPR[j] = num_TP/num_positive
  FPR[j] = num_FP/num_negitive
}
plot(FPR,TPR,main = paste("ROC curve of a model including only allelies suggested by lasso with an AUR of",AUC(FPR,TPR)) ,xlab = "False Positive Rate", ylab = "True Positive Rate",type = "l")



lasso_suggested_model = bayesglm(ND.score ~ .,
                                 data = lasso_suggests,
                                 family = binomial,
                                 prior.mean.for.intercept = indipendent_norm_priors$mu[1],
                                 prior.scale.for.intercept = indipendent_norm_priors$Sd[1],
                                 prior.mean = lasso_suggests_priors$mu,
                                 prior.scale = lasso_suggests_priors$Sd,
                                 prior.df.for.intercept = Inf,
                                 prior.df = Inf,
                                 maxit = 10000
                                 )
full_model = bayesglm(ND.score ~ .,
                      data = Big[-c(i),],
                      family = binomial,
                      prior.mean.for.intercept = indipendent_norm_priors$mu[1],
                      prior.scale.for.intercept = indipendent_norm_priors$Sd[1],
                      prior.mean = indipendent_norm_priors$mu[2:15],
                      prior.scale = indipendent_norm_priors$Sd[2:15],
                      prior.df.for.intercept = Inf,
                      prior.df = Inf,
                      maxit = 10000
                      )

lasso_suggested_model_sims = sim(lasso_suggested_model,n.sims = 100000)                      
lasso_suggested_model_sims_coef = coef(lasso_suggested_model_sims)
Summery_lasso_suggested_model = matrix(nrow = length(lasso_suggested_model_sims_coef[1,]),ncol = 3)
for (i in 1:length(lasso_suggested_model_sims_coef[1,])) {
  Summery_lasso_suggested_model[i,1] = mean(lasso_suggested_model_sims_coef[,i])
  Summery_lasso_suggested_model[i,2] = median(lasso_suggested_model_sims_coef[,i])
  Summery_lasso_suggested_model[i,3] = sqrt(var(lasso_suggested_model_sims_coef[,i]))
}
row.names(Summery_lasso_suggested_model) <- names(lasso_suggested_model_sims_coef)
colnames(Summery_lasso_suggested_model) <- c("Mean","Median","Varience")
write.csv2(Summery_lasso_suggested_model,file = "/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/lasso_postirior_beta.txt",eol = "\\\\ \n")



full_model_sims = sim(full_model, n.sims = 100000)
full_model_sims_coef = coef(full_model_sims)
Summary_of_full_model = matrix(nrow = length(full_model_sims_coef[1,]),ncol = 3)
for (i in 1:length(full_model_sims_coef[1,])) {
  Summary_of_full_model[i,1] = mean(full_model_sims_coef[,i])
  Summary_of_full_model[i,2] = median(full_model_sims_coef[,i])
  Summary_of_full_model[i,3] = sqrt(var(full_model_sims_coef[,i]))

}
row.names(Summary_of_full_model) <- names(full_model_sims_coef[1,])
colnames(Summary_of_full_model) <- c("Mean","Median","Varience")
Summary_of_full_model = data.frame(Summary_of_full_model)
View(Summary_of_full_model)
write.csv2(Summary_of_full_model,file="/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/postirior_beta.txt",sep="&",eol = "\\\\ \n")
Summary_of_full_model_odds_ratio = apply(Summary_of_full_model,MARGIN = 2,FUN=exp)
write.csv2(Summary_of_full_model_odds_ratio[,1],file = "/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/odds_ratios.csv",eol = "\\\\ \n")
View(Summary_of_full_model_odds_ratio)


summary(full_model_sims_coef)














###### SCRATCH for parallel lambda #####
Get_AUC_lasso<-function(l,ac=0.01){
  require(glmnet)
  require(DescTools)
  raw_data <- read.csv("/home/nick/Desktop/CSUMB/Classes/Spring_2021/Data/Original_22_People.csv", sep=",")
  Big <- raw_data[,c(18,19,23,24,30,38,49,50,54,55,56,61,62,63,65)]
  #define info
x_vals <- Big[,names(Big) != rep("ND.score",length(names(Big)))]
num_positive = sum(Big$ND.score == 1)
num_negitive = sum(Big$ND.score == 0)
thresholds = seq(0,1,ac)
TPR = FPR = rep(NA,length(thresholds))



for (k in 1:length(thresholds)) {
  num_TP = 0
  num_FP = 0
  for(c1 in 1:length(Big$ND.score)){
    ###c1 = 1
    working_model = glmnet(x_vals[-c(c1),]
                           ,Big$ND.score[-c(c1)]
                           ,lambda = l
                           ,alpha = 1
                           ,family = "binomial")
    
    
    beta0_to_beta_14 = c(working_model$a0,as.vector(working_model$beta)) 
    X0_to_X14 = c(1,x_vals[c1,])
    N = as.numeric(X0_to_X14) * as.numeric(as.vector(beta0_to_beta_14))
    E = exp(sum(N))
    prob = E/(1+E)
    if(prob > thresholds[k]){
      #we think they have an ND (1)
      num_TP = num_TP + Big$ND.score[c1]
      num_FP = num_FP + (1-Big$ND.score[c1])
    }# else {
    #we think they are nurotypical (0) so do nothing
    #}
  }
  TPR[k] = num_TP / num_positive
  FPR[k] = num_FP / num_negitive
}
return( AUC(TPR,FPR))}

Get_AUC_lasso(l=0.2,ac=0.00001)



lambda_poss = seq(0.1499,0.1598,0.001)
library(doParallel)
library(foreach)
myCluster <- makeCluster(10, type = "PSOCK")
registerDoParallel(myCluster)
AUC_lam = foreach(i = 1:length(lambda_poss), .combine = 'c')%dopar%{Get_AUC_lasso(lambda_poss[i],ac=0.00001)}
stopCluster(myCluster)
plot(main="Predictive Ability of Lasso Models with a Given Lambda Perameter",x=lambda_poss[1:length(AUC_lam)],y=AUC_lam, ylab = "Area under the Receiver Operating Curve", xlab = "Lambda penalty")
length(lambda_poss)
