library(nolabel)
library(data.table)
library(Matrix)
library(pROC)

#Hyperparameters
lambda <- 10
gamma <- .01

set.seed(123456)

#The dataset is provided in the form of (label, bag, features)
#with the label being still the binary labels.
#This is useful in case of testing the model
heart.sample <- heart.data[sample(nrow(heart.data)),] #shuffle

test.fold <- 1:(floor(nrow(heart.sample)/5))
testset   <- heart.sample[test.fold,]
trainset  <- heart.sample[-test.fold,]

N <- length(unique(trainset$bag)) #count the bags into the trainset

#Logistic regression - the labels are the binary ones, not proportions - Oracle
trainset  <- heart.sample[-test.fold,]
w.lr      <- logistic.regression(trainset,lambda)
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lr))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)


#Cast label in proportions
for (bag in unique(trainset$bag)) {
  id <- which(trainset$bag==bag)
  trainset$label[id] <- rep(mean((trainset$label[id]+1)/2), length(id))
}

#Mean Map
w.mm      <- mean.map(trainset,lambda)
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.mm))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)


#Laplacian Mean Map with similarity v^G
#Functions for building the Laplacian matrix are into laplacian.mean.map.R
L         <- laplacian(similarity="G,s", trainset, N, sigma=1)
w.lmm     <- laplacian.mean.map(trainset, lambda, gamma, L = L)
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lmm))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)


#Laplacian Mean Map with similarity v^NC
L         <- laplacian(similarity="NC", trainset, N)
w.lmm     <- laplacian.mean.map(trainset, lambda, gamma, L = L)
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lmm))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)


#Alternating Mean Map started with MM
w.amm     <- alternating.mean.map(trainset, lambda=lambda, init="1")
w.amm     <- w.amm$theta #the algorithm returns a structure that contains also the number of step until termination
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.amm))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)

#Alternating Mean Map started with LMM with similarity v^{G,s}
L         <- laplacian(similarity="G,s", trainset, N, sigma=10)
w.amm     <- alternating.mean.map(trainset, lambda=lambda, init="1", L=L, gamma=gamma, minmax = T)
w.amm     <- w.amm$theta #the algorithm returns a structure that contains also the number of step until termination
test.X    <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.amm))
test.auc  <- auc((testset$label+1)/2, c(test.pred))
print(test.auc)
