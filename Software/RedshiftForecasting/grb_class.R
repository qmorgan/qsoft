#############################
## classification routines ##
#############################

#### CART ####
# Author: Joey Richards
# Description: Wraps around rpart to do CART classification.  This version of
# CART outputs a posterior probability based upon the relative fraction of
# each class that appears in each leaf. 
#
# Notes: 
# Having trouble making the "loss" keyword work; however, using the "prior" 
# keyword allows us to penalize low-as-high.  Something is screwed up with the 
# rpart loss functionality.
#
# ToDo:
# Modify code to take into account the "weights" option; should be similar
# functionality to the "prior" or "loss" functions.  Try to find the link
# between these!  Weights (or prior) is the only thing to tune here.

rpart.cv = function(x,y,nfolds=5,method="gini",loss=NULL,prior=NULL,seed=sample(1:10^5,1)){
  require(rpart)
  set.seed(seed)

  n = length(y)
  p = length(table(y))
  folds = sample(1:nfolds,n,replace=TRUE)
  predictions = matrix(0,nrow=n,ncol=p)

  # default loss function
  if(is.null(loss)){
    loss = matrix(1,p,p)
    diag(loss) = rep(0,p)}

  # default prior: prop. to observed class rates
  if(is.null(prior)){
    prior = table(y)/n
  }

  prior = prior / sum(prior)
    
  for(ii in 1:nfolds){
    print(paste("fold",ii,"of",nfolds))
    leaveout = which(folds==ii)
    # fit tree
    y.in = y[-leaveout]
    tree.fit = rpart(y.in~.,data=data.frame(x[-leaveout,]),parms=list(split=method,loss=loss,prior=prior),control=rpart.control(minsplit=2,minbucket=4,cp=.001,xval=10))
    if(dim(tree.fit$cptable)[2] < 4) {
      print("tree is a stump, skip pruning")
      tree.prune  = tree.fit }
    else {
      tree.prune = prune(tree.fit,cp=tree.fit$cptable[which.min(tree.fit$cptable[,"xerror"]),"CP"]) }
    predictions[leaveout,] = predict(tree.prune,newdata=data.frame(x[leaveout,]),type="prob")
  }
  pred = levels(y)[apply(predictions,1,which.max)]
  pred = factor(pred,levels=levels(y))
  confmat = table(pred,y)
  err.rate = 1-sum(diag(confmat))/n
  return(list(predclass=pred,predprob=predictions,confmat=confmat,err.rate=err.rate))
}


#### Random Forest (party) ####
# Author: Joey Richards
# Description: wraps around the "party" version of random forest which can 
# properly deal with missing data.  
#
# ToDo: 
# Random forest should work better than CART - just need to figure out how to 
# properly take weights into account.  Here there are two things to tune - 
# the weights and the total number of allowed trees.

rfc.cv = function(x,y,nfolds=5,testset=NULL,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
  # don't train on any of the data in testset
  # this is to use in the hierarchical classifier
  require(party)
  set.seed(seed)
  
  n = length(y)
  p = length(table(y))
  folds = sample(1:nfolds,n,replace=TRUE)
  predictions = matrix(0,nrow=n,ncol=p)

  if(is.null(mtry)){
    mtry = ceiling(sqrt(dim(x)[2]))
  }

  for(ii in 1:nfolds){
    print(paste("fold",ii,"of",nfolds))
    leaveout = which(folds==ii)
    train = cbind(y[-union(leaveout,testset)],x[-union(leaveout,testset),])
    test = cbind(y[leaveout],x[leaveout,])
    names(train)[1] = names(test)[1] = "y"
    rf.tmp = cforest(y~.,data=train,weights=weights[-leaveout],controls=cforest_classical(mtry=mtry,ntree=n.trees))
    predictions[leaveout,] = matrix(unlist(treeresponse(rf.tmp,newdata=test)),length(leaveout),p,byrow=T)
  }
  pred = levels(y)[apply(predictions,1,which.max)]
  pred = factor(pred,levels=levels(y))
  confmat = table(pred,y)
  err.rate = 1-sum(diag(confmat))/n

  return(list(predclass=pred,predprob=predictions,confmat=confmat,err.rate=err.rate))
}

# GRB high-z classification

# read in GRB functions
source('./algorithm1/algorithm1.R')

# read in data:
library(foreign)

##
## get the data in the desired form
##
filename = './algorithm1/uvot_no_error.arff'
data1 = read.arff(filename)
Z = data1$Z
data1$triggerid_str = NULL # meaningless feature in the data sets
data1 = removeErrors(data1)
data1 = cleanData(data1,4) # define above 4 as high, below 4 as low

### run rpart classifier
classes = data1[,1]
features = data1[,-1]

confmats = list()
# prior.high = seq(0.1,0.9,0.05)
# for(ii in 1:length(prior.high)){
#   tree.out = rpart.cv(features,classes,prior=c(1-prior.high[ii],prior.high[ii]),nfolds=10,seed=1)
#   confmats[[ii]]=tree.out$confmat
# }

# Documentation on rpart commands - mayo.edu/hsr/techrpt/61.pdf
# run RF classifier
# you will need to add the weights argument here
carttest = rpart.cv(features,classes,prior=c(0.45,0.55),nfolds=10,seed=1)
forest_order = NULL # save the probabilities-order output from random forests
forest_res = NULL # save the raw-probabilities output from random forests
weights_try = seq(1,101,10)
print(weights_try)
for(whigh in weights_try) {
	weights_vec = 1*(data1$class == "low") + whigh*(data1$class == "high")
	weights_vec = length(weights_vec) * weights_vec / sum(weights_vec)
	foresttest = rfc.cv(features,classes,nfolds=10,weights=weights_vec,seed=1)
	forest_res = cbind(forest_res,foresttest$predprob[,1])
	forest_order = cbind(forest_order,order(foresttest$predprob[,1]))
}

# Do bumps plot
color_vec = array(1,dim=length(data1$class))
lwd_vec = array(1,dim=length(data1$class))
for(i in seq(1,length(data1$class))) {
	if(data1$class[i] == "high") {
			color_vec[i] = 1
			lwd_vec[i] = 3
	} else {
			color_vec[i] = ceiling(1+Z[i])
			# Z = 0-1: red
			# Z = 1-2: green
			# Z = 2-3: blue
			# Z = 3-4: cyan
			# Z = 4-8: Black
	}
}
# plot here---as well as in eps file below
parcoord(forest_res,color_vec,lwd=lwd_vec,var.label=TRUE)

# write output
write(forest_res,"forest_probs_pred.txt") # write forest_res vector to text file
postscript(file="forest_probs_pred_bumps.eps",width=10,height=10) # save bumps plot
parcoord(forest_res,color_vec,lwd=lwd_vec,var.label=TRUE)
dev.off()

