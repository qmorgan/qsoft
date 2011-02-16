########################
## random forest code ##
########################
# Authors: T. Broderick, J. Long, A. Morgan, J. Richards


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

# Note: to remove global variables from namespace, use rm(list = ls(all = TRUE))

library(MASS)

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

## OUTDATED -  Problems with our dataset.  Use rforest.cv
# rfc.cv = function(x,y,nfolds=5,folds=NULL,testset=NULL,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
#   # don't train on any of the data in testset
#   # this is to use in the hierarchical classifier
#   require(party)
#   set.seed(seed)
#   
#   n = length(y)
#   p = length(table(y))
#   if(is.null(folds)){
#     folds = sample(1:nfolds,n,replace=TRUE)
#   }
#   else{
#     nfolds = length(unique(folds))
#   }
#   predictions = matrix(0,nrow=n,ncol=p)
# 
#   if(is.null(mtry)){
#     mtry = ceiling(sqrt(dim(x)[2]))
#   }
# 
#   for(ii in 1:nfolds){
#     print(paste("fold",ii,"of",nfolds))
#     leaveout = which(folds==ii)
#     train = cbind(y[-union(leaveout,testset)],x[-union(leaveout,testset),])
#     test = cbind(y[leaveout],x[leaveout,])
#     names(train)[1] = names(test)[1] = "y"
#     rf.tmp = cforest(y~.,data=train,weights=weights[-leaveout],controls=cforest_classical(mtry=mtry,ntree=n.trees))
#     predictions[leaveout,] = matrix(unlist(treeresponse(rf.tmp,newdata=test)),length(leaveout),p,byrow=T)
#   }
#   pred = levels(y)[apply(predictions,1,which.max)]
#   pred = factor(pred,levels=levels(y))
#   confmat = table(pred,y)
#   err.rate = 1-sum(diag(confmat))/n
# 
#   return(list(predclass=pred,predprob=predictions,confmat=confmat,err.rate=err.rate))
# }

##########################################
####### GRB high-z classification #######
##########################################

####### read in GRB functions#######
source('./algorithm1/algorithm1.R')

####### read in data: #######
library(foreign)
library(fields)

read_data = function(filename='./Data/GRB_short+outliers+noZ_removed_reduced.arff',high_cutoff=4){
   data1 = read.arff(filename)
   Z = data1$Z
   ####### define above high_cutoff as high, below as low $ ####### 
   num_high = length(Z[Z >= high_cutoff])
   num_low = length(Z[Z < high_cutoff])
   data1$triggerid_str = NULL # meaningless feature in the data sets
   data1 = removeErrors(data1)
   data1 = cleanData(data1,high_cutoff) 
   classes = data1[,1]
   features = data1[,-1]
   confmats = list()
   print(paste("Read in:",filename))
   print(paste("with a high-z cutoff value of",high_cutoff))
   return(list(data1=data1,num_high=num_high,num_low=num_low,classes=classes,features=features,confmats=confmats,Z=Z,high_cutoff=high_cutoff))
}


####### run rpart CART classifier ####### 
test_cart = function(data_obj=NULL,seed=1){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not defined; using default values")
      data_obj = read_data()
   }
   ###########################################################################
   # prior.high = seq(0.1,0.9,0.05)
   # for(ii in 1:length(prior.high)){
   #   tree.out = rpart.cv(features,classes,prior=c(1-prior.high[ii],prior.high[ii]),nfolds=10,seed=1)
   #   confmats[[ii]]=tree.out$confmat
   # }
   carttest = rpart.cv(data_obj$features,data_obj$classes,prior=c(0.45,0.55),nfolds=10,seed=1)
   # Documentation on rpart commands - mayo.edu/hsr/techrpt/61.pdf
}

####### run Random Forest classifier over a vector of weights ####### 
test_random_forest_weights = function(data_obj=NULL,log_weights_try=seq(0,5,0.5),seed=1,stratified=FALSE){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not defined; using default values")
      data_obj = read_data()
   }
   ###########################################################################
	forest_res = NULL # save the raw-probabilities output from random forests

	weights_try = 10^log_weights_try
	Nweights = length(weights_try)
	for(nweight in seq(1,Nweights)) {
   		print(paste("*** weight",nweight,"of",Nweights,"***"))
   		whigh = weights_try[nweight]
   		weights_vec = 1*(data_obj$data1$class == "low") + whigh*(data_obj$data1$class == "high")
   		#	weights_vec = length(weights_vec) * weights_vec / sum(weights_vec) # REMOVE? Causing problems
                # stratified folds (for high-z bursts)
                # The number of folds is determined by the number of high z bursts
                # As of 2-9-11, this should no longer be important; could go down to regular CV?
         if(stratified==TRUE){
            n.high = sum(data_obj$classes=="high")
            folds = sample(1:n.high,length(data_obj$classes),replace=TRUE)
            folds[data_obj$classes=="high"]=1:n.high
            folds[data_obj$classes=="low"]=rep(1:n.high,length.out=sum(data_obj$classes=="low"))       
   		   foresttest = forest.cv(data_obj$features,data_obj$classes,folds=folds,weights=weights_vec,seed=seed)
   		}
   		else{
   		   foresttest = forest.cv(data_obj$features,data_obj$classes,nfolds=10,weights=weights_vec,seed=seed)
      	}
   		forest_res = cbind(forest_res,foresttest)
	}
   
	return(forest_res)
}

noisify_residuals = function(forest_res){
   numZ = dim(forest_res)[1]
   numTrials = dim(forest_res)[2]
   forest_res_rand = forest_res + matrix(runif(numTrials*numZ,1E-6,1E-5),numZ,numTrials)
   return(forest_res_rand)
}

order_residuals = function(forest_res,reverse=FALSE){
   forest_res_ordered=apply(forest_res,2,rank)
   if(reverse==TRUE){
      numinstances = length(forest_res[,1])
      forest_res_ordered = forest_res_ordered*-1 + numinstances + 1
   }
   return(forest_res_ordered)
}

####### smooth weighted random forest classifiers over a number of seeds ####### 
smooth_random_forest_weights = function(data_obj=NULL,log_weights_try=seq(0,5,0.5),Nseeds=10,results_dir="redshift-output"){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
	for(nseed in seq(1,Nseeds)) {
   		print(paste("@@@@@@ seed",nseed,"of",Nseeds,"@@@@@@"))
		forest_res_loc = test_random_forest_weights(data_obj=data_obj,log_weights_try=log_weights_try,seed=nseed) # result for this seed
		
		# save results to file
		results_file = paste(results_dir,"/",nseed,".txt",sep="")
		write(t(forest_res_loc), results_file, append=FALSE, ncolumns=length(log_weights_try))
	}
   return(forest_res_loc)
 }


extract_stats = function(data_obj=NULL, log_weights_try=seq(0,5,0.5), forest_res_dir="smooth_weights_reduced"){
   ##### If data object is not defined, create the default data object ######
   ##### Results are then stored in the fres list within the data object ####
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
	# collect files in directory
	file_list = dir(forest_res_dir)
	
	# read in files
	Nweights = length(log_weights_try)
	# initialize matrix of zeros to average across seeds
	res_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	ord_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	objective_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	
	####### Grab info
	high_cutoff=data_obj$high_cutoff
	Nweights = length(log_weights_try)
	Nz = length(data_obj$Z)
	Nhigh = sum(data_obj$Z > high_cutoff)
	alpha_vec=seq(0,Nz-1)/(Nz-1)
	Nalpha = length(alpha_vec)
	col = rainbow(Nalpha)
	random_guess_vec = (data_obj$num_high)/Nz*alpha_vec

	#######
	num_of_seeds = length(file_list)
	big_forest_res = array(0,dim=c(Nz,Nweights,num_of_seeds))
	# get ready to loop over all seeds
	nseed = 1
	for(file in file_list){
		forest_res_loc = as.matrix(read.table(paste(forest_res_dir,"/",file,sep=""), header=FALSE))
		# forest_res_loc_ordered = order_residuals(forest_res_loc)
		# forest_res = (1./nseed) * forest_res_loc + ((nseed-1.)/nseed) * forest_res
		big_forest_res[,,nseed] = forest_res_loc
		nseed = nseed + 1
	}
	
	# the following should be the same as the old forest res function
	for(weight_index in seq(1,Nweights)){
	   res_avg_over_seeds[,weight_index] = rowSums(big_forest_res[,weight_index,])/num_of_seeds
	   res_ordered = order_residuals(noisify_residuals(big_forest_res[,weight_index,]))
	   ord_avg_over_seeds[,weight_index] = rowSums(res_ordered)/num_of_seeds
	   
	   # populate the objective_avg_over_seeds array
	   # each row in this array is an ALPHA*Nz value (not a grb instance like the others)
	   # if we observe alpha% of events, how many high-as-high/actual_high can we expect?
	   for(nalpha in seq(1,Nalpha)) {
   		alpha = alpha_vec[nalpha]
   		rand_guess = random_guess_vec[nalpha]
   		# find number of follow ups for this alpha
   		Nfollow = floor(Nz*alpha)
   		# find high-z (>4) that are in the alpha threshold
   		Nfound_loc = colSums((as.matrix(data_obj$Z)%*%matrix(1,1,num_of_seeds) > high_cutoff)&(res_ordered<=Nfollow))
   		avg_frac_found= sum(Nfound_loc / (1.*Nhigh))/num_of_seeds
   		objective_avg_over_seeds[nalpha,weight_index] = avg_frac_found
		}
		
	   
	}
	
	colnames(res_avg_over_seeds)=paste(log_weights_try)
	
	forest_res_obj = list(full=big_forest_res,avg_over_seeds=res_avg_over_seeds, ord_avg_over_seeds = ord_avg_over_seeds, objective_avg_over_seeds = objective_avg_over_seeds)
   data_obj$fres = forest_res_obj
	
	return(data_obj)
}

####### makes objective function plot ####### 
# forest_order is ordered using "order_residuals" function with high-z at low rank numbers
make_obj_fcn_plot = function(forest_order,data_obj=NULL,alpha_vec=seq(0.1,0.9,0.1),log_weights_try=seq(0,5,0.5),imagefile="objective_fcn.pdf"){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
      high_cutoff=data_obj$high_cutoff
   	Nweights = length(log_weights_try)
   	Nz = length(data_obj$Z)
   	Nhigh = sum(data_obj$Z > high_cutoff)
	pdf(file=imagefile,width=12,height=8) # save obj plot
	plot(x = c(min(log_weights_try), max(log_weights_try)), y = c(0,1), xlim = c(min(log_weights_try), max(log_weights_try)), ylim=c(0,1), pch="") # initialize plot))
	Nalpha = length(alpha_vec)
	col = rainbow(Nalpha)

	random_guess_vec = (data_obj$num_high)/Nz*alpha_vec
	for(nalpha in seq(1,Nalpha)) {
		alpha = alpha_vec[nalpha]
		rand_guess = random_guess_vec[nalpha]
		# find number of follow ups for this alpha
		Nfollow = floor(Nz*alpha)
		# find high-z (>4) that are in the alpha threshold
		Nfound_loc = colSums((as.matrix(data_obj$Z)%*%matrix(1,1,Nweights) > high_cutoff)&(forest_order<=Nfollow))
		frac_found_loc = Nfound_loc / (1.*Nhigh)
		
		## Dumb Classifier - based soley on UVOT detection ##
      uvd=data_obj$data1$uvot_detection
      N_uvot_yes = length(uvd[uvd=='yes'])
      N_uvot_no = length(uvd[uvd=='no'])
      base_objective_uvot_no = 13./17.
      base_objective_uvot_yes = 4./17.
      objective_cutoff = N_uvot_no/(N_uvot_no+N_uvot_yes) 
      if(N_uvot_yes == 0 && N_uvot_no ==0){
         print('UVOT DATA NOT INCLUDED IN SAMPLE')
         print('Using 49/151 = 0.3245033 as objective cutoff for dumb classifier')
         objective_cutoff = 0.3245033
      }
      if(alpha <= objective_cutoff){
         dumb_objective = base_objective_uvot_no*alpha/objective_cutoff
      }
      if(alpha > objective_cutoff){
         dumb_objective= base_objective_uvot_no+((alpha-objective_cutoff)*base_objective_uvot_yes)/(1-objective_cutoff)
      }
   	## End dumb Classifier
		
		lines(log_weights_try,frac_found_loc,lty=1,lwd=4,col=col[nalpha])
		yline(alpha,lty=1,col=col[nalpha])
		yline(dumb_objective,lty=3,lwd=3,col=col[nalpha])
		
		
	}
	legend(4.8,0.19,c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"), cex=0.6,col=col,lty=1:1)
	
	dev.off()
}

####### makes bumps plot, writes it to an image file, and saves the data that made it in a text file ####### 
make_bumps_plot = function(res_avg_over_seeds,data_obj=NULL,xlabel="log high-z weight",ylabel=expression(widehat(alpha)),n_colors=128,z_width=3,imagefile="forest_probs_pred_bumps.pdf",textfile="forest_probs_pred.txt"){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   hc = data_obj$high_cutoff
   ###########################################################################
   ####### PLOTTING - here and in EPS ####### 
   #######  Set up color for bumps plot ####### 
   color_vec = array(1,dim=length(data_obj$data1$class))
   lwd_vec = array(1,dim=length(data_obj$data1$class))
   for(i in seq(1,length(data_obj$data1$class))) {
   	if(data_obj$data1$class[i] == "high") {
   			color_vec[i] = 1
   			lwd_vec[i] = z_width # High-z bursts are thicker
   	} else {
   			color_vec[i] = ceiling(1+data_obj$Z[i])
   	}
   }

   logz=log10(1+data_obj$Z)


   ####### The next few lines are for coloring ####### 
   ####### Break up log(1+z) into 64 color bins ####### 
   col.vec=ceiling((logz-min(logz))/(max(logz)-min(logz))*n_colors)
   for(n in 1:length(col.vec)) {
      if(col.vec[n] > 1) {
         col.vec[n] = col.vec[n]
      } else {
         col.vec[n] = 1
      }
   }
   tc = tim.colors(n_colors)

   pdf(file=imagefile,width=12,height=8) # save bumps plot
   
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(10,1), heights=c(2,2)) # make a separate plot for colorbar
par(mar=c(4,2,0,0))
   parcoord(res_avg_over_seeds,lwd=lwd_vec,var.label=TRUE,col=tc[col.vec])
   title(xlab=xlabel,cex.lab=1.25,mgp=c(2.5,1,0)) # axis labels
   title(ylab=ylabel,cex.lab=1.5,mgp=c(.25,1,0)) # yaxis
par(mar=c(4,0,0,1))
   plot(1, type="n", axes=F, xlab="z (log\n scale)", ylab="",xlim=c(-1,1),ylim=c(-1,1),mgp=c(1,1,0),cex.lab=1.25) # empty plot for colorbar
   colorbar.plot(0,0,strip=seq(min(logz),max(logz),length.out=n_colors),col=tc,horizontal=FALSE,strip.width=.6,strip.length=7.25) # plot colorbar
   text(0,-1,signif(10^min(logz)-1,2)) # add min and max to colorbar
   text(0,1,signif(10^max(logz)-1,2))
   abline(h= log10(hc+1)/(max(logz)-min(logz))/1.9,lwd=4) # plot z=4 cutoff (the 1.9 is a hack)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z > ",hc),pos=3)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z < ",hc),pos=1)

   dev.off()
   
   # write res_avg_over_seeds vector to text file
   write(res_avg_over_seeds,textfile)
}

## OUTDATED?
# # Calculate the number of GRBs we are allowing ourselves to follow-up
# forest_run = function(data_obj=NULL,nfolds=10,alpha=0.3,mtry=NULL,weight=61,seed=1,n.trees=500){
#    ##### If data object is not defined, create the default data object ######
#    if(is.null(data_obj)){
#       print("data_obj not specified; using default values")
#       data_obj = read_data()
#    }
#    ###########################################################################
#    if(is.null(mtry)){
#      mtry = ceiling(sqrt(dim(data_obj$features)[2]))
#    }
#    print(mtry)
#    num_of_grbs = length(data_obj$Z)
# 
#    weights_vec = 1*(data_obj$data1$class == "low") + weight*(data_obj$data1$class == "high")
#    # the following might screw things up
#    #weights_vec = length(weights_vec) * weights_vec / sum(weights_vec)
#    # ff
# #   foresttest = rfc.cv(data_obj$features,data_obj$classes,nfolds=nfolds,weights=weights_vec,seed=seed,mtry=mtry)
#                 # stratified folds (for high-z bursts)
#    n.high = sum(data_obj$classes=="high")
#    folds = sample(1:n.high,length(data_obj$classes),replace=TRUE)
#    folds[data_obj$classes=="high"]=1:n.high
#    foresttest = rfc.cv(data_obj$features,data_obj$classes,folds=folds,weights=weights_vec,seed=seed,mtry=mtry)
# 
#    probs = foresttest$predprob[,1]
#    num_to_follow = ceiling(alpha*num_of_grbs)
#    # Grab the probability above which a GRB is considered low for following up
#    # a given percentage of bursts
#    mid_prob=sort(foresttest$predprob[,1])[num_to_follow]
#    # Grab the array of Zs which are less than mid_prob 
#    high_array=sort(data_obj$Z[foresttest$predprob[,1] <= mid_prob])
#    # Calculate our objective function and confusion matrices
#    num_actually_high=length(high_array[high_array >= data_obj$high_cutoff])
#    pct_actually_high=num_actually_high/num_of_grbs
#    objective = num_actually_high/num_to_follow
# 
#    high_as_high = num_actually_high
#    high_as_low = data_obj$num_high - num_actually_high
#    low_as_high = num_to_follow - num_actually_high
#    low_as_low = data_obj$num_low - low_as_high
#    print(high_as_high)
#    print(low_as_high)
#    print(low_as_low)
#    print(high_as_low)
#    return(list(objective=objective,high_as_high=high_as_high,high_as_low=high_as_low,low_as_low=low_as_low,low_as_high=low_as_high))
# }

# fit and return forest on training data (using optimal tuning parameters)
forest.fit = function(x,y,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
  require(party)
  set.seed(seed)
  
  n = length(y)
  p = length(table(y))
  if(is.null(mtry)){ # default for mtry
    mtry = ceiling(sqrt(dim(x)[2]))
  }
  train = cbind(y,x) # set up data to read into cforest
  names(train)[1] = "y"
  # fit random forest
  rf.fit = cforest(y~.,data=train,weights=weights,controls=cforest_classical(mtry=mtry,ntree=n.trees))

  return(rf.fit)
}

# Predict redshift (alpha-hat) for new GRBs
forest.pred = function(forest,xnew){
  xnew = as.data.frame(xnew)
  n.new = dim(xnew)[1]
  n.old = length(predict(forest))
  # predictions for training data (to compute alpha-hats)
  pred.train = matrix(unlist(treeresponse(forest)),n.old,2,byrow=T) # CV this?
  # predict post. probs. for new data, with input forest
  predictions = matrix(unlist(treeresponse(forest,newdata=xnew)),n.new,2,byrow=T)
  alpha.hat = NULL # compute alpha-hat values
  for(ii in 1:n.new){
    alpha.hat = c(alpha.hat, sum(predictions[ii,2]< pred.train[,2])/n.old)
  }

  return(list(alpha.hat = alpha.hat,prob.high=predictions[,2],prob.low=predictions[,1]))
}

forest.cv = function(x,y,nfolds=10,folds=NULL,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
  require(party)
  set.seed(seed)
  
  n = length(y)
  p = length(table(y))
  if(is.null(folds)){
    folds = sample(1:nfolds,n,replace=TRUE)
  }
  else{
    nfolds = length(unique(folds))
  }
  predictions = rep(0,p)

  if(is.null(mtry)){
    mtry = ceiling(sqrt(dim(x)[2]))
  }

  for(ii in 1:nfolds){
    print(paste("fold",ii,"of",nfolds))
    leaveout = which(folds==ii)
    rf.tmp = forest.fit(x[-leaveout,],y[-leaveout],mtry=mtry,weights=weights[-leaveout],n.trees=n.trees,seed=seed)
    predictions[leaveout] = forest.pred(rf.tmp,x[leaveout,])$alpha.hat
  }
  return(predictions)
}


# roc curves
# 1. true class is n length vector of high / low
# 2. prediction_matrix is n x c matrix of probabilities (of low, this needs
#    to be thought about, otherwise curve goes down)
make_roc_curve = function(true_class,prediction_matrix,curve_colors=NULL,filename="roc_curve.pdf"){

  # check to make sure data is in correct form
  require(ROCR)
  if(!is.factor(true_class)){
    print("true_class must be a factor")
    return(0)
  }
  if(!is.matrix(prediction_matrix)){
    print("prediction_matrix must be a matrix")
    return(0)
  }
  if(nrow(prediction_matrix) != length(true_class)){
    print("the number of obs you are predicting != number of true classes")
    return(0)
  }
  
  # if colors not specified, get some
  if(is.null(curve_colors)){
    curve_colors = 1:ncol(prediction_matrix)
  }

  # use functions in ROCR package to make objects for plotting
  true_class = matrix(rep(true_class,ncol(prediction_matrix)),nrow=nrow(prediction_matrix),byrow=F)
  pred <- prediction(prediction_matrix,true_class,label.ordering=c("low","high"))
  # see ROCR user guide on CRAN p.2 for definition of "tpr", "pcfall", ect.
  performance1 = performance(pred,"tpr","pcfall")

  # make the plot
  pdf(filename)
  plot(performance1,col=as.list(curve_colors),xlab="False Discovery Rate = Contamination = False High / Total Predicted High",ylab="True Positive Rate = Efficiency = True High / Total Actual High)",main="ROC Curve for Classifiers")
  dev.off()
}

efficiency_vs_alpha = function(data_obj,weight_index=5,imagefile='test'){
   # Take the fifth weight for now 
   avg_obj = data_obj$fres$objective_avg_over_seeds[,weight_index]
   Zlen_1 = length(avg_obj) - 1
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   alpha_tries = seq(0,1,1/Zlen_1)
   pdf(imagefile)
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up: "~alpha), ylab=expression('Fraction of high (z>4) GRBs observed'), pch="") # initialize plot))
   title(main=expression("Efficiency vs"~alpha), sub=data_obj$data_string)
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   lines(alpha_tries, avg_obj, lty=1, lwd=4, col='red')
   dev.off()
   
}

multiple_efficiency_vs_alpha = function(data_obj_list,weight_index=5,imagefile='./Plots/ROC_multi.pdf'){
   pdf(imagefile)
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up: "~alpha), ylab=expression('Fraction of high (z>4) GRBs observed'), pch="") # initialize plot))
   title(main=expression("Efficiency vs"~alpha))
   n_curves = length(data_obj_list)
   col = rainbow(n_curves)
   name_list = c()
   curve_index = 1
   for(data_obj_name in ls(data_obj_list)){
      print(data_obj_name)
      data_obj = get(data_obj_name,pos=data_obj_list)
      avg_obj = data_obj$fres$objective_avg_over_seeds[,weight_index]
      Zlen_1 = length(avg_obj) - 1
      print(Zlen_1)
      alpha_tries = seq(0,1,1/Zlen_1)
      lines(alpha_tries, avg_obj, lty=1, lwd=4, col=col[curve_index])
      curve_index = curve_index + 1
      name_list = c(name_list,data_obj$data_string)
   }
   # now print diagonal line for reference
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   
   print(name_list)
   legend(0.5,0.39, name_list, cex=1.2,col=col,lty=1)
	
   dev.off()
}

# Wrapper to make all representative plots for a given dataset
make_forest_plots = function(data_string="reduced",generate_data=FALSE){
   # generate_data will re-do the smooth_random_forest_weights function, which takes a while
   data_filename = paste("./Data/GRB_short+outliers+noZ_removed_",data_string,".arff",sep="")
   data_results_dir = paste("smooth_weights_",data_string,sep="")
   obj_func_name = paste("./Plots/objective_fcn_",data_string,".pdf",sep="")
   bumps_pred_plot_name = paste("./Plots/forest_pred_bumps_",data_string,".pdf",sep="")
   bumps_rand_plot_name = paste("./Plots/forest_rand_bumps_",data_string,".pdf",sep="")
   bumps_plot_name = paste("./Plots/forest_order_bumps_",data_string,".pdf",sep="")
   bumps_pred_text_name = paste("forest_pred_bumps_",data_string,".txt",sep="")
   bumps_text_name = paste("forest_order_bumps_",data_string,".txt",sep="")
   roc_plot_name = paste("./Plots/ROC_",data_string,".pdf",sep="")
   mydata = read_data(filename=data_filename,high_cutoff=4)
   mydata$data_string = data_string
   if(generate_data == TRUE){
      alphas.cv = smooth_random_forest_weights(data_obj = mydata,results_dir=data_results_dir)
   }
   mydata = extract_stats(data_obj = mydata, forest_res_dir=data_results_dir)
   efficiency_vs_alpha(mydata,imagefile=roc_plot_name)
   fres = mydata$fres$avg_over_seeds
   make_bumps_plot(fres, data_obj = mydata, imagefile=bumps_pred_plot_name,textfile=bumps_pred_text_name)
   fres_rand = noisify_residuals(fres)
   make_bumps_plot(fres_rand, data_obj = mydata, imagefile=bumps_rand_plot_name)
   fres_ordered = order_residuals(fres_rand)
   make_obj_fcn_plot(fres_ordered, data_obj = mydata, imagefile=obj_func_name)
   fres_ordered = order_residuals(fres_rand,reverse=TRUE)
   make_bumps_plot(fres_ordered, data_obj = mydata, ylabel=paste(expression(widehat(alpha)),' (Normalized)'),imagefile=bumps_plot_name,textfile=bumps_text_name)
   ## make_bumps_plot(alphas.cv)
 }
 
make_efficiency_plots = function(generate_data=FALSE, data_string_list=list('reduced','UVOTonly','UVOTandZpred','Nat_Zprediction','Full','reduced_nozpredict')){
   roc_plot_name = paste("./Plots/ROC_Multi.pdf",sep="")
   
   curve_index = 1
   data_obj_list = list()
   
   for(data_string in data_string_list){
      data_filename = paste("./Data/GRB_short+outliers+noZ_removed_",data_string,".arff",sep="")
      data_results_dir = paste("smooth_weights_",data_string,sep="")
      
      mydata = read_data(filename=data_filename,high_cutoff=4)
      mydata$data_string = data_string
      print(data_string)
      if(generate_data == TRUE){
         alphas.cv = smooth_random_forest_weights(data_obj = mydata,results_dir=data_results_dir)
      }
      objname = paste('o',curve_index,sep="")
      mydata = extract_stats(data_obj = mydata, forest_res_dir=data_results_dir)
     # data_obj_list$data_string = mydata
      # assign(data_string, mydata)
      # append(data_obj_list,mydata)
      # WOW this is poor programming.  but it's late and I couldn't get it to work otherwise.
      if(curve_index==1){data_obj_list$l1 = mydata}
      if(curve_index==2){data_obj_list$l2 = mydata}
      if(curve_index==3){data_obj_list$l3 = mydata}
      if(curve_index==4){data_obj_list$l4 = mydata}
      if(curve_index==5){data_obj_list$l5 = mydata}
      if(curve_index==6){data_obj_list$l6 = mydata}
      if(curve_index==7){data_obj_list$l7 = mydata}
      if(curve_index==8){data_obj_list$l8 = mydata}
      curve_index = curve_index + 1
      print(curve_index)
   }
   multiple_efficiency_vs_alpha(data_obj_list)
}