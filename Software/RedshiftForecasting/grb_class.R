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

##########################################
####### GRB high-z classification #######
##########################################

####### read in GRB functions#######
source('./algorithm1/algorithm1.R')

####### read in data: #######
library(foreign)
library(fields)

read_data = function(filename='./algorithm1/uvot_no_error.arff',high_cutoff=4){
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
test_random_forest_weights = function(data_obj=NULL,log_weights_try=seq(0,5,0.5),seed=1){
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
   		foresttest = rfc.cv(data_obj$features,data_obj$classes,nfolds=10,weights=weights_vec,seed=seed)
   		forest_res = cbind(forest_res,foresttest$predprob[,1])
	}
   
	return(forest_res)
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
}

####### extract composite stats from the random forest seeds output ######
extract_stats = function(data_obj=NULL, log_weights_try=seq(0,5,0.5), forest_res_dir="redshift-output"){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
	# collect files in directory
	file_list = dir(forest_res_dir)
	
	# read in files
	Nweights = length(log_weights_try)
	forest_res = matrix(0,length(data_obj$Z),Nweights) # average across seeds
	nseed = 1
	for(file in file_list){
		forest_res_loc = as.matrix(read.table(paste(forest_res_dir,"/",file,sep=""), header=FALSE))
		forest_res = (1./nseed) * forest_res_loc + ((nseed-1.)/nseed) * forest_res
		nseed = nseed + 1
	}
	
	colnames(forest_res)=paste(log_weights_try)
	return(forest_res)
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
      if(alpha <= objective_cutoff){
         dumb_objective = base_objective_uvot_no*alpha/objective_cutoff
      }
      if(alpha > objective_cutoff){
         dumb_objective= base_objective_uvot_no+((alpha-objective_cutoff)*base_objective_uvot_yes)/(1-objective_cutoff)
      }
      print(dumb_objective)
   	## End dumb Classifier
		
		lines(log_weights_try,frac_found_loc,lty=1,lwd=2,col=col[nalpha])
		yline(alpha,lty=1,col=col[nalpha])
		yline(dumb_objective,lty=3,lwd=2,col=col[nalpha])
	}
	dev.off()
}

####### makes bumps plot, writes it to an image file, and saves the data that made it in a text file ####### 
make_bumps_plot = function(forest_res,data_obj=NULL,n_colors=128,z_width=3,imagefile="forest_probs_pred_bumps.pdf",textfile="forest_probs_pred.txt"){
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

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(10,1), heights=c(2,2)) # make a separate plot for colorbar
par(mar=c(4,2,0,0))
   parcoord(forest_res,lwd=lwd_vec,var.label=TRUE,col=tc[col.vec])
   title(xlab="log high-z weight",cex.lab=1.25,mgp=c(2.5,1,0)) # axis labels
   title(ylab="Pr(low-z GRB)",cex.lab=1.25,mgp=c(.25,1,0)) # yaxis
par(mar=c(4,0,0,1))
   plot(1, type="n", axes=F, xlab="z (log\n scale)", ylab="",xlim=c(-1,1),ylim=c(-1,1),mgp=c(1,1,0),cex.lab=1.25) # empty plot for colorbar
   colorbar.plot(0,0,strip=seq(min(logz),max(logz),length.out=n_colors),col=tc,horizontal=FALSE,strip.width=.6,strip.length=7.25) # plot colorbar
   text(0,-1,signif(10^min(logz)-1,2)) # add min and max to colorbar
   text(0,1,signif(10^max(logz)-1,2))
   abline(h= log10(hc+1)/(max(logz)-min(logz))/1.9,lwd=4) # plot z=4 cutoff (the 1.9 is a hack)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z > ",hc),pos=3)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z < ",hc),pos=1)
   
   pdf(file=imagefile,width=12,height=8) # save bumps plot
   
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(10,1), heights=c(2,2)) # make a separate plot for colorbar
par(mar=c(4,2,0,0))
   parcoord(forest_res,lwd=lwd_vec,var.label=TRUE,col=tc[col.vec])
   title(xlab="log high-z weight",cex.lab=1.25,mgp=c(2.5,1,0)) # axis labels
   title(ylab="Pr(low-z GRB)",cex.lab=1.25,mgp=c(.25,1,0)) # yaxis
par(mar=c(4,0,0,1))
   plot(1, type="n", axes=F, xlab="z (log\n scale)", ylab="",xlim=c(-1,1),ylim=c(-1,1),mgp=c(1,1,0),cex.lab=1.25) # empty plot for colorbar
   colorbar.plot(0,0,strip=seq(min(logz),max(logz),length.out=n_colors),col=tc,horizontal=FALSE,strip.width=.6,strip.length=7.25) # plot colorbar
   text(0,-1,signif(10^min(logz)-1,2)) # add min and max to colorbar
   text(0,1,signif(10^max(logz)-1,2))
   abline(h= log10(hc+1)/(max(logz)-min(logz))/1.9,lwd=4) # plot z=4 cutoff (the 1.9 is a hack)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z > ",hc),pos=3)
   text(0,log10(hc+1)/(max(logz)-min(logz))/1.9,paste("z < ",hc),pos=1)

   dev.off()
   
   # write forest_res vector to text file
   write(forest_res,"forest_probs_pred.txt")
}


# Calculate the number of GRBs we are allowing ourselves to follow-up
forest_run = function(data_obj=NULL,nfolds=10,alpha=0.3,mtry=NULL,weight=61,seed=1,n.trees=500){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
   if(is.null(mtry)){
     mtry = ceiling(sqrt(dim(data_obj$features)[2]))
   }
   print(mtry)
   num_of_grbs = length(data_obj$Z)

   weights_vec = 1*(data_obj$data1$class == "low") + weight*(data_obj$data1$class == "high")
   # the following might screw things up
   #weights_vec = length(weights_vec) * weights_vec / sum(weights_vec)
   # ff
   foresttest = rfc.cv(data_obj$features,data_obj$classes,nfolds=nfolds,weights=weights_vec,seed=seed,mtry=mtry)

   probs = foresttest$predprob[,1]
   num_to_follow = ceiling(alpha*num_of_grbs)
   # Grab the probability above which a GRB is considered low for following up
   # a given percentage of bursts
   mid_prob=sort(foresttest$predprob[,1])[num_to_follow]
   # Grab the array of Zs which are less than mid_prob 
   high_array=sort(data_obj$Z[foresttest$predprob[,1] <= mid_prob])
   # Calculate our objective function and confusion matrices
   num_actually_high=length(high_array[high_array >= data_obj$high_cutoff])
   pct_actually_high=num_actually_high/num_of_grbs
   objective = num_actually_high/num_to_follow

   high_as_high = num_actually_high
   high_as_low = data_obj$num_high - num_actually_high
   low_as_high = num_to_follow - num_actually_high
   low_as_low = data_obj$num_low - low_as_high
   print(high_as_high)
   print(low_as_high)
   print(low_as_low)
   print(high_as_low)
   return(list(objective=objective,high_as_high=high_as_high,high_as_low=high_as_low,low_as_low=low_as_low,low_as_high=low_as_high))
}




