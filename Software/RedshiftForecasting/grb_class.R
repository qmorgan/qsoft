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
library(randomForest)

read_data = function(filename='./Data/GRB_short+outliers+noZ_removed_reduced.arff',high_cutoff=4,impute=TRUE){
   data1 = read.arff(filename)
   if(impute == TRUE){
     data1 = na.roughfix(data1)
   }
   Z = data1$Z
   ####### define above high_cutoff as high, below as low $ ####### 
   num_high = length(Z[Z >= high_cutoff])
   num_low = length(Z[Z < high_cutoff])
   triggerids = data1$triggerid_str
   data1$triggerid_str = NULL # meaningless feature in the data sets
   data1 = removeErrors(data1)
   data1 = cleanData(data1,high_cutoff) 
   classes = data1[,1]
   features = data1[,-1]
   confmats = list()
   print(paste("Read in:",filename))
   print(paste("with a high-z cutoff value of",high_cutoff))
   return(list(data1=data1,filename=filename,num_high=num_high,num_low=num_low,classes=classes,features=features,confmats=confmats,Z=Z,triggerids=triggerids,high_cutoff=high_cutoff))
}

##### adds useless features to an arff file - very simple
add_useless_features = function(input_filename='./Data/GRB_short+outliers+noZ_removed_UVOTonly.arff', output_filename=NULL, number_useless_features=1){
  # need read.arff and write.arff
  require(foreign)

  # if the output file name is not specified, create name
  if(is.null(output_filename)){
    if( grepl(".arff",input_filename,fixed=TRUE) ){
      to_add = paste("_num_useless",number_useless_features,".arff",
                     sep="")
      output_filename = sub(".arff",to_add,input_filename,fixed=TRUE)
    }
    else {
      print('input file must be .arff')
      return()
    }
  }

  # get data, add useless features, write data to arff  
  input_data = read.arff(input_filename)
  useless_features = matrix(runif(nrow(input_data) * number_useless_features),nrow=nrow(input_data))
  colnames(useless_features) = paste("useless",1:number_useless_features,sep="")
  input_data = data.frame(input_data,useless_features)
  write.arff(input_data,output_filename)
  
}

remake_all_useless = function(input_filename='./Data/GRB_short+outliers+noZ_removed_UVOTonly.arff'){
   # 
   add_useless_features(number_useless_features=10)
   add_useless_features(number_useless_features=20)
   add_useless_features(number_useless_features=30)
   add_useless_features(number_useless_features=40)
   add_useless_features(number_useless_features=50)
   add_useless_features(number_useless_features=60)
   add_useless_features(number_useless_features=70)
   add_useless_features(number_useless_features=80)
   add_useless_features(number_useless_features=90)
   add_useless_features(number_useless_features=1)
   add_useless_features(number_useless_features=2)
   add_useless_features(number_useless_features=4)
   add_useless_features(number_useless_features=8)
   add_useless_features(number_useless_features=16)
   add_useless_features(number_useless_features=32)
   add_useless_features(number_useless_features=64)
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
test_random_forest_weights = function(data_obj=NULL,log_weights_try=seq(-1,1,0.2),seed=1,stratified=FALSE){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not defined; using default values")
      data_obj = read_data()
   }
   ###########################################################################
	forest_res = NULL # save the raw-probabilities output from random forests

   # weights_try = c(NA,10^log_weights_try)
   weights_try = 10^log_weights_try
   weights_try[1] = NA
	print(weights_try)
	Nweights = length(weights_try)
	for(nweight in seq(1,Nweights)) {
   		print(paste("*** weight",nweight,"of",Nweights,"***"))
   		whigh = weights_try[nweight]
   		if(!is.na(whigh)){
   		   weights_vec = 1*(data_obj$data1$class == "low") + whigh*(data_obj$data1$class == "high")
   		}
   		else{weights_vec = NULL}
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
# Add a very small random number to the residuals to avoid duplicates when ordering
   numZ = dim(forest_res)[1]
   numTrials = dim(forest_res)[2]
   forest_res_rand = forest_res + matrix(runif(numTrials*numZ,1E-6,1E-5),numZ,numTrials)
   return(forest_res_rand)
}

order_residuals = function(forest_res,reverse=FALSE){
# Return the index order of the residuals from low to high
   forest_res_ordered=apply(forest_res,2,rank)
   if(reverse==TRUE){
      numinstances = length(forest_res[,1])
      forest_res_ordered = forest_res_ordered*-1 + numinstances + 1
   }
   return(forest_res_ordered)
}

####### smooth weighted random forest classifiers over a number of seeds ####### 
smooth_random_forest_weights = function(data_obj=NULL,log_weights_try=seq(-1,1,0.2),first_seed=1,Nseeds=10,results_dir="redshift-output", redo_useless=FALSE){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   # Create the results directory if it doesn't exist
   if(file.exists(results_dir) == FALSE){
      dir.create(results_dir)
   }
   ###########################################################################
	for(nseed in seq(first_seed,Nseeds)) {
   		print(paste("@@@@@@ seed",nseed,"of",Nseeds,"@@@@@@"))
   	
   	# Remake the useless datasets for better randomness
   	if(redo_useless == TRUE){
   	   remake_all_useless()
   	   data_obj=read_data(filename=data_obj$filename)
   	}
   	
		forest_res_loc = test_random_forest_weights(data_obj=data_obj,log_weights_try=log_weights_try,seed=nseed) # result for this seed
		
		# save results to file
		results_file = paste(results_dir,"/",nseed,".txt",sep="")
		write(t(forest_res_loc), results_file, append=FALSE, ncolumns=length(log_weights_try))
	}
   return(forest_res_loc)
 }


extract_stats = function(data_obj=NULL, forest_res_dir="./smooth_weights_results/smooth_weights_reduced_4"){
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
	testloc = as.matrix(read.table(paste(forest_res_dir,"/",'1.txt',sep="")))
	Nweights = dim(testloc)[2]
	# initialize matrix of zeros to average across seeds
	res_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	res_sd_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	res_median_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	ord_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights) 
	objective_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	objective_stdev_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	purity_avg_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	purity_stdev_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	purity_upper_perc_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	purity_lower_perc_over_seeds = matrix(0,length(data_obj$Z),Nweights)
	
	####### Grab info
	high_cutoff=data_obj$high_cutoff

	Nz = length(data_obj$Z)
	Nhigh = sum(data_obj$Z > high_cutoff)
	alpha_vec=seq(0,Nz-1)/(Nz-1)
	#alpha_vec=seq(1,Nz)/(Nz)
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
	   # DOUBLE CHECK THAT THIS IS CORRECT WAY TO DO ROW MEDIAN AND SD - it is.
	   res_median_over_seeds[,weight_index] = apply(big_forest_res[,weight_index,],1,median)
	   res_sd_over_seeds[,weight_index] = apply(big_forest_res[,weight_index,],1,sd)
	   # res_ordered is an (N_bursts x N_seeds) array; each column is the ranking of each burst against the others from 1 to N_bursts
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
         # Old method: Changed because this code was opaque to me and I couldn't extract the stdev
         # Nfound_loc = colSums((as.matrix(data_obj$Z)%*%matrix(1,1,num_of_seeds) > high_cutoff)&(res_ordered<=Nfollow))
         # avg_frac_found= sum(Nfound_loc / (1.*Nhigh))/num_of_seeds
   		
   		
   		# printing for testing purposes
         # if(weight_index==5 & nalpha == 35){
         #    as_high_all_seeds = (res_ordered<=Nfollow)
         #    high_as_high_all_seeds = as_high_all_seeds[data_obj$Z >= high_cutoff,1:num_of_seeds]
         #    summed_as_high_all_seeds = rowSums(as_high_all_seeds)
         #    summed_high_as_high = summed_as_high_all_seeds[data_obj$Z >= high_cutoff]
         #    average_high_as_high = sum(summed_high_as_high/1700)
         #    mean_as_high_all_seeds = apply(as_high_all_seeds,1,mean)
         #    mean_high_as_high_all_seeds = apply(high_as_high_all_seeds,1,mean) # [1xNalpha array]
         #    mean_high_as_high_all_seeds_dim2 = apply(high_as_high_all_seeds,2,mean) # [1xNseeds array]
         #    print(mean_high_as_high_all_seeds_dim2)
         #             
         #             # the following 3 are the same
         #    print(mean_as_high_all_seeds[data_obj$Z >= high_cutoff])
         #    print(mean_high_as_high_all_seeds)
         #    print(summed_high_as_high/num_of_seeds)
         #    
         #    # the following 3 are the same; though we can get the std from the last
         #    print(average_high_as_high)
         #    print(avg_frac_found)
         #    print(mean(mean_high_as_high_all_seeds_dim2))
         #    print(sd(mean_high_as_high_all_seeds_dim2))
         #    
         # }
         
   		# Find bursts classified as high
   		as_high_all_seeds = (res_ordered<=Nfollow)
         # Grab the ones that were actually high
   		high_as_high_all_seeds = as_high_all_seeds[data_obj$Z >= high_cutoff,1:num_of_seeds]
         # average these 17 bursts together
   		mean_high_as_high_all_seeds_dim2 = apply(high_as_high_all_seeds,2,mean) # [1xNseeds array]
         # Take the mean of these 100 averages
		   avg_frac_found = mean(mean_high_as_high_all_seeds_dim2)
         # Take the stdev of these 100 averages
		   stdev_frac_found = sd(mean_high_as_high_all_seeds_dim2)
		   
   		objective_stdev_over_seeds[nalpha,weight_index] = stdev_frac_found
   		objective_avg_over_seeds[nalpha,weight_index] = avg_frac_found
   		
   		sum_high_as_high_all_seeds_dim2 = apply(high_as_high_all_seeds,2,sum) # [1xNseeds array]
   		avg_sum_found = mean(sum_high_as_high_all_seeds_dim2)
		   stdev_sum_found = sd(sum_high_as_high_all_seeds_dim2)
   		
   		sum_upper=quantile(sum_high_as_high_all_seeds_dim2, c(0.84))
   		sum_lower=quantile(sum_high_as_high_all_seeds_dim2, c(0.16))
   		
   		if(weight_index==5 & nalpha == 9){
   		   print(sum_high_as_high_all_seeds_dim2)
   		   print(sum_upper)
   		   print(sum_lower)
   		   print(avg_sum_found)
   		   print(stdev_sum_found)
   		   print(avg_sum_found/Nfollow)
		   }
   		
   		purity_upper_perc_over_seeds[nalpha,weight_index]=sum_upper/Nfollow
   		purity_lower_perc_over_seeds[nalpha,weight_index]=sum_lower/Nfollow
   		purity_avg_over_seeds[nalpha,weight_index] = avg_sum_found/Nfollow
   		purity_stdev_over_seeds[nalpha,weight_index] = stdev_sum_found/Nfollow
		   
		}
		
	   
	}
	# is this necessary?
	#colnames(res_avg_over_seeds)=paste(log_weights_try)
	
	forest_res_obj = {list(full=big_forest_res,avg_over_seeds=res_avg_over_seeds, 
	   sd_over_seeds=res_sd_over_seeds, median_over_seeds=res_median_over_seeds, 
	   ord_avg_over_seeds = ord_avg_over_seeds, objective_avg_over_seeds = objective_avg_over_seeds,
	   objective_stdev_over_seeds=objective_stdev_over_seeds,
	   purity_avg_over_seeds=purity_avg_over_seeds,purity_stdev_over_seeds=purity_stdev_over_seeds,
	   purity_lower_perc_over_seeds=purity_lower_perc_over_seeds,purity_upper_perc_over_seeds=purity_upper_perc_over_seeds)}
   data_obj$fres = forest_res_obj
	
	return(data_obj)
}

pred_new_data = function(data_obj_train=NULL,data_obj_test=NULL,plot=TRUE){
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj_train)){
      print("data_obj_test not specified; using default values")
      data_obj_train = read_data()
      data_obj_train = add_forest_to_obj(data_obj_train,log_weight_try=0)
   }
   if(is.null(data_obj_test)){
      print("data_obj_test not specified; using default values")
      data_obj_test = read_data(filename='./Data/GRB_short+outliers+Z_removed_reduced.arff')
   }
   ###########################################################################
	
	pred_vals = forest.pred(data_obj_train$forest,data_obj_test$features)
	data_obj_test$pred_vals = pred_vals
	
	if(!is.null(data_obj_test$triggerids)){
	   data_obj_test$pred_vals$triggerids = data_obj_test$triggerids
   }
	
	if(plot==TRUE){
	   Zlen_1 = length(pred_vals$alpha.hat) - 1
	   alpha_try_array = c(0:Zlen_1)/Zlen_1
	   ordered_alpha_hats = sort(pred_vals$alpha.hat)
	   frac_followed_up = ordered_alpha_hats*0.0
	   count = 0
	   for(alpha in alpha_try_array){
	      count = count+1
	      frac_followed_up[count] = sum(ordered_alpha_hats < alpha)/Zlen_1
	   }
	   
	   imagefile=paste('Plots/PopulationUnknown.pdf',sep="")
	   pdf(imagefile)
   	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), ylab=expression("Fraction of GRBs Followed Up (q^hat < q)"), xlab=expression('q'), pch="") # initialize plot))
      title(main=expression("Performance on bursts with unknown Z"), sub=data_obj_test$data_string)
      lines(alpha_try_array,alpha_try_array,lty=2,lwd=1)
      lines(alpha_try_array,frac_followed_up,lty=1,lwd=2)
      dev.off()
	}
	
	textfile='PopulationUnknown.txt'
	write.table(data_obj_test$pred_vals,textfile)
	
	return(data_obj_test)
}

add_forest_to_obj = function(data_obj=NULL,log_weight_try=2){
   # This function adds the full forest to the data object for later usage on 
   # new data.
    
   ##### If data object is not defined, create the default data object ######
   if(is.null(data_obj)){
      print("data_obj not specified; using default values")
      data_obj = read_data()
   }
   ###########################################################################
   if(!is.null(log_weight_try)){
      weight_try = 10^log_weight_try
      weights_vec = 1*(data_obj$data1$class == "low") + weight_try*(data_obj$data1$class == "high")
	}
	else{weights_vec=NULL}
	forest = forest.fit(data_obj$features,data_obj$classes,mtry=NULL,weights=weights_vec,n.trees=500,seed=sample(1:10^5,1))
   data_obj$forest = forest
   
   return(data_obj)
}



####### makes objective function plot ####### 
# forest_order is ordered using "order_residuals" function with high-z at low rank numbers
make_obj_fcn_plot = function(forest_order,data_obj=NULL,alpha_vec=seq(0.1,0.9,0.1),log_weights_try=seq(-1,1,0.2),imagefile="objective_fcn.pdf"){
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
	legend("bottomright", c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"), cex=0.6,col=col,lty=1:1)
	
	dev.off()
}

####### makes bumps plot, writes it to an image file, and saves the data that made it in a text file ####### 
make_bumps_plot = function(res_avg_over_seeds,data_obj=NULL,xlabel="log high-z weight",ylabel=expression(widehat(Q)),n_colors=128,z_width=3,imagefile="forest_probs_pred_bumps.pdf",textfile="forest_probs_pred.txt"){
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
   
   # Add the proper weight column names if we know them
   if(!is.null(data_obj$log_weights_try)){
      colnames(res_avg_over_seeds)=seq(-1,1,0.2)
   }
   
   pdf(file=imagefile,width=12,height=8) # save bumps plot
   
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(10,1), heights=c(2,2)) # make a separate plot for colorbar
par(mar=c(4,2,0,0))
   parcoord(res_avg_over_seeds,lwd=lwd_vec,var.label=FALSE,col=tc[col.vec])
   title(xlab=xlabel,cex.lab=1.25,mgp=c(2.5,1,0)) # axis labels
   title(ylab=ylabel,cex.lab=1.5,mgp=c(.25,1,0)) # yaxis
# Colorbar stuff
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


# fit and return forest on training data (using optimal tuning parameters)
cforest.fit = function(x,y,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
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
cforest.pred = function(forest,xnew){
  xnew = as.data.frame(xnew)
  n.new = dim(xnew)[1]
  n.old = length(predict(forest))
  # predictions for training data (to compute alpha-hats)
  pred.train = matrix(unlist(treeresponse(forest)),n.old,2,byrow=T) # CV this?
  # predict post. probs. for new data, with input forest
  predictions = matrix(unlist(treeresponse(forest,newdata=xnew)),n.new,2,byrow=T)
#  print(xnew$bat_is_rate_trig)[1]
  alpha.hat = NULL # compute alpha-hat values
  for(ii in 1:n.new){
    alpha.hat = c(alpha.hat, sum(predictions[ii,2]< pred.train[,2])/n.old)
  }

  return(list(alpha.hat = alpha.hat,prob.high=predictions[,2],prob.low=predictions[,1]))
}



# SAME AS cforest.fit but for the randomForest package
# fit and return forest on training data (using optimal tuning parameters)
forest.fit = function(x,y,mtry=NULL,weights=NULL,n.trees=500,seed=sample(1:10^5,1)){
  require(randomForest)
  set.seed(seed)

  # convert to randomForest weights
  # this solution is a bit hackish JPL
  if(!is.null(weights)){
    if(sum(weights != 1) > 0){
      high_weight = weights[weights != 1][1]
      weights = c(1,high_weight)
    }
    else {
      weights = c(1,1)
    }
  }
  
  n = length(y)
  p = length(table(y))
  if(is.null(mtry)){ # default for mtry
    mtry = ceiling(sqrt(dim(x)[2]))
  }
  train = cbind(y,x) # set up data to read into cforest
  names(train)[1] = "y"
  # fit random forest
  rf.fit = randomForest(y~.,data=train,
    classwt=c("low"=weights[1],"high"=weights[2]))
  return(rf.fit)
}



# USE FOR randomForest, use cforest.pred for party
# Predict redshift (alpha-hat) for test GRBs
# and returns probabilities of HIGH and LOW
# for test GRBs (currently not used)
forest.pred = function(forest,xnew){
  # xnew in correct format, should be already
  xnew = as.data.frame(xnew)
  n.new = nrow(xnew)
  n.old = length(forest$y)
  # predictions for training data and test
  pred.train = predict(forest,type="prob")
  predictions = predict(forest,newdata=xnew,type="prob")  
  # compute alpha.hats for the training
  alpha.hat = vapply(predictions[,2],
    function(x) { sum(x < pred.train[,2]) / n.old},0)
  # return everything as a list
  return(list(alpha.hat = alpha.hat,
              prob.high=predictions[,2],
              prob.low=predictions[,1]))
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

purity_vs_alpha = function(data_obj,weight_index=5,imagefile='test.pdf'){
   # Take the fifth weight for now 
   pdf(imagefile)
   
   high_cutoff = data_obj$high_cutoff
   newylab=paste("Percent of observed GRBs that are high z (z > ",high_cutoff,")",sep="")
   
   avg_pur = data_obj$fres$purity_avg_over_seeds[,weight_index]
   
   Zlen_1 = length(avg_pur) - 1
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   alpha_tries = seq(0,1,1/Zlen_1)
   base_purity = c(0:Zlen_1)*0+data_obj$num_high/(data_obj$num_low + data_obj$num_high)
   
   # avg_pur_high = data_obj$fres$purity_upper_perc_over_seeds[,weight_index]
   # avg_pur_low = data_obj$fres$purity_lower_perc_over_seeds[,weight_index]
   # 
   avg_pur_high = avg_pur + data_obj$fres$purity_stdev_over_seeds[,weight_index]
   avg_pur_low = avg_pur - data_obj$fres$purity_stdev_over_seeds[,weight_index]
   
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up (Normalized)"), ylab=newylab, pch="") # initialize plot))
   title(main=expression("Purity"), sub=data_obj$data_string)
   lines(alpha_tries, avg_pur, lty=1, lwd=4, col='red')
   xx = c(alpha_tries, rev(alpha_tries))
   yy = c(avg_pur_high, rev(avg_pur_low))
   col=rainbow(1)
   col_alpha = rainbow(1,alpha=0.3)
   polygon(xx,yy, col=col_alpha)
   lines(alpha_try_array,base_purity,lty=1,lwd=2)
   
   dev.off()
   
}



efficiency_vs_alpha = function(data_obj,weight_index=5,imagefile='test.pdf'){
   # Take the fifth weight for now 
   high_cutoff = data_obj$high_cutoff
   newylab=paste("Fraction of high (z > ",high_cutoff,") GRBs observed",sep="")
   
   avg_obj = data_obj$fres$objective_avg_over_seeds[,weight_index]
   avg_obj_high = avg_obj + data_obj$fres$objective_stdev_over_seeds[,weight_index]
   avg_obj_low = avg_obj - data_obj$fres$objective_stdev_over_seeds[,weight_index]
   
   Zlen_1 = length(avg_obj) - 1
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   alpha_tries = seq(0,1,1/Zlen_1)
   pdf(imagefile)
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up (Normalized)"), ylab=newylab, pch="") # initialize plot))
   title(main=expression("Efficiency"), sub=data_obj$data_string)
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   lines(alpha_tries, avg_obj, lty=1, lwd=4, col='red')
   xx = c(alpha_tries, rev(alpha_tries))
   yy = c(avg_obj_high, rev(avg_obj_low))
   col=rainbow(1)
   col_alpha = rainbow(1,alpha=0.3)
   polygon(xx,yy, col=col_alpha)
   
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   lines(alpha_tries, avg_obj, lty=1, lwd=4, col=col)
   
   dev.off()
   
}

efficiency_vs_purity = function(data_obj,weight_index=5,imagefile='test'){
   # Take the fifth weight for now 
   high_cutoff = data_obj$high_cutoff
   newxlab=paste("Fraction of high (z > ",high_cutoff,") GRBs observed",sep="")
   newylab=paste("Percent of observed GRBs that are high z (z > ",high_cutoff,")",sep="")
   
   avg_obj = data_obj$fres$objective_avg_over_seeds[,weight_index]
   avg_obj_high = avg_obj + data_obj$fres$objective_stdev_over_seeds[,weight_index]
   avg_obj_low = avg_obj - data_obj$fres$objective_stdev_over_seeds[,weight_index]
   
   avg_pur = data_obj$fres$purity_avg_over_seeds[,weight_index]
   
   
   Zlen_1 = length(avg_obj) - 1
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   alpha_tries = seq(0,1,1/Zlen_1)
   pdf(imagefile)
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=newxlab, ylab=newylab, pch="") # initialize plot))
   title(main=expression("Efficiency vs Purity"), sub=data_obj$data_string)
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   lines(avg_obj, avg_pur, lty=1, lwd=4, col='red')
   xx = c(alpha_tries, rev(alpha_tries))
   yy = c(avg_obj_high, rev(avg_obj_low))
   col=rainbow(1)
   col_alpha = rainbow(1,alpha=0.3)
#   polygon(xx,yy, col=col_alpha)
   
#   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   
#   dev.off()
   
}

multiple_efficiency_vs_alpha = function(data_obj_list,weight_index=1,ploterr=FALSE,imagefile='./Plots/ROC_multi.pdf'){
   pdf(imagefile)
   data_obj_1 = get(ls(data_obj_list)[1],pos=data_obj_list)
   high_cutoff = data_obj_1$high_cutoff
   newylab=paste("Fraction of high (z > ",high_cutoff,") GRBs observed",sep="")
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up (Normalized)"), ylab=newylab, pch="") # initialize plot))
   title(main=expression("Efficiency"))
   n_curves = length(data_obj_list)
   col = rainbow(n_curves)
   col_alpha = rainbow(n_curves,alpha=0.3)
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
      
      if(ploterr == TRUE){
         avg_obj_high = avg_obj + data_obj$fres$objective_stdev_over_seeds[,weight_index]
         avg_obj_low = avg_obj - data_obj$fres$objective_stdev_over_seeds[,weight_index]
         xx = c(alpha_tries, rev(alpha_tries))
         yy = c(avg_obj_high, rev(avg_obj_low))
         polygon(xx,yy, col=col_alpha[curve_index])         
      }
      curve_index = curve_index + 1
      name_list = c(name_list,data_obj$data_string)
   }
   # now print diagonal line for reference
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   
   print(name_list)
   legend("bottomright", name_list, cex=1.2,col=col,lty=1,lwd=4)
	
   dev.off()
}

multiple_efficiency_vs_alpha_weights = function(data_obj,weight_index_list=seq(1,11),ploterr=FALSE,imagefile='./Plots/ROC_multi_weights.pdf'){
   pdf(imagefile)
   data_obj_1 = data_obj
   high_cutoff = data_obj_1$high_cutoff
   newylab=paste("Fraction of high (z > ",high_cutoff,") GRBs observed",sep="")
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up (Normalized)"), ylab=newylab, pch="") # initialize plot))
   title(main=expression("Efficiency"))
   n_curves = length(weight_index_list)
   col = rainbow(n_curves)
   col_alpha = rainbow(n_curves,alpha=0.3)
   name_list = c()
   curve_index = 1
   for(weight_index in weight_index_list){
      print(paste('weight',weight_index))
      avg_obj = data_obj$fres$objective_avg_over_seeds[,weight_index]
      Zlen_1 = length(avg_obj) - 1
      print(Zlen_1)
      alpha_tries = seq(0,1,1/Zlen_1)
      lines(alpha_tries, avg_obj, lty=1, lwd=4, col=col[curve_index])
      
      if(ploterr == TRUE){
         avg_obj_high = avg_obj + data_obj$fres$objective_stdev_over_seeds[,weight_index]
         avg_obj_low = avg_obj - data_obj$fres$objective_stdev_over_seeds[,weight_index]
         xx = c(alpha_tries, rev(alpha_tries))
         yy = c(avg_obj_high, rev(avg_obj_low))
         polygon(xx,yy, col=col_alpha[curve_index])         
      }
      curve_index = curve_index + 1
      name_list = c(name_list,paste('weight',weight_index))
   }
   # now print diagonal line for reference
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   lines(alpha_try_array,alpha_try_array,lty=1,lwd=2)
   
   print(name_list)
   legend("bottomright", name_list, cex=1.2,col=col,lty=1,lwd=4)
	
   dev.off()
}

multiple_purity_vs_alpha_weights = function(data_obj,weight_index_list=seq(1,11),ploterr=FALSE,imagefile='./Plots/purity_multi_weights.pdf'){
   pdf(imagefile)
   data_obj_1 = data_obj
   high_cutoff = data_obj_1$high_cutoff
   newylab=paste("Percent of observed GRBs that are high z (z > ",high_cutoff,")",sep="")
	plot(x = c(0,1), y = c(0,1), xlim = c(0,1), ylim=c(0,1), xlab=expression("Fraction of GRBs Followed Up (Normalized)"), ylab=newylab, pch="") # initialize plot))
   title(main=expression("Purity"))
   n_curves = length(weight_index_list)
   col = rainbow(n_curves)
   col_alpha = rainbow(n_curves,alpha=0.3)
   name_list = c()
   curve_index = 1
   
   for(weight_index in weight_index_list){
      print(paste('weight',weight_index))
      avg_pur = data_obj$fres$purity_avg_over_seeds[,weight_index]
      Zlen_1 = length(avg_pur) - 1
      print(Zlen_1)
      base_purity = c(0:Zlen_1)*0+data_obj$num_high/(data_obj$num_low + data_obj$num_high)
      
      alpha_tries = seq(0,1,1/Zlen_1)
      lines(alpha_tries, avg_pur, lty=1, lwd=4, col=col[curve_index])
      
      if(ploterr == TRUE){
         avg_pur_high = avg_pur + data_obj$fres$purity_stdev_over_seeds[,weight_index]
         avg_pur_low = avg_pur - data_obj$fres$purity_stdev_over_seeds[,weight_index]
         xx = c(alpha_tries, rev(alpha_tries))
         yy = c(avg_pur_high, rev(avg_pur_low))
         polygon(xx,yy, col=col_alpha[curve_index])         
      }
      curve_index = curve_index + 1
      name_list = c(name_list,paste('weight',weight_index))
   }
   # now print diagonal line for reference
   alpha_try_array = c(0:Zlen_1)/Zlen_1
   lines(alpha_try_array,base_purity,lty=1,lwd=2)
   
   print(name_list)
   legend("topright", name_list, cex=1.2,col=col,lty=1,lwd=4)
	
   dev.off()
}

# Wrapper to make all representative plots for a given dataset
make_forest_plots = function(data_string="reduced",generate_data=FALSE, log_weights_try=seq(-1,1,0.2), Nseeds=10, roc_weight=11,redo_useless=FALSE, high_cutoff=4){
   # generate_data will re-do the smooth_random_forest_weights function, which takes a while
   data_filename = paste("./Data/GRB_short+outliers+noZ_removed_",data_string,".arff",sep="")
   data_results_dir = paste("./smooth_weights_results/smooth_weights_",data_string,"_",high_cutoff,sep="")
   obj_func_name = paste("./Plots/objective_fcn_",data_string,"_",high_cutoff,".pdf",sep="")
   bumps_pred_plot_name = paste("./Plots/forest_pred_bumps_",data_string,"_",high_cutoff,".pdf",sep="")
   bumps_rand_plot_name = paste("./Plots/forest_rand_bumps_",data_string,"_",high_cutoff,".pdf",sep="")
   bumps_plot_name = paste("./Plots/forest_order_bumps_",data_string,"_",high_cutoff,".pdf",sep="")
   bumps_pred_text_name = paste("forest_pred_bumps_",data_string,"_",high_cutoff,".txt",sep="")
   bumps_text_name = paste("forest_order_bumps_",data_string,"_",high_cutoff,".txt",sep="")
   roc_plot_name = paste("./Plots/ROC_",data_string,"_",high_cutoff,".pdf",sep="")
   purity_plot_name = paste("./Plots/Purity_",data_string,"_",high_cutoff,".pdf",sep="")
   multiple_purity_weight_plot_name = paste("./Plots/purity_multi_weights_",data_string,"_",high_cutoff,".pdf",sep="")
   multiple_efficiency_weight_plot_name = paste("./Plots/ROC_multi_weights_",data_string,"_",high_cutoff,".pdf",sep="")
   mydata = read_data(filename=data_filename,high_cutoff=high_cutoff)
   mydata$data_string = data_string
   mydata$log_weights_try=log_weights_try
   if(generate_data == TRUE){
      alphas.cv = smooth_random_forest_weights(data_obj = mydata,results_dir=data_results_dir,log_weights_try=log_weights_try,Nseeds=Nseeds,redo_useless=redo_useless)
   }
   mydata = extract_stats(data_obj = mydata, forest_res_dir=data_results_dir)
   purity_vs_alpha(mydata,imagefile=purity_plot_name,weight_index=roc_weight)
   efficiency_vs_alpha(mydata,imagefile=roc_plot_name,weight_index=roc_weight)
   if(redo_useless==FALSE){
      fres = mydata$fres$avg_over_seeds
      make_bumps_plot(fres, data_obj = mydata, imagefile=bumps_pred_plot_name,textfile=bumps_pred_text_name)
      fres_rand = noisify_residuals(fres)
      make_bumps_plot(fres_rand, data_obj = mydata, imagefile=bumps_rand_plot_name)
      fres_ordered = order_residuals(fres_rand)
      make_obj_fcn_plot(fres_ordered, data_obj = mydata,log_weights_try=log_weights_try, imagefile=obj_func_name)
      fres_ordered = order_residuals(fres_rand,reverse=TRUE)
      make_bumps_plot(fres_ordered, data_obj = mydata, ylabel=paste(expression(widehat(Q)),' (Normalized)'),imagefile=bumps_plot_name,textfile=bumps_text_name)
   }
   multiple_purity_vs_alpha_weights(data_obj=mydata,weight_index_list=seq(1,length(log_weights_try)),imagefile=multiple_purity_weight_plot_name)
   multiple_efficiency_vs_alpha_weights(data_obj=mydata,weight_index_list=seq(1,length(log_weights_try)),imagefile=multiple_efficiency_weight_plot_name)

   ## make_bumps_plot(alphas.cv)
 }
 
make_efficiency_plots = function(generate_data=FALSE, data_string_list=list('reduced','UVOTonly','UVOTandZpred','Nat_Zprediction','reduced_nozpredict','Full'), log_weights_try=seq(-1,1,0.2),roc_weight=11, ploterr=FALSE, high_cutoff=4){
   roc_plot_name = paste("./Plots/ROC_Multi.pdf",sep="")
   
   curve_index = 1
   data_obj_list = list()
   
   for(data_string in data_string_list){
      data_filename = paste("./Data/GRB_short+outliers+noZ_removed_",data_string,".arff",sep="")
      data_results_dir = paste("./smooth_weights_results/smooth_weights_",data_string,"_",high_cutoff,sep="")
      
      mydata = read_data(filename=data_filename,high_cutoff=high_cutoff)
      mydata$log_weights_try=log_weights_try
      mydata$data_string = data_string
      print(data_string)
      if(generate_data == TRUE){
         alphas.cv = smooth_random_forest_weights(data_obj = mydata,results_dir=data_results_dir,log_weights_try=log_weights_try)
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
   multiple_efficiency_vs_alpha(data_obj_list, weight_index=roc_weight,ploterr=ploterr)
}

# make_efficiency_plots(data_string_list=list('UVOTonly','UVOTonly_num_useless10','UVOTonly_num_useless30','UVOTonly_num_useless50','UVOTonly_num_useless70','UVOTonly_num_useless90'),log_weights_try=seq(0.8,1.0,0.2), roc_weight=1)
# make_efficiency_plots(data_string_list=list('UVOTonly','UVOTonly_num_useless1','UVOTonly_num_useless2','UVOTonly_num_useless4','UVOTonly_num_useless8','UVOTonly_num_useless16','UVOTonly_num_useless32'),log_weights_try=seq(-1.0,1.0,0.2), roc_weight=8)


make_all_plots = function(generate_data=FALSE,Nseeds=10,high_cutoff=4){
   make_forest_plots(data_string='reduced',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
   make_forest_plots(data_string='UVOTandZpred',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
   make_forest_plots(data_string='UVOTonly',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
 #  make_forest_plots(data_string='Nat_Zprediction',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
   make_forest_plots(data_string='reduced_nozpredict',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
   make_forest_plots(data_string='reduced_allzpredict',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)
   make_forest_plots(data_string='Full',generate_data=generate_data,Nseeds=Nseeds,high_cutoff=high_cutoff)

}

really_make_all_plots = function(){
 #  make_all_plots(generate_data=TRUE,Nseeds=100,high_cutoff=3.5)
   make_all_plots(generate_data=TRUE,Nseeds=100,high_cutoff=3.0)
}
make_all_useless_plots = function(generate_data=FALSE,Nseeds=10,log_weights_try=seq(-1.0,1.0,0.2),high_cutoff=4,roc_weight=1){
   make_forest_plots(data_string='UVOTonly_num_useless1',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)
   make_forest_plots(data_string='UVOTonly_num_useless2',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)
   make_forest_plots(data_string='UVOTonly_num_useless4',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)
   make_forest_plots(data_string='UVOTonly_num_useless8',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)
   make_forest_plots(data_string='UVOTonly_num_useless16',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)
   make_forest_plots(data_string='UVOTonly_num_useless32',log_weights_try=log_weights_try,generate_data=generate_data,Nseeds=Nseeds, roc_weight=roc_weight,redo_useless=TRUE)

}
