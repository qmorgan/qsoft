###
###
### for testing algorithm1
### 
###

source('algorithm1.R')

###
### Important: data1 should now be dataframe with first column ''class'' a factor with two 
### levels. Remaining columns are features, which may be continuous and/or categorical, 
### missingness okay. The actual redshift (numerical value) should be stored in the variable
### Z. The order of Z should match the observation order in data1.
###




##
## functions for generating observations, will be called by later sections of code
##

# assists generateObs in generating observations from the clean training distribution
getCoord = function(point){
	 uniDraw = runif(1)
	 obs = c(0,0,0)
	 if(point == 1){
	 	  obs[1] = 1*(uniDraw < .9)
		  obs[2] = runif(1)
		  obs[3] = 1 + runif(1)
	 }
	 if(point == 2){
	 	  obs[1] = 1*(uniDraw < .1)
		  obs[2] = 1 + runif(1)
		  obs[3] = 1 + runif(1)
	 }
	 if(point == 3){
	 	  obs[1] = 1*(uniDraw < .5)
		  obs[2] = 2 * runif(1)
		  obs[3] = runif(1)
	 }
	 return(obs)
}


# generates observations according to desired distribution
generateObs = function(nobs){
	    partition = sample(c(1,2,3),size=nobs,replace=T,prob=c(.45,.45,.1))
	    data1 = vapply(partition, getCoord,c(0,0,0))
	    data1 = t(data1)
	    data1 = as.data.frame(data1)
	    names(data1) = c("Class","X","Y")
	    class1 = rep("Black",nobs)
	    class1[data1[,1] == 1] = "Red"
	    data1$Class = class1
	    return(data1)
}





# generates observations according to desired distribution
generateObsFrac = function(nobs){
	class_levels = as.character(c('low','high'))
	class = factor(rep("high",nobs),levels=class_levels)
	X1 = runif(nobs)
	X2 = runif(nobs)
	success = runif(nobs)
	class[success > X1] = "low"
	data1 = as.data.frame(class)
	data1$X1 = X1
	data1$X2 = X2
	return(data1)
}


data1 = generateObsFrac(10000)


colors1 = rep('black',nrow(data1))
colors1[data1$class == "high"] = "red"
par(mfcol=c(1,1))
plot(data1$X1,data1$X2,col=colors1,xlab="X",ylab="Y")





##
## should run algorithm1 on data and then make the tree
##
# heatMap(data1)
# implement(data1)




prior_chosen = algorithm1(data1,10,(1:9) / 10,.9)
print(prior_chosen)


prior_chosen = .9
fit1 = rpart(class ~ .,parms=list(prior=c(1-prior_chosen,prior_chosen)),method="class",data=data1)


par(mfcol=c(1,2))
plot(data1$X1,data1$X2,col=colors1,xlab="X",ylab="Y")
plot(fit1,margin=.1)
text(fit1)









##
## simulate with 2 features, 1 useless, the other the probability of high, see how
## line changes (mostly test algorithm 1's performance)
##





##
## simulate with p features, p-1 are useless, 1 is binary and quite useful, 120 obs, GRB like
##
