#############
############# these pieces of code show how the loss and priors are related in rpart
#############

# p 112 from Breiman "The priors can be used to adjust the individual misclassification cost in any desired direction." - this is exactly what we want to do

set.seed(250)
nobs = 500
X1 = runif(nobs)
X2 = runif(nobs)
classes = factor(rep("1",nobs),levels=c("1","2"))
classes[runif(nobs) > X2] = "2" # X2 is prob of class 2

datatest = as.data.frame(cbind(X1,X2))
datatest$class = classes


######### now for the demo

fit0 = rpart(class ~ X1 + X2,method='class',data=datatest)

# fit 1 and fit 2 are the same because of how their priors / loss interact
prior = c(.5,.5)
fit1 = rpart(class ~ X1 + X2,method='class',parms=list(prior=prior),data=datatest)

prior = c(.8,.2)
loss = matrix(c(0,1,4,0),ncol=2)
fit2 = rpart(class ~ X1 + X2,method='class',parms=list(prior=prior,loss=loss),data=datatest)



##########
########## an example with the data
##########


set.seed(250)
library('foreign')
library('rpart')
library('plotrix')
library('xtable')



##
## get the data in the desired form
##
data1 = read.arff('070710_shortremoved_NoZremoved.arff')
Z = data1$Z
data1 = removeErrors(data1)
data1 = cleanData(data1,4)

data1 = subset(data1,select=(!(names(data1) %in% c("CHI2_PC","CHI2_WT","CHI2_PC_LATE"))))




##
## fit1 and fit2 are the same because priors and loss work out to give the same
## pi twidle ( see page 9 of guide to rpart)
##

fit0 = rpart(class~.,method='class',data=data1)

prior = c(.5,.5)
fit1 = rpart(class~.,method='class',parms=list(prior=prior),data=data1)

prior = c(.8,.2)
loss = matrix(c(0,1,4,0),ncol=2)
fit2 = rpart(class ~ .,method='class',parms=list(prior=prior,loss=loss),data=data1)


prior = c(.8,.2)
fit3 = rpart(class~.,method='class',parms=list(prior=prior),data=data1)





