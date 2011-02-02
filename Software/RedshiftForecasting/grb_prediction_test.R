## tests grb prediction algorithm
##
## by JWR
## created Feb 1, 2011
##

# remove everything from R's memory
rm(list=ls(all=TRUE))

library('foreign')
source('grb_class.R')
source('algorithm1/algorithm1.R')


# parameters to random forest function
weight = 10^1 # how much are high bursts weighted
alpha = .3 # fraction to follow up on

data_obj = read_data()
features = data_obj$features
classes = data_obj$classes
wts = ifelse(classes=="high",weight,1)

# build random forest on all features
rf = forest.fit(features,classes,mtry=NULL,weights=wts,n.trees=500)
# note: you can save using save(), load using load()
# predict for new data (here new data is just the old data)
pred = forest.pred(rf,features)

##### do a test for 25 "new" bursts
validate = 1:25

rf = forest.fit(features[-validate,],classes[-validate],mtry=NULL,weights=wts[-validate],n.trees=500)
pred.val = forest.pred(rf,features[validate,])

pdf("grb_pred_test.pdf")
plot(pred.val$alpha.hat,pred.val$prob.high,col=classes[validate],pch=19)
dev.off()
