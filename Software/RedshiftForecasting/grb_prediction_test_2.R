# To Load new GRBs into the classifier and test them:
library('foreign')
source('grb_class.R')
source('algorithm1/algorithm1.R')

data_obj = read_data()
features = data_obj$features
classes = data_obj$classes
weight = 10^1 # how much are high bursts weighted
wts = ifelse(classes=="high",weight,1)
rf = forest.fit(features,classes,mtry=NULL,weights=wts,n.trees=500)
newdata = read_data(filename='./algorithm1/testdata.arff')
nbursts = dim(newdata$features)[1]
allfeatures = rbind(newdata$features,features)
newfeatures = data.frame(allfeatures[1:nbursts,])
pred.val = forest.pred(rf,newfeatures)