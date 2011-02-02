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
newfeatures = newdata$features
nbursts = dim(newfeatures)[1]
# this is the fix: have to loop thru the columns and ensure that the
# factors have the same levels as the original data...ugh!
for(jj in 1:dim(newfeatures)[2]){
  if(is.factor(newfeatures[,jj])){
    cat("Feature",names(newfeatures)[jj],"is factor\n")
    newfeatures[,jj] = factor(newfeatures[,jj],levels=levels(features[,jj]))
  }
}

pred.val = forest.pred(rf,newfeatures)

