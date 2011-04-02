#######
####### rewrite of forest.pred for implementation
####### using randomForest, had to remove a few
####### functions associated with the party
####### package such as 'treeresponse'
#######
####### by James Long
####### date April 2, 2011
#######


library('randomForest')


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



####
#### test forest.pred on iris dataset
####

## only use 2 classes . . .
iris = iris[1:100,]
iris$Species = as.factor(as.character(iris$Species))
# set up test and training
test = (1:nrow(iris)) %in% sample(1:nrow(iris),10)
train = !test
iris.train = iris[train,]
iris.test = iris[test,1:4]
# construct classifier
rf = randomForest(Species ~ .,data=iris.train)
# get alpha-hats and probability of highs
alpha.hats = forest.pred(rf,iris.test)


