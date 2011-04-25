########## 
########## 
########## COMPUTE FEATURE IMPORTANCE FROM RF CLASSIFIER 
########## USES FUNCTIONS IN grb_class.R - PRETTY SIMPLE 
##########
########## by James Long 
########## date: 04/24/2011 
########## 


source('grb_class.R')

data1 = read_data()

rf.fit = randomForest(class~.,data=data1$data1,
  strata=data1$data1$class,sampsize=c(8,8),importance=TRUE)

rf.fit
importance(rf.fit)
varImpPlot(rf.fit,type=1)


plot(data1$data1$PROB_Z_GT_4,(as.numeric(data1$data1$uvot_detection) + rnorm(nrow(data1$data1),sd=.1)),col=data1$data1$class,pch=as.numeric(data1$data1$class))


table(data1$data1$PROB_Z_GT_4 > .13,data1$data1$class)
table(data1$data1$uvot_detection,data1$data1$class)

R> library(randomForest)
randomForest 4.5-23
Type rfNews() to see new features/changes/bug fixes.
R> f <- factor(sample(1:4, nrow(iris), replace=TRUE)) 
R> rf1 <- randomForest(iris[1:4], iris[[5]], strata=f, sampsize=rep(5, nlevels(f)))
R> rf1 
