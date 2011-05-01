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

bs_features = matrix(rnorm(nrow(data1$data1)*20),ncol=20)

data1$data1 = data.frame(data1$data1,bs_features)
dim(data1$data1)

data1$data1$PROB_Z_GT_4 = NULL

rf.fit = randomForest(class~.,data=data1$data1,
  strata=data1$data1$class,sampsize=c(7,7),importance=TRUE)


  # change uvot to numeric
  data1$data1$uvot_detection = as.numeric(
    data1$data1$uvot_detection) + rnorm(nrow(data1$data1),sd=.0001)


rf.fit
importance(rf.fit)
varImpPlot(rf.fit,type=1)


plot(data1$data1$PROB_Z_GT_4,(as.numeric(data1$data1$uvot_detection) + rnorm(nrow(data1$data1),sd=.1)),col=data1$data1$class,pch=as.numeric(data1$data1$class))


table(data1$data1$PROB_Z_GT_4 > .13,data1$data1$class)
table(data1$data1$uvot_detection,data1$data1$class)


  data1 = read_data()
  rf.fit = randomForest(class~.,data=data1$data1,
    classwt=c("low"=1,"high"=1000),importance=TRUE)
  rf.fit



# get the data
high_weights = c(.5,1,5,20,50,100)
for(high_weight in high_weights){
  data1 = read_data()
  bs_features = matrix(rnorm(nrow(data1$data1)*20),ncol=20)
  data1$data1 = data.frame(data1$data1,bs_features)
  dim(data1$data1)

  # using weights
  pdf(paste('RFIMP',high_weight,'.pdf',sep=""))
  par(mfcol=c(1,2))
  rf.fit = randomForest(class~.,data=data1$data1,
    classwt=c("low"=1,"high"=high_weight),importance=TRUE)
  rf.fit
  importance(rf.fit)
  varImpPlot(rf.fit,type=2,main=paste("UVOT Discrete",high_weight))

  # change uvot to numeric
  data1$data1$uvot_detection = as.numeric(
    data1$data1$uvot_detection) + rnorm(nrow(data1$data1),sd=.0001)


  rf.fit = randomForest(class~.,data=data1$data1,
    classwt=c("low"=1,"high"=high_weight),importance=TRUE)
  importance(rf.fit)
  varImpPlot(rf.fit,type=2,main=paste("UVOT Continuous",high_weight))
  dev.off()
}






R> library(randomForest)
randomForest 4.5-23
Type rfNews() to see new features/changes/bug fixes.
R> f <- factor(sample(1:4, nrow(iris), replace=TRUE)) 
R> rf1 <- randomForest(iris[1:4], iris[[5]], strata=f, sampsize=rep(5, nlevels(f)))
R> rf1 


