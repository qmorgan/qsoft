## tests how performance changes as more useless features are added to data
##
## by James Long
## date Dec 22, 2010
##

library('foreign')

source('grb_class.R')
source('algorithm1/algorithm1.R')

filename='algorithm1/uvot_no_error.arff'
data1 = read.arff(filename)
Z = data1$Z
data1 = cleanData(data1,4) # define above 4 as high, below 4 as low
data1 = subset(data1,select=((names(data1) %in% c("class","uvot_detection"))))


## get some numbers about the data
print('number of features being used:')
print(length(data1) - 1)
print('number of observations:')
print(nrow(data1))
print('number of high grbs:')
print(sum(data1$class=="high"))
print('number of low grbs:')
print(sum(data1$class=="low"))
print('names of features:')
print(names(data1))
table(data1$class,data1$uvot_detection)


implement(data1,Z)


####
#### have to define inputs for forest_run
####
weight = 10^5
classes = data1$class
high_cutoff = 4
num_high = sum(classes == "high")
num_low = sum(classes == "low")
useless_features = matrix(rnorm(length(Z)*10),ncol=10)
features = data.frame(data1$uvot_detection,useless_features)




useless = 10:70
result = rep(0,length(useless))
for(i in 1:length(result)){
  print(paste(i,"of",length(result),"runs"))
  useless_features = matrix(rnorm(length(Z)*useless[i]),ncol=useless[i])
  features = data.frame(data1$uvot_detection,useless_features)
  output = forest_run()
  result[i] = output$high_as_high / (output$high_as_high + output$high_as_low)
}



pdf('performance_vs_useless_features_log.pdf')
linear.fit = lm(result ~ useless)
plot(useless,result,ylim=c(0,max(max(result),13/15)),xlab="Number of Useless Features",ylab="Fraction of Highs Classified as High",main="RF Performance vs. Number of Useless Features: Weight 7 on High")
abline(h=13/17,col='grey')
abline(linear.fit,col='grey')
dev.off()

