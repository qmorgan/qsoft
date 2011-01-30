## tests how performance changes as more useless features are added to data
##
## by James Long
## created Dec 22, 2010, updated Jan 28, 2011
##

# remove everything from R's memory
rm(list=ls(all=TRUE))

library('foreign')
source('grb_class.R')
source('algorithm1/algorithm1.R')


# parameters to random forest function
weight = 10^1 # how much are high bursts weighted
alpha = .3 # fraction to follow up on

# numbers of useless features to try
useless = 1:50

data_obj = read_data()
features = data_obj$features
result = rep(0,length(useless))

### so the good features we are using are, could try
### replacing these with just uvot
names(data_obj$features)

### at each iteration we add 1 more useless feature
### and output high as high / (high as high + high as low)
### to the results vector
for(i in 1:length(result)){
  print(paste(i,"of",length(result),"runs"))
  useless_features = matrix(rnorm(length(data_obj$Z)*useless[i]),ncol=useless[i])
  data_obj$features = data.frame(features,useless_features)
  output = forest_run(data_obj=data_obj,weight=weight,alpha=alpha)
  result[i] = output$high_as_high / (output$high_as_high + output$high_as_low)
}


### plot the results to see how performance degrades
plot_title = paste('performance_vs_useless_features_weight_',weight,'.pdf',sep="")
pdf(plot_title)
linear.fit = lm(result ~ useless)
main_title = paste("RF Performance vs. Number of Useless Features: Weight",weight,"on High")
plot(useless,result,ylim=c(0,max(max(result),13/15)),xlab="Number of Useless Features",ylab="Fraction of Highs Classified as High",main=main_title)
abline(h=13/17,col='grey') # if we follow up on all uvot no, alpha > 49 / 151 in order to do this
abline(linear.fit,col='grey')
dev.off()

