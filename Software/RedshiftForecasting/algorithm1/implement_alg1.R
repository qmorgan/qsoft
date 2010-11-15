########
########
######## Get GRB data in form so we can run through algorithm1.R
########
########
########
set.seed(250)
library('foreign')
library('rpart')
library('plotrix')
library('xtable')
source('algorithm1.R')


##
## get the data in the desired form
##
filename = '070710_shortremoved_NoZremoved_outremoved_late_proccessed.arff'
data1 = read.arff(filename)
Z = data1$Z
data1$triggerid_str = NULL # meaningless feature in the data sets
data1 = removeErrors(data1)
data1 = cleanData(data1,4) # define above 4 as high, below 4 as low

# these features are removed because they are formatted weirdly, according to Adam they
# probably are not important in prediction anyway
data1 = subset(data1,select=(!(names(data1) %in% c("CHI2_PC","CHI2_WT","CHI2_PC_LATE"))))

# uncomment the following line if want code to run fast (vastly reduced # of features)
# this is used mostly for testing
# data1 = subset(data1,select=((names(data1) %in% c("class","A","FLX_PC_LATE","wh_mag_isupper"))))
data1 = subset(data1,select=((names(data1) %in% c("class","wh_mag_isupper"))))



# features names Adam thinks are important are stored in vector good_features
# recall there are 4 sets of features we get, say F1,F2,F3, and F4 where F1<F2<F3<F4
# (we get more features as time passes). if we are analyzing feature
# set F3 we can analyze the features that are in F3 and ''good_features''
#good_features = c("class","A","B","EP0","FL","FLX_PC_LATE","GAM_PC","MAX_SNR","NH_PC","T90","bat_image_signif","bat_image_peak","bat_trigger_dur","bat_is_rate_trig","v_mag_isupper")
# uncomment line below if you only want to use Adam's recommended feature set
#data1 = subset(data1,select=((names(data1) %in% good_features)))


###
### Important: data1 should now be dataframe with first column ''class'' a factor with two 
### levels. Remaining columns are features, which may be continuous and/or categorical, 
### missingness okay. The actual redshift (numerical value) should be stored in the variable
### Z. The order of Z should match the observation order in data1.
###


###
### some info about this run for the user to see
###
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





# get a heat map
#heatmap(data1)
# get everything else




implement(data1,Z)






prior_alpha1 = algorithm1(data1,10,(1:9) / 10,.1)
prior_alpha3 = algorithm1(data1,10,(1:9) / 10,.3)
prior_alpha5 = algorithm1(data1,10,(1:9) / 10,.5)


#### generate trees on entire data using a few different priors, if the tree is 
#### this is optional and separate for all the machinery of algorithm1
fit1 = rpart(class ~ .,parms=list(prior=c(1-prior_alpha1,prior_alpha1)),method="class",data = data1)
fit2 = rpart(class ~ .,parms=list(prior=c(1-prior_alpha3,prior_alpha3)),method="class",data = data1)
fit3 = rpart(class ~ .,parms=list(prior=c(1-prior_alpha5,prior_alpha5)),method="class",data = data1)


par(mfcol=c(1,1))
pdf('alpha_1.pdf')
plot(fit1,margin=.1)
text(fit1,use.n=T,pretty=F)
dev.off()

pdf('alpha_3.pdf')
plot(fit2,margin=.1)
text(fit2,use.n=T,pretty=F)
dev.off()

pdf('alpha_5.pdf')
plot(fit3,margin=.1)
text(fit3,use.n=T,pretty=F)
dev.off()





if(length(fit1$cptable) > 3){
	bestRow = which.min(fit1$cptable[,4])
	cp = fit1$cptable[bestRow,1]
	fit1 = prune(fit1,cp=cp)
}
fit1 = rpart(class ~ .,parms=list(prior=c(.7,.3)),method="class",data = data1)
if(length(fit1$cptable) > 3){
	bestRow = which.min(fit1$cptable[,4])
	cp = fit1$cptable[bestRow,1]
	fit1 = prune(fit1,cp=cp)
}







