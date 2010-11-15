data1 = list()

####
#### 1st set of features
####

data1[[1]] = read.arff("070710_shortremoved_NoZremoved_outremoved_bat_prompt.arff")
data1[[1]] = removeErrors(data1[[1]])
data1[[1]] = cleanData(data1[[1]],4)
data1[[1]]$triggerid_str = NULL


####
#### 2nd set of features
####

data1[[2]] = read.arff("070710_shortremoved_NoZremoved_outremoved_nfi_prompt.arff")
data1[[2]] = removeErrors(data1[[2]])
data1[[2]] = cleanData(data1[[2]],4)
data1[[2]]$triggerid_str = NULL




####
#### 3rd set of features
####

data1[[3]] = read.arff("070710_shortremoved_NoZremoved_outremoved_nfi_prompt.arff")
data1[[3]] = removeErrors(data1[[3]])
data1[[3]] = cleanData(data1[[3]],4)
data1[[3]]$triggerid_str = NULL





####
#### 4 set of features
####

data1[[4]] = read.arff("070710_shortremoved_NoZremoved_outremoved_late_proccessed.arff")

# make data1 look nice
data1[[4]] = removeErrors(data1[[4]])
data1[[4]] = cleanData(data1[[4]],4)
data1[[4]]$triggerid_str = NULL
data1[[4]]$CHI2_PC = NULL
data1[[4]]$CHI2_WT = NULL
data1[[4]]$CHI2_PC_LATE = NULL







# levels of prior
PRIOR = (1:19) / 20
trees = list()

for(j in 1:length(data1)){
	for(i in 1:length(PRIOR)){
      		fit1 = rpart(class ~ .,parms=list(prior=c(PRIOR[i],1-PRIOR[i])),method="class",data = data1[[j]])
      		if(length(fit1$cptable) > 3){
        		bestRow = which.min(fit1$cptable[,4])
        		cp = fit1$cptable[bestRow,1]
        		fit1 = prune(fit1,cp=cp)
      		}	
      		trees[[j]][[i]] = fit1
	}
}




