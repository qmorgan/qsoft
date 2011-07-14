########## 
########## 
########## CODE FOR INTERPRETING HOW RF IS WORKIN 
########## ANAYSIS ON NO REDSHIFT DATA AND TRAINING 
##########
########## by James Long 
########## date: 5/23/2011 
########## 


rm(list = ls(all = TRUE))
library('rpart')
library("RColorBrewer")
source('grb_class.R')

## get alpha hats for no redshift data
calib = read.table('Calib_testdata.txt')
summary(calib)

## get features for no redshift data
data1 = read_data('Data/GRB_short+outliers+Z_removed_reduced.arff')
summary(data1)
data1$Z

## merge alpha hats with features
head(data1$features)
data1$features$triggerids = data1$triggerids

## 212 new bursts with all features
alpha.hat.features = merge(data1$features,calib)
nrow(alpha.hat.features)
head(alpha.hat.features)

alpha.hat.features$uvot_detection = as.numeric(alpha.hat.features$uvot_detection) - 1


tc = two.colors(n=nrow(alpha.hat.features),start="red",end="blue",middle="white")
scatter = rnorm(nrow(alpha.hat.features),sd=.1)
ys = as.numeric(alpha.hat.features$uvot_detection) + scatter
xs = alpha.hat.features$PROB_Z_GT_4
ymax = max(ys) * 1.02
ymin = min(ys) * .95
xmax = max(xs) * 1.02
xmin = min(xs) * .95


## look at how uvot and nat's feature determine q.hat
pdf('Plots/uvotNatZcolorByQhat.pdf')
plot(xs,ys,ylim=c(ymin,ymax),xlim=c(xmin,xmax),ylab="uvot (0=no, 1=yes)",xlab="PROB_Z_GT_4",col=tc[rank(alpha.hat.features$alpha.hat)])
colorbar.plot(.01 + xmax,(ymax + ymin)/2,strip=seq(min(alpha.hat.features$alpha.hat),max(alpha.hat.features$alpha.hat),length.out=length(tc)),col=tc,horizontal=FALSE,strip.width=.02,strip.length=1) # 
text(xmax-.02,ymax,paste("q.hat=",round(max(alpha.hat.features$alpha.hat),2)))
text(xmax-.02,ymin,paste("q.hat=",round(min(alpha.hat.features$alpha.hat),2)))
dev.off()

alpha.hat.features[order(alpha.hat.features$alpha.hat),c("alpha.hat","uvot_detection","PROB_Z_GT_4")]


pdf('Plots/uvotQhat.pdf')
plot(alpha.hat.features$alpha.hat,alpha.hat.features$uvot_detection + rnorm(nrow(alpha.hat.features),sd=.1),xlab="q hat",ylab="uvot")
dev.off()

d1 = density(alpha.hat.features$alpha.hat[alpha.hat.features$uvot_detection == 0])
d2 = density(alpha.hat.features$alpha.hat[alpha.hat.features$uvot_detection == 1])
ymax = max(d1$y,d2$y)

pdf('Plots/qHatdensityUVOTyesUVOTno.pdf')
plot(d1,col='blue',lty=2,xlab="q hat",main="")
lines(d2,col='orange')
legend(-.1,1.5,c("uvot=NO","uvot=yes"),col=c("blue","orange"),lty=c(2,1))
dev.off()

## construct CART tree on new and old data
feat.noz = data1$features
feat.noz$class = "noZ"
feat.noz$triggerids = NULL

data.train = read_data()
feat.z = data.train$features
feat.z$class = "Z"

both = rbind(feat.noz,feat.z)
both$class = as.factor(both$class)

tree.hasz = rpart(class ~ .,data=both)
tree.hasz
plot(tree.hasz,margin=.1)
text(tree.hasz,use.n=TRUE)

plot(both$B,both$A,col=both$class)



both$A = NULL
both$B = NULL
tree.hasz = rpart(class ~ .,data=both)
tree.hasz
plot(tree.hasz,margin=.1)
text(tree.hasz,use.n=TRUE)




