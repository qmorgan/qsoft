# fit a regression spline by minimizing GCV over the 
# degrees of freedom; use natural spline basis and return fit errors
library(splines)

splinetest = function(str='PROMPT_I',method='force'){
   #ignoring the first datapoint.. 
   if(str=='PROMPT_I'){
      maglist = c(15.642, 15.903, 16.259, 16.293, 16.386, 16.233, 16.253, 16.211, 16.166, 16.105, 16.144, 16.084, 15.995, 15.841, 15.979, 16.118, 16.243, 16.332, 16.459, 16.626, 16.882, 17.425, 17.755, 17.952, 18.26, 18.425, 18.644, 18.834, 18.946, 19.185, 19.389, 19.347, 20.008)
      magerrlist = c(0.061, 0.133, 0.175, 0.063, 0.069, 0.063, 0.112, 0.038, 0.042, 0.039, 0.036, 0.026, 0.024, 0.017, 0.018, 0.028, 0.02, 0.028, 0.028, 0.026, 0.031, 0.035, 0.044, 0.054, 0.061, 0.061, 0.08, 0.094, 0.115, 0.145, 0.169, 0.107, 0.317)
      tmidlist = c(62.208, 78.624, 93.312, 118.368, 146.88, 177.12, 202.176, 239.328, 288.576, 338.688, 388.80, 458.784, 550.368, 847.584, 1029.024, 1209.6, 1391.04, 1571.616, 1753.056, 1977.696, 2337.984, 3523.392, 4178.304, 4879.872, 5641.056, 6498.144, 7437.312, 8589.024, 9838.368, 11219.04, 12553.92, 14552.352, 16818.624)
   }
   
   if(str=='PAIRITEL_J'){
      maglist=c(13.6448, 14.2054, 14.5029, 14.7127, 14.4663, 14.5284, 14.4899, 14.4754, 14.3893, 14.2512, 14.23795, 14.203, 14.23225, 14.2286, 14.2353, 14.3229, 14.34645, 14.4171, 14.5291, 14.5703, 14.6332, 14.76545, 14.7765, 14.81475, 14.8298, 14.9631, 15.1151, 15.113, 15.3478, 15.6572, 15.7122, 15.7086, 16.14175, 16.5721, 16.84695, 17.347, 17.9644, 18.3544)
      magerrlist=c(0.0467224656447, 0.0796821931026, 0.0777198176215, 0.0862736900416, 0.0764255198218, 0.0743254394469, 0.0533391093159, 0.0525050750699, 0.0475727370972, 0.0448481794319, 0.0446858473088, 0.0413034177889, 0.0427611404985, 0.0432800556163, 0.0432035338533, 0.0444671158308, 0.0461222458038, 0.0479519872238, 0.0521054664421, 0.0507909913333, 0.0549084634088, 0.057527706502, 0.0590073469066, 0.0592619993242, 0.0639372897778, 0.0544676188237, 0.060448414827, 0.0607066742052, 0.0543182790159, 0.0699630052287, 0.101028993618, 0.0947614641131, 0.0917413582103, 0.0907513951089, 0.141951243185, 0.127075272216, 0.160167075531, 0.276049466344)
      tmidlist= c(65.5081150000032, 104.238114000192, 139.488200999712, 178.153240999968, 215.83830699993598, 253.578403000224, 311.12853000019203, 405.828696999744, 501.563774999808, 578.523967000416, 657.69901699968, 735.339140000352, 810.719209999584, 887.1046629984, 964.2394199980799, 1042.3645989984, 1119.7897390022401, 1195.7102140032, 1274.874822, 1351.5098870016, 1427.32004300064, 1504.4102400038398, 1582.0451889984, 1659.2253690009602, 1736.66542400352, 1870.30567799712, 2060.7510579984, 2250.71137799712, 2578.4669650032, 3038.3426560032, 3364.40309100192, 3556.40335699968, 3880.9639010016003, 5090.755774000319, 7186.378907001601, 10125.947663971201, 15645.7141150176, 22881.0319680096)
   }
   
   if(str=='PAIRITEL_K'){
      maglist=c(11.324, 11.6216, 11.86225, 11.8657, 12.0163, 12.0487, 12.1347, 12.05075, 12.002, 11.9073, 11.8568, 11.8904, 11.8885, 11.952, 11.9561, 11.984, 12.1988, 12.3072, 12.3329, 12.3852, 12.40845, 12.46195, 12.6021, 12.6171, 12.6581, 12.692, 12.75125, 12.95375, 13.0718, 13.4646, 13.6336, 13.4988, 13.92465, 14.29325, 14.5291, 15.2455, 15.60365)
      magerrlist=c(0.0504600537282, 0.0641346365201, 0.0622045384518, 0.0590804508857, 0.063505866629, 0.0633843446928, 0.0490441372373, 0.04848058265, 0.0454873238284, 0.0456998612331, 0.0459437317991, 0.0443466881129, 0.0435791155053, 0.0472164511477, 0.0452409794978, 0.04478503843, 0.0494026202085, 0.0529472096362, 0.0575971423091, 0.0526224549169, 0.0544129381042, 0.0560037503704, 0.0597131251192, 0.0620757043101, 0.0651452256356, 0.0517351053514, 0.0526948623281, 0.0595292176535, 0.054695103153, 0.0714847964689, 0.101420986407, 0.088968712733, 0.0876827734703, 0.083793101736, 0.125934668153, 0.129540540017, 0.119588108183)
      tmidlist=c(65.5081150000032, 104.238114000192, 139.488200999712, 178.153240999968, 215.83830699993598, 253.578403000224, 311.12853000019203, 405.828696999744, 501.563774999808, 578.523967000416, 657.69901699968, 735.339140000352, 810.719209999584, 887.1046629984, 964.2394199980799, 1042.3645989984, 1119.7897390022401, 1195.7102140032, 1274.874822, 1351.5098870016, 1427.32004300064, 1504.4102400038398, 1582.0451889984, 1659.2253690009602, 1736.66542400352, 1870.30567799712, 2060.7510579984, 2250.71137799712, 2578.4669650032, 3038.3426560032, 3364.40309100192, 3556.40335699968, 3880.9639010016003, 5090.755774000319, 7186.378907001601, 10125.947663971201, 15645.7141150176)
   }
     
   # should be doing the spline fit in log space
   logtmidlist = log10(tmidlist)
   logouttmid=seq(1.8,4.1,0.01)
   ret = regrSplineTextout(logtmidlist,maglist,magerrlist,logouttmid,maxdf=20,method=method)
   
   # not sure maxdf is doing what i intended.. 
   plot(logtmidlist,maglist*-1)
   points(logtmidlist,maglist*-1-magerrlist,pch=20)
   points(logtmidlist,maglist*-1+magerrlist,pch=20)
   
   lines(logouttmid,ret$fit*-1)
   lines(logouttmid,ret$fit*-1-ret$fiterr)
   lines(logouttmid,ret$fit*-1+ret$fiterr)
   
   return(ret)
   
}

# setting up store environment
q_dir = Sys.getenv('Q_DIR')
store_dir = paste(q_dir,'store/',sep='')

regrSplineTextout = function(maxdf=floor(length(x)/2),method="gcv",outbase='regrSpline'){
   
   # Reads in text outputted to the store directory 
   
   magpath=paste(store_dir,'regrSpline_y_in.txt',sep='')
   toutpath=paste(store_dir,'regrSpline_x_in.txt',sep='')
   
   magmat = as.matrix(read.table(magpath),header=FALSE)
   x = magmat[,1] #tmids
   y = magmat[,2] #magnitudes
   yErr = magmat[,3] #magerrs 
   
   toutmat = as.matrix(read.table(toutpath),header=FALSE)
   xgrid = toutmat[,1] #desired output times
   
   fitlist = regrSpline(x,y,yErr,xgrid,maxdf=maxdf,method=method)
   
   fitfile=paste(store_dir,outbase,'_fit.txt',sep='')
	fiterrfile=paste(store_dir,outbase,'_fiterr.txt',sep='')
	
	write.table(fitlist$fit,fitfile)
	write.table(fitlist$fiterr,fiterrfile)
   #return(fitlist)
}


regrSpline = function(x,y,yErr,xgrid,maxdf=floor(length(x)/2),method="gcv"){
  ####
  #### INPUT ####
  # x - input array of x values (predictors)
  # y - input array of y values (response)
  # yErr - array of errors on input y values
  # xgrid - grid of x values for output fit
  # maxdf - maximum number of spline degrees of freedom (knots)
  # method - method to choose optimal spline fit. Either "gcv" or "BIC" or "force"
   # "force" just uses the maxdf as the # of df in the spline
  ####
  #### OUTPUT ####
  # fit - optimal spline model, evaluated on xgrid
  # fiterr - error on model output, evaluated on xgrid
  # df - optimal spline degrees of freedom
  # score - array of GCV or BIC scores
  
  n = length(x)
  
  score = NULL
  # max 50 splines
  df = 2:min(50,max(2,maxdf)) # grid of df to minimize over
  for(ii in 1:length(df)){
    # spline basis at each df
    #Generate a Basis Matrix for Natural Cubic Splines
    #Generate the B-spline basis matrix for a natural cubic spline. Usage: 
    # ns(x, df = NULL, knots = NULL, intercept = FALSE,
    #         Boundary.knots = range(x))
    spbasis = ns(x,df=df[ii])

    # fit spline model
    # ‘lm’ is used to fit linear models.  It can be used to carry out
    # regression, single stratum analysis of variance and analysis of
    # covariance (although ‘aov’ may provide a more convenient interface
    # for these).
    spfit = lm(y ~ spbasis,weights=1/yErr^2)

    if(method=="gcv"){
      # compute generalized cross-validation score
      L =  spbasis %*% solve(t(spbasis)%*%spbasis)%*%t(spbasis)
      nu = sum(diag(L))
      score = c(score,1/n * sum(((y-predict(spfit))/(yErr*(1-nu/n)))^2))
    }
    
    if(method=="BIC"){
      # compute Bayesian Information Criterion
      score = c(score, sum(log(2*pi*yErr^2) + (y-predict(spfit))^2/(2*yErr))+length(spfit$coefficients)*log(n))
    }
  }

  # determine optimal model (min GCV or min BIC)
  if(method=="gcv"){
    dfopt = df[which.min(score)]
    print(paste("Optimal # of splines:",dfopt))
  }
  if(method=="BIC"){
    dfopt = df[which.min(score)]
  }
  if(method=='force'){
     dfopt= maxdf
  }
  # regress on the optimal model
  fit = lm(y ~ ns(x,dfopt),weights=1/yErr^2)
  # generate predictions and standard errors for optimal model on xgrid
  ygrid = predict(fit, data.frame(x=xgrid),se.fit=T)
  if(n==3){
    ygrid$se.fit = rep(max(yErr),length(xgrid))
  }
  
  return(list(fit = ygrid$fit, fiterr = ygrid$se.fit, df=dfopt, score = score))
}

