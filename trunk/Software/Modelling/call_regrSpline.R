library(splines)

# setting up store environment
q_dir = Sys.getenv('Q_DIR')
load_dir = paste(q_dir,'trunk/Software/Modelling/regrSpline.R',sep='')

source(load_dir)

# setting up store environment
q_dir = Sys.getenv('Q_DIR')
store_dir = paste(q_dir,'store/',sep='')

regrSplineTextout()
