import numpy as np
import os
import sys

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
splinedir = os.environ.get("Q_DIR") + '/trunk/Software/Modelling/'
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def qSpline(xvals,yvals,yerrvals,xgrid,allow_out_of_bounds=False):
    
    if not allow_out_of_bounds:
        if min(xgrid) < min(xvals):
            raise ValueError('Desired xgrid out of bounds - min value too low')
        if max(xgrid) > max(xvals):
            raise ValueError('Desired xgrid out of bounds - max value too high')
        
    
    fitpath = storepath + 'regrSpline_fit.txt'
    fiterrpath = storepath + 'regrSpline_fiterr.txt'
    write_ypath = storepath + 'regrSpline_y_in.txt'
    write_xpath = storepath + 'regrSpline_x_in.txt'
    
    ### WRITE OUT THE DATA TO FILE FOR R TO READ IN
    xout=file(write_xpath,'w')
    for val in xgrid:
        outstr = str(val) + '\n'
        xout.write(outstr)
    xout.close()
    
    yout=file(write_ypath,'w')
    ind = 0
    for xval in xvals:
        yval = yvals[ind]
        yeval = yerrvals[ind]
        outstr = str(xval) + ' ' + str(yval) + ' ' + str(yeval) + '\n'
        yout.write(outstr)
        ind += 1
    yout.close()
    
    ### CALL THE R FUNCTION TO CREATE THE SPLINE OUTPUT
    functionstr='R < %scall_regrSpline.R --no-save' % splinedir
    os.system(functionstr)
    
    
    # read in the new lines 
    f=file(fitpath,'r')
    g=file(fiterrpath,'r')
    
    # assuming line format of
    # '"10" 17.719301369 \n'
    newylist=[]
    for line in f.readlines():
        linesplit = line.split()
        try:
            newylist.append(float(linesplit[1]))
        except:
            print "skipping line %s" % line
    f.close()
    
    
    spline_model_errlist=[]
    for line in g.readlines():
        linesplit = line.split()
        try:
            spline_model_errlist.append(float(linesplit[1]))
        except:
            print "skipping line %s" % line
    g.close()
    
    newyarr = np.array(newylist)
    spline_model_errarr = np.array(spline_model_errlist)
    
    return newyarr, spline_model_errarr