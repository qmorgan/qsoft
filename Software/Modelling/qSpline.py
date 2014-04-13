import numpy as np
import os
import sys
from matplotlib import rc
import matplotlib.pyplot as plt

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
splinedir = os.environ.get("Q_DIR") + '/Software/Modelling/'
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def qSpline(xvals,yvals,yerrvals,xgrid,allow_out_of_bounds=False,plot=False):
    
    if not allow_out_of_bounds:
        if min(xgrid) < min(xvals):
            raise ValueError('Desired xgrid out of bounds - min value too low')
        if max(xgrid) > max(xvals):
            raise ValueError('Desired xgrid out of bounds - max value too high')
        
    newyarr, spline_model_errarr = _callRegrSpline(xvals,yvals,yerrvals,xgrid)
    
    if plot:
        qSplinePlot(xvals,yvals,yerrvals)

    return newyarr, spline_model_errarr

def qSplinePlot(xvals,yvals,yerrvals,fig=None,ax_index=None,inverse_y=False,
    inverse_x=False,xlabel='',ylabel='',x_max=None,color='black',fmt='o'):
    '''Use the min and max xvals to generate a plot showing what the continuous 
    plot would be
    
    x_max: only plot up to a particular x value
    
    color: color of the model line and data points
    fmt: shape of the data points
    '''
    plotindices = None
    rc('font', family='Times New Roman')
    if not fig:
        fig=plt.figure()
        ax=fig.add_axes([0.1,0.1,0.8,0.8])
    else:
        if ax_index==None:
            ax_index=0
        ax=fig.get_axes()[ax_index]
    
    print "Rerunning spline fit for plot"
    
    xmin = xvals.min()
    if not x_max: # if not explicitly defined, use all data
        xmax = xvals.max()
    else: #zoom in on a particualr value
        xmax = x_max
        plotindices = np.nonzero(xvals<x_max)[0] # indices; this array should be in order
        # xvals = xvals[plotindices]
        # yvals = yvals[plotindices]
        # yerrvals = yerrvals[plotindices]
        
    model_xgrid = np.linspace(xmin,xmax,num=500)
    modelyarr, modelyerrarr = _callRegrSpline(xvals,yvals,yerrvals,model_xgrid)
    

    # plot model
    ax.plot(model_xgrid,modelyarr,lw=2,color=color)
    ax.plot(model_xgrid,modelyarr-modelyerrarr,lw=1,color=color,alpha=0.5)
    ax.plot(model_xgrid,modelyarr+modelyerrarr,lw=1,color=color,alpha=0.5)
    # plot data
    
    # if we redifined what we want to plot, grab the relevant data values
    if plotindices != None: 
        xvalsplot = xvals[plotindices]
        yvalsplot = yvals[plotindices]
        yerrvalsplot = yerrvals[plotindices]
    else: 
        xvalsplot = xvals
        yvalsplot = yvals
        yerrvalsplot = yerrvals
        
    ax.errorbar(xvalsplot,yvalsplot,yerr=yerrvalsplot,color=color,fmt=fmt)
    
    ax.set_ylabel(ylabel, size=20)
    ax.set_xlabel(xlabel, size=20)
    
    if inverse_y:
        ax.set_ylim(ax.get_ylim()[1],ax.get_ylim()[0])
    if inverse_x:
        ax.set_xlim(ax.get_xlim()[1],ax.get_xlim()[0])
    
    fig.show()
    
    
def _callRegrSpline(xvals,yvals,yerrvals,xgrid):
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