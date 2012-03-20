import scipy
import scipy.integrate
import scipy.optimize
import numpy
import pylab
import os, sys
import matplotlib.pyplot as plt

class Param:
    '''Parameter for model fitting.
    
    Grabbed Cdoe from http://www.scipy.org/Cookbook/FittingData
    
    '''
    def __init__(self,value,name='unknown'):
        self.value = value
        self.uncertainty = None
        self.initial_value = value # this will be unchanged after the fit
        self.name=name
    
    def set(self,value):
        self.value = value
    
    def __call__(self):
        # When paramter is called, return itself as the value
        return self.value


def fit(function, parameters, y, yerr, x = None, return_covar=False):
    '''Fit performs a simple least-squares fit on a function.  To use:
    
    Give initial paramaters:
        mu = Param(7)
        sigma = Param(3)
        height = Param(5)
        
    Define your function:
        def f(x): return height() * exp(-((x-mu())/sigma())**2)
    
    Make simulated data:
        xvals=np.arange(100)
        gaussian = lambda x: 3*exp(-(30-x)**2/20.)
        ydata = gaussian(xvals)
        ydata = scipy.randn(100)*.1+ydata #adding noise
        yerr = np.zeros(100)+.1 #array of uncertainties
    
    
    
    Fit the function (provided 'data' is an array with the data to fit):
        fit(f, [mu, sigma, height], ydata, yerr, xvals)
    
    Plot the fitted model over the data if desired
        simxvals = np.arange(10000)/100. # 10000 points from 0-100
        plot(simxvals,f(simxvals))
        scatter(xvals,ydata)
    
    Your input parameters will now have the fitted values assigned to them.
    
    '''
    # This errfunc assumes a uniform error for each datapoint!
    def errfunc(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return (y - function(x))/yerr

    #errfunc = lambda x, y, err: (y-function(x))/err
    
    print "Fitting %i free parameters on %i datapoints...\n" % (len(parameters),len(y))
    
    # If no x axis given, set x array as integral steps starting with 
    # zero for the length of the y array.  For instance, if 
    # y = array([1.32,2.15,3.01,3.92]), then x will be x = array([0,1,2,3])
    if x is None: x = numpy.arange(y.shape[0]) 
    paraminit = [param() for param in parameters]
    #print "Parameter list:", paraminit
    #print "error function:", errfunc
    fitout = scipy.optimize.leastsq(errfunc, paraminit, full_output=1)
    #print fitout
    
    # paramfinal=[param() for param in parameters]
    paramfinal = fitout[0]
    covarmatrix = fitout[1]
    print 'Final Parameters:', paramfinal
    print 'Covariance Matrix', covarmatrix
    print ''
    # If paramfinal is not an array, make it one to avoid scripting errors
    if not isinstance(paramfinal,numpy.ndarray):
        paramfinal=numpy.array([paramfinal])
        
    # Calculate the chi-squared value, and the reduced chi-squared
    # http://mail.scipy.org/pipermail/scipy-user/2005-June/004632.html    
    chi2 = sum(numpy.power(errfunc(paramfinal), 2)) 
    degrees_of_freedom = y.shape[0] - len(paramfinal)
    chi2r = chi2/degrees_of_freedom
     
    print "chi^2/dof", chi2, '/', degrees_of_freedom
    print "reduced chi^2 = %f  \n" % (chi2r)
        
    retdict = {'parameters':parameters,'covarmatrix':covarmatrix,'chi2':chi2,'dof':degrees_of_freedom}
    
    count=0
    fitstrlist=[]
    for param in retdict['parameters']:
        uncertainty = numpy.sqrt(retdict['covarmatrix'].diagonal()[count])
        param.uncertainty = uncertainty # update the uncertainty in the param object
        fitstr = '%s: %.2f +/- %.2f' % (param.name, param.value, uncertainty) 
        print fitstr
        fitstrlist.append(fitstr)
        count += 1
    retdict.update({'strings':fitstrlist})
    
    if return_covar:
        return covarmatrix
    else:
        return retdict # return the dictonary of outputs

def test_fit():
    import numpy as np
    import matplotlib.pyplot as plt
    #Give initial paramaters:
    mu = Param(20,name='mu')
    sigma = Param(4,name='sigma')
    height = Param(5,name='height')
    
    #Define your function:
    def f(x): 
        # here we are using proof of concept of a multi-input function e.g. y=f(x,t)
        # 'x' is a list of tuples, zipped with the zip function (see below)
        # we unzip them and then evaluate. 
        # e.g., x=[(0, 0.0),(1, 0.0),(2, 0.0),(3, 0.0),(4, 0.0),(5, 0.0)]
        xval,zero = zip(*x) # unzip the x feature vector
        return height() * np.exp(-((xval-mu())/sigma())**2) + zero
    
    #Make simulated data:
    xvals=np.arange(100)
    zeros = np.zeros(100)
    zipxvals = zip(xvals,zeros)
    
    gaussian = lambda x: 3*np.exp(-(30-x)**2/20.)
    # true values: mu=30, height=3, sigma=sqrt(20)=4.472
    ydata = gaussian(xvals)
    ydata = scipy.randn(100)*.1+ydata #adding noise
    yerr = np.zeros(100)+.1 #array of uncertainties

    #Fit the function (provided 'data' is an array with the data to fit):
    fit(f, [mu, sigma, height], ydata, yerr, zipxvals)

    #Plot the fitted model over the data if desired
    simxvals = np.arange(10000)/100. # 10000 points from 0-100
    simzeros = np.zeros(len(simxvals))
    zipsimxvals = zip(simxvals,simzeros)
    
    fig2=plt.figure()
    ax=fig2.add_axes([0.1,0.1,0.8,0.8])
    ax.plot(simxvals,f(zipsimxvals))
    
    ax.scatter(xvals,ydata)
    fig2.show()