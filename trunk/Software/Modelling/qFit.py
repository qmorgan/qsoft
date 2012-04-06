import scipy
import scipy.integrate
import scipy.optimize
import numpy as np
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
        self.covindex = None # the corresponding index in the covariance matrix
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
    if x is None: x = np.arange(y.shape[0]) 
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
    if not isinstance(paramfinal,np.ndarray):
        paramfinal=np.array([paramfinal])
        
    # Calculate the chi-squared value, and the reduced chi-squared
    # http://mail.scipy.org/pipermail/scipy-user/2005-June/004632.html    
    chi2 = sum(np.power(errfunc(paramfinal), 2)) 
    degrees_of_freedom = y.shape[0] - len(paramfinal)
    chi2r = chi2/degrees_of_freedom
     
    print "chi^2/dof", chi2, '/', degrees_of_freedom
    print "reduced chi^2 = %f  \n" % (chi2r)
        
    retdict = {'parameters':parameters,'covarmatrix':covarmatrix,'chi2':chi2,'dof':degrees_of_freedom}
    
    count=0
    fitstrlist=[]
    for param in retdict['parameters']:
        param.covindex = count # corresponding index in the covmatrix
        uncertainty = np.sqrt(retdict['covarmatrix'].diagonal()[count])
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

def plot_marg_from_fitdict(fitdict,paramnames):
    '''
    Given a fit dictionary from fit(), and a tuple of parameter names (from param.name)
    get the covariance matrix and plot the marginalization
    e.g.
    paramnames = ('Av_1','beta_1')
    '''
    allvalues = np.zeros(len(fitdict['parameters']))
    indices = [-1,-1]
    values = [0,0]
    names = ['a1','a2']
    covmat = fitdict['covarmatrix']
    count = 0
    for param in fitdict['parameters']:
        allvalues[count] = param.value
        count +=1
        if paramnames[0] == param.name:
            indices[0] = param.covindex
            names[0] = param.name
            values[0] = param.value
        if paramnames[1] == param.name:
            indices[1] = param.covindex
            names[1] = param.name
            values[1] = param.value
    
    # print indices
    # return allvalues
    ret = plot_marginalization(covmat=covmat,indices=indices,names=names,values=values)
    return ret
    
def plot_marginalization(covmat=None,indices=None,names=None,values=None):
    
    if covmat == None: # default just for illustration
        covmat=np.matrix([[ 5.29626719,  0.57454987, -0.73125854],
            [ 0.57454987,  1.16079146, -0.28095744],
            [-0.73125854, -0.28095744,  0.23075755]])
    # Suppose we dont care much about some parameters and want to explore 
    # the uncertainties involved in just two. We can marginalize over the 
    # other parameters by first extracting the relevant values from the 
    # covariance matrix above to form a new 2x2 covariance matrix:
    if indices == None: # default for illustration
        ind1 = 1
        ind2 = 2
    else:
        ind1 = indices[0]
        ind2 = indices[1]
    if names == None: # default for names
        names = ['a1','a2']
    
    
    unc_1 = np.sqrt(covmat[ind1,ind1])
    unc_2 = np.sqrt(covmat[ind2,ind2])
    
    # slice the covariance matrix
    cov_i = np.matrix([[covmat[ind1,ind1],covmat[ind1,ind2]],
                        [covmat[ind2,ind1],covmat[ind2,ind2]]])
                        
    # And we invert this to get a new curvature matrix:
    curv_i = cov_i.getI()
    
    #Now, from this curvature matrix we can write (c.f. Eq. 9.3 of Aficionados):
    # \[
    # \Delta \chi^2_{\mathbf{a_i}} = \delta \mathbf{a_i^T} \cdot 
                # [\alpha_\chi]_\mathbf{i} \cdot \delta \mathbf{a_i},
    # \]
    # where $\delta \mathbf{a_i^T} = [\delta a_1 \; \delta a_2]$.
    # 
    # To visualize these $\chi^2$ values, we create a 256 by 256 grid of 
    #  $\delta a_1, \delta a_2$ values and calculate the $\chi^2$ for each.
    
    scale = max((unc_1,unc_2))*3.5
    
    dx=np.linspace(-1*scale,scale,256)
    dy=np.linspace(-1*scale,scale,256)
    delta_chi_sq = np.zeros((256,256))

    # calculate grid of delta_chi_squared:
    x_ind = 0
    for dx_i in dx:
        y_ind = 0
        for dy_j in dy:
            delta_a = np.matrix([[dx_i],[dy_j]])
            delta_chi_sq[x_ind,y_ind] = delta_a.getT()*curv_i*delta_a
            y_ind +=1
        x_ind += 1

    
    # Now with this data grid, we can plot the contours corresponding to 
    # $\Delta \chi^2_{\mathbf{a_i}} = 1.0$ and $\Delta \chi^2_{\mathbf{a_i}} = 2.3.$ 
    # The latter value is where 68.3\% of measured values should lie inside for 
    # two degrees of freedom.  
    # $\Delta \chi^2_{\mathbf{a_i}} = 6.17.$ is where 95.4% of measured values
    # should lie for 2 degrees of freedom (corresponding to 4sigma for 1dof)
    # 9.21 is 99% confidence
    levels = np.array([1.0,2.3,6.17,9.21])
    
    import matplotlib.pyplot as plt
    
    
    CS = plt.contour(dx,dy,delta_chi_sq,levels)
    plt.clabel(CS, inline=1, fontsize=10)
    yname = '$\delta %s$' % names[0]
    xname = '$\delta %s$' % names[1]
    plt.ylabel(yname)
    plt.xlabel(xname)
    

    # marginalize over two parameters and plot the corresponding lines; 
    # these values are equivalent to the uncertainties from the diagonals 
    # of the original covariance matrix above.  It can be seen that they bound 
    # the error ellipse corresponding to$\Delta \chi^2 = 1.0,$ as they should
    # (hence the term 'marginalization' - this contour, projected into the 
    # margins, gives the uncertainty for a single parameter of interest).

    
    plt.axvline(x=unc_2,linestyle='dashed')
    plt.axvline(x=-1*unc_2,linestyle='dashed')
    
    plt.axhline(y=unc_1,linestyle='dashed')
    plt.axhline(y=-1*unc_1,linestyle='dashed')
    
    # 
    if values != None:
        xlim = plt.xlim() + values[1]
        ylim = plt.ylim() + values[0]
        plt.twiny()
        plt.twinx()
        plt.xlim(xlim)
        plt.ylim(ylim)
        yname = '$%s$' % names[0]
        xname = '$%s$' % names[1]
        plt.ylabel(yname)
        plt.xlabel(xname)
    
    return (dx,dy,delta_chi_sq,levels)
    
    
def test_fit():
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
        # x needs to be a single parameter because of the way the errfunc in fit() is defined 
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