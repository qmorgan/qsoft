#!/usr/bin/env python
# encoding: utf-8
"""
optiphot.py
Author: Adam Morgan
Last Updated: June 25, 2009
    Cleaned up code, comments.
"""

import scipy
import scipy.integrate
import scipy.optimize
import numpy
import pylab
import os, sys

class LCModel:
    '''Represents a Model Lightcurve'''
    def __init__(self,model):
        self.model = model
        print "Initialized model class"
    
    def integrate(self,tstart,tstop,approx=False):
        '''Integrate the model lightcurve'''
        statement = "Integral of Model from t=%f to t=%f" % (tstart,tstop)
        if approx:
            statement = 'Approximate ' + statement
            texp = tstop - tstart
            tmid = tstart + texp
            # Model evaulated at the midpoint of the exposure times integration time
            integrated_model = self.model(tmid) * texp
        else:
            integrated_model = scipy.integrate.quad(self.model,tstart,tstop)
        print statement
        print integrated_model
        return integrated_model
    
    def simdata(self,tstart=1.21e6,tstop=2.42e6,exptime=1200,cadence=1.21e6, skynoise=0.05*10**-26):  
        # All times in seconds, flux in cgs
        numpoints=(tstop+1-tstart)/cadence
        # Assume same exptime for each exposure for now.
        exptime= numpy.zeros(numpoints) + exptime  
        x_data = numpy.linspace(tstart,tstop,numpoints)
        # We probably don't need this to run in here; will be faster without.  
        # Not needed for header.
        y_data = self.model(x_data)  
        # Make a log-log plot of the simulated data
        pylab.loglog(x_data,y_data,'ro')
        x_data_stop = x_data+exptime
        # Think of some better way to simulate noise.. convert to "counts"
        # skycounts = (skynoise+0.001*scipy.randn(numpoints))*exptime 
        skycounts = (skynoise*exptime)
        header = {'numexp':numpoints,'tstart':x_data,'tstop':x_data_stop, \
                    'exptime':exptime,'skycounts':skycounts}
        return header
    
    def plot(self, modelname="NoModelName"):
        '''Plots the Model Lightcurve; Optimized for built-in RSne'''
        t = numpy.arange(0.1, 10000.0+0.1, 0.1)
        modelname = 'RSNe Models'
        pylab.loglog(t,self.model(t))
        pylab.xlabel('Time (days)')
        pylab.ylabel('Flux Density (mJy)')
        pylab.title(modelname)
        pylab.ylim(0.005,50)
        pylab.xlim(1,10000)
        pylab.show()
    

class GRB(LCModel):
    '''Represents a GRB Model Lightcurve'''
    def __init__(self,R_1,t_1,alpha,name="NoName"):
        self.R_1 = R_1
        self.t_1 = t_1
        self.alpha = alpha
        self.name = name
        self.model = lambda tt: R_1*(tt/t_1)**alpha
        print '(Initialized GRB Model: %s)' % self.name
        

class GRB2(LCModel):
    '''Represents a GRB model via Beuermann (1999)
    
    Two Power-law function
    
    F(t) = (F_1^-n + F_2^-n)^(-1/n)
    where F_i = k_i*t^(-\alpha_i), n > 0
    
    F_1 = F_2 at the transition time t=t_*
    
    5 free parameters (k_1,a_1,k_2,a_2,n)
    '''
    def __init__(self,k_1,a_1,k_2,a_2,n,name="noName"):
        self.k_1 = k_1
        self.k_2 = k_2
        self.a_1 = a_1
        self.a_2 = a_2
        self.n = n
        self.name = name
        F_1 = lambda tt: (self.k_1*tt**(-1*a_1))
        F_2 = lambda tt: (self.k_2*tt**(-1*a_2))
        self.model = lambda tt: ((F_1(tt)**(-1*self.n) + (F_2(tt)**(-1*self.n)))**(-1/self.n))
        print '(Initialized GRB Model: %s)' % self.name
        
class SN(LCModel):
    '''Represents a SN Model Lightcurve'''
    def __init__(self,k_1,k_2,freq,alpha,beta,delta,name="NoName"):
        self.k_1 = k_1
        self.k_2 = k_2
        self.freq = freq
        self.alpha = alpha
        self.beta = beta
        self.delta = delta
        self.name = name
        self.model = lambda tt: (k_1*((freq/5)**alpha)*tt**beta)* \
                                numpy.e**-(k_2*((freq/5)**(-2.1))*tt**delta)
        print '(Initialized SN Model: %s)' % self.name
    

class SNcgs(LCModel):
    '''Represents a SN Model Lightcurve in CGS Units'''
    def __init__(self,k_1,k_2,freq,alpha,beta,delta,name="NoName"):
        k_1 = k_1*10**-26 # convert from mJy to erg/(s*cm^2*Hz)
        self.k_1 = k_1
        self.k_2 = k_2
        self.freq = freq
        self.alpha = alpha
        self.beta = beta
        self.delta = delta
        self.name = name
        # Convert tt in seconds into days
        self.model = lambda tt: (k_1*((freq/5)**alpha)*(tt/86400)**beta)* \
                                numpy.e**-(k_2*((freq/5)**(-2.1))*\
                                (tt/86400)**delta)
        print '(Initialized SN Model: %s)' % self.name
    

class DataFile:
    '''Manipulates the data from various sources'''
    def __init__(self,filename,filetype="NoType"):
        self.filename = filename
        self.filetype = filetype
        if not os.path.exists(self.filename):
            #Change this to call an AssignName function
            print 'The file "%s" does not exist. Exiting..' % self.filename
            sys.exit(1)
        elif self.filetype == "NoType":
            self.AssignType()
            # If not a MIRIAD data set, it cannot be a directory                    
        if os.path.isdir(self.filename) and self.filetype != 'MIRIAD':
            print '"%s" is a directory. Exiting..' %self.filename
            sys.exit(1)
        print '(Initialized DataFile %s of Type %s)' % (self.filename, self.filetype)
    
    def AssignType(self):
        running = True
        typelist = ['FITS','ASCII','MIRIAD']
        while running:
            self.filetype = raw_input("Enter File Type (FITS, ASCII, etc.): ")
            if self.filetype == "quit": 
                print "Exiting..."
                sys.exit(1)
            elif self.filetype != "quit":
                try:
                    typelist.index(self.filetype) # Check if acceptable filetype
                    print "%s has been labeled as type %s" % (self.filename,self.filetype)
                    running = False
                except ValueError:
                    print "Filetype not supported.  Type 'quit' or acceptable filetypes:", typelist
    

class Param:
    '''Parameter for model fitting.
    
    Grabbed Cdoe from http://www.scipy.org/Cookbook/FittingData
    
    '''
    def __init__(self,value):
        self.value = value
    
    def set(self,value):
        self.value = value
    
    def __call__(self):
        # When paramter is called, return itself as the value
        return self.value
    

# ==================
# = Top Level Code =
# ==================
# Here are some Sample Lightcurves:
sn1983n = SN(3300,301,5,-1.08,-1.55,-2.53,"SN1983n")
sn1983n_cgs = SNcgs(3300,301,5,-1.08,-1.55,-2.53,"SN1983n_cgs")
sn1984l = SN(352,301,5,-1.15,-1.56,-2.59,"SN1984l")
sn1990b = SN(177,12400,5,-1.07,-1.24,-2.83,"SN1990b")
grb050525a = GRB(30.938,120.9655,-1.12,"GRB050525a")
testgrb = GRB(3,50,-1.2,"TestGRB")

# **MAKE INTO DEFINITION OF LCMODEL**
def fit(function, parameters, y, yerr, x = None, return_covar=False):
    '''Fit performs a simle least-squares fit on a function.  To use:
    
    Give initial paramaters:
        mu = Param(7)
        sigma = Param(3)
        height = Param(5)
        
    Define your function:
        def f(x): return height() * exp(-((x-mu())/sigma())**2)
        
    Fit the function (provided 'data' is an array with the data to fit):
        fit(f, [mu, sigma, height], data)
    
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
    print 'errors are', numpy.sqrt(covarmatrix[0,0]), numpy.sqrt(covarmatrix[1,1])
    # If paramfinal is not an array, make it one to avoid scripting errors
    if not isinstance(paramfinal,numpy.ndarray):
        paramfinal=numpy.array([paramfinal])
        
    # Calculate the chi-squared value, and the reduced chi-squared
    # http://mail.scipy.org/pipermail/scipy-user/2005-June/004632.html    
    chi2 = sum(numpy.power(errfunc(paramfinal), 2))  
    print "chi^2", chi2
    degrees_of_freedom = y.shape[0] - len(paramfinal)
    chi2r = chi2/degrees_of_freedom
    print "reduced chi^2 = %f on %f degrees of freedom" % (chi2r,degrees_of_freedom)
    if return_covar:
        return covarmatrix
    else:
        return paramfinal

# ==========================================
# = PORTING OLD CODE. MAKE OBJECT ORIENTY. =
# ==========================================

def models2n(sumto,model_input,images_input):
    '''Calculate the summed Signal-to-Noise for a given model lightcurve
    and lightcurve data.  Will sum all data up to "sumto"

    '''
    # If we want to use homegrown simulated data for the input, grab it now.
    if images_input == "simdata":
        header = model_input.simdata()
    else: header = getheaderinfo(images_input)
    
    s2narray = numpy.zeros((header['numexp'],7))
    
    s2narray[:,0] = header['tstart']
    s2narray[:,1] = header['tstop']
    s2narray[:,3] = header['skycounts']
    
    # Calculate c_model_i
    the_model = model_input
    
    for i in range(0,sumto):   # have to use forloop because integrating
        tstart = s2narray[i,0]
        tstop = s2narray[i,1]
    
        c_model_i = LCModel.integrate(the_model,tstart,tstop)
    
        # Now put c_model_i in an array 
        s2narray[i,2] = c_model_i[0]
        
    sum_c_model = sum(s2narray[0:sumto,2])
    
    # Calculate P_i
    
    s2narray[:,5] = s2narray[:,2]/sum_c_model  
    print "sum of P_i should equal 1: ", sum(s2narray[0:sumto,5])
    
    # Calculate weights  
    s2narray[:,6] = s2narray[:,5]**2/(s2narray[:,2]+s2narray[:,3])
    print s2narray[:,6]
    
    # Calculate Noise
    s2narray[:,4] = s2narray[:,2] + s2narray[:,3]  # Noise^2_i
    sum_noise_squared = sum(s2narray[0:sumto,4])
    
    # Calculate unweighted S/N
    sums2n = sum_c_model/numpy.sqrt(sum_noise_squared)
    print sums2n
    
    # Calculate weighted s/n
    numer = s2narray[:,5]*s2narray[:,2]/(s2narray[:,2]+s2narray[:,3])
    denomer = s2narray[:,5]**2/(s2narray[:,2]+s2narray[:,3])
    wtsums2n = sum(numer[0:sumto])/numpy.sqrt(sum(denomer[0:sumto]))
    print wtsums2n
    
    # Cool, that seems to work. Return vector with summed S/N
    sumvec = numpy.zeros((1,2))
    sumvec[0,0] = sums2n
    sumvec[0,1] = wtsums2n
    return sumvec

def sums2nloop(model_input,images_input,iterations):
    '''Run through a loop of calculating the summed S/N for a set of data.
    
    This is useful if you'd like to see how well we'd do as we keep adding
    more and more data together in order to attempt to extract a signal.
    
    '''
    if images_input == "simdata":
        header = model_input.simdata()
    else: header = getheaderinfo(images_input)
    
    # iterations = header['numexp']
    
    snvec = numpy.zeros((iterations,2))
    
    for i in range(0,iterations):
        snvec[i] = models2n(i+1,model_input,images_input)
    
    print snvec
    
    # GRAB EXPOSURE TIMES
    timearr = header['tstart'][0:iterations] # should be error bar from tstart to tstop
    
    # Plot the results
    pylab.semilogx(timearr,snvec[:,1],linestyle='',marker='o')
    pylab.semilogx(timearr,snvec[:,0],linestyle='',marker='x')
    pylab.xlim(50,200000)
    pylab.xlabel('Time')
    pylab.ylabel('Summed S/N (not a lightcurve!)')
    plottitle = 'Model: ' + model_input.name + ', Data: ' + images_input
    pylab.title(plottitle)

def getheaderinfo(images_input):
    '''Grab the relevent information from the images
    
    THIS NEEDS TO BE DEVELOPED MORE, AND NEEDS TO BE PUT IN A CLASS.
    FOR NOW, JUST FAKE IT and PROVIDE WITH FAKE 050525a DATA
    
    '''
    if images_input == '050525a':  
        print "THIS IS A GRB"
        
        numexp = 19
        met = 138672172.8  # MET of burst
        sky_area=113.094  # square arcseconds
        
        #[met_tstart, met_tstop, R_sky/arcsec^2]        
        sarray = numpy.array([[138672251.042, 138672350.824, 0.016],\
                        [138672426.146, 138672435.911, 0.012],\
                        [138672510.570, 138672520.334, 0.015],\
                        [138672594.640, 138672604.417, 0.014],\
                        [138672679.241, 138672689.006, 0.012],\
                        [138672763.610, 138672773.386, 0.013],\
                        [138672848.023, 138672857.799, 0.013],\
                        [138672932.480, 138672942.245, 0.016],\
                        [138673017.003, 138673026.778, 0.019],\
                        [138673101.570, 138673111.335, 0.024],\
                        [138673366.614, 138673465.536, 0.052],\
                        [138690672.489, 138690825.791, 0.046],\
                        [138694046.555, 138694627.140, 0.013],\
                        [138707436.241, 138708186.369, 0.025],\
                        [138719002.042, 138719324.258, 0.020],\
                        [138730582.066, 138731333.639, 0.021],\
                        [138742156.330, 138742907.032, 0.022],\
                        [138754206.432, 138754479.938, 0.041],\
                        [138758451.217, 138759350.996, 0.013]])
        
        tstart = sarray[:,0] - met
        tstop = sarray[:,1] - met
        exptime = tstop - tstart
        skycounts = (sarray[:,2]*sky_area)*exptime   #c_sky_i
    elif images_input == 'testsn':
        print "TEST SN"
        tstart = tstart - 70
        tstop = tstop - 70
    else: 
        print "Cannot parse file.  Exiting."
        sys.exit(1)
        
    header = {'numexp':numexp,'tstart':tstart,'tstop':tstop, \
            'exptime':exptime,'skycounts':skycounts}
    return header

# The end.
