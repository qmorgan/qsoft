import matplotlib.pyplot as plt
import numpy as np
import copy
from MiscBin.q import filtcheck
from MiscBin.q import mag2flux
from MiscBin import q
from MiscBin import qObs
from Phot import qFit
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
from Phot import PhotParse

class Extinction:
    '''Represents an extinction law'''
    def __init__(self):
        pass

        
    def UnreddenFlux(self,flux):
        '''        
        flux - calibrated flux vector, same number of elements as self.wave
        '''
        if not hasattr(self,'curve'):
            raise Exception("Need to run EvalCurve first")
        if np.isscalar(flux):
            flux = np.array([flux]) # convert the scalar to a one element array
        self.flux = flux
        
        if not len(self.flux) == len(self.wave):
            raise ValueError("Flux and wavelength arrays are not equal")
        
        self.funred = self.flux * 10.**(0.4*self.curve)
    
    
    def plotModel(self,fig=None,color='blue',show=True,loglog=False,xlim=None):    
        if not hasattr(self,'curve'):
            raise Exception("Need to run EvalCurve first")
        if not hasattr(self,'funred'):
            raise Exception("Need to apply the reddening with UnreddenFlux")
        if not len(self.flux) == len(self.wave):
            raise ValueError("Flux and wavelength arrays are not equal")
        
        if not fig:
            fig=plt.figure()
            ax=fig.add_axes([0.1,0.1,0.8,0.8])
        else:
            ax=fig.get_axes()[0]
        ax.plot(self.wave,self.funred,linewidth=2,color=color)
        ax.plot(self.wave,self.flux,linewidth=1,color=color,linestyle='dashed')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux Density')
        if loglog:
            ax.loglog()
        if xlim:
            ax.set_xlim(xlim)
        if show:
            plt.show()
        else:
            return(fig)
    
        
class FM(Extinction):
    def __init__(self,R_V,c1,c2,c3,c4,gamma,x0):
        '''
        Using the Fitzpatrick (1999) parameterization
        
        The R-dependent Galactic extinction curve is that of Fitzpatrick & Massa
        (Fitzpatrick, 1999, PASP, 111, 63; astro-ph/9809387 ).    
        Parameterization is valid from the IR to the far-UV (3.5 microns to 0.1 
        microns).    UV extinction curve is extrapolated down to 912 Angstroms.
        
        R_V - scalar specifying the ratio of total to selective extinction
                 R(V) = A(V) / E(B - V).    If not specified, then R = 3.1
                 Extreme values of R(V) range from 2.3 to 5.3
        
        The following five input keyword parameters allow the user to customize
        the adopted extinction curve.    For example, see Clayton et al. (2003,
        ApJ, 588, 871) for examples of these parameters in different interstellar
        environments.
        
        x0 - Centroid of 2200 A bump in microns 
        gamma - Width of 2200 A bump in microns 
        c3 - Strength of the 2200 A bump
        c4 - FUV curvature 
        c2 - Slope of the linear UV extinction component 
        c1 - Intercept of the linear UV extinction component              
        
        NOTES:
         (1) The following comparisons between the FM curve and that of Cardelli,
             Clayton, & Mathis (1989), (see ccm_unred.pro):
             (a) - In the UV, the FM and CCM curves are similar for R < 4.0, but
                   diverge for larger R
             (b) - In the optical region, the FM more closely matches the
                   monochromatic extinction, especially near the R band.
         (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
                 can be represented with the FM curve, if the proper value of 
                 R(V) is supplied.     
                                                            
        Based of the IDL astronomer user library FM_UNRED.pro by Wayne Landsman
        '''
        self.R_V=R_V
        self.c1=c1
        self.c2=c2
        self.c3=c3
        self.c4=c4
        self.gamma=gamma
        self.x0=x0
    
    def EvalCurve(self,wave,ebv):
        '''Evaluate the extinction curve for a given array of wavelengths.
        
        wave - wavelength vector (Angstroms)

        ebv  - color excess E(B-V), scalar.  If a negative EBV is supplied,
                 then fluxes will be reddened rather than dereddened.
        
        note A_V = R_V*E(B-V)
        '''
        self.ebv=ebv
        if np.isscalar(wave):
            wave = np.array([wave]) # convert the scalar to a one element array
        self.wave=wave

        # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function
        # and R-dependent coefficients
        x = 10000./ self.wave               # ; Convert to inverse microns 
        curve = x*0.
        
        xcutuv = 10000.0/2700.0 #UV Cutoff
        xspluv = 10000.0/np.array([2700.0,2600.0])
        
        # Find the indices of where the input wave vector is in the UV and 
        # where it is in the optical/ir
        
        iuv = np.nonzero(x>=xcutuv)[0] # indices in the UV
        iopir = np.nonzero(x<xcutuv)[0] # indicies in the optical/ir
        N_UV = len(iuv) # number of UV elements in the wave vector
        N_opir = len(iopir) # number of optical/ir elements in the wave vector
        
        if N_UV > 0:
            xuv = np.append(xspluv,x[iuv])
        else:
            xuv = xspluv
        
        xuvLT5p9ind = np.nonzero(xuv<5.9)[0] # indices where xuv < 5.9 
        xuvGT5p9 = copy.copy(xuv)
        xuvGT5p9[xuvLT5p9ind] = 5.9 # replacing all elements less than 5.9 with 5.9
        
        yuv = self.c1 + self.c2*xuv
        yuv = yuv + self.c3*xuv**2/((xuv**2 - self.x0**2)**2 + (xuv*self.gamma)**2)
        yuv = yuv + self.c4*(0.5392*((xuvGT5p9)-5.9)**2 + 0.05644*((xuvGT5p9)-5.9)**3)
        yuv = yuv + self.R_V
        
        yspluv = yuv[0:2]
        
        if N_UV > 0:
            curve[iuv] = yuv[2:] #remove UV spline points from final curve if there are UV points
        
        # setting up the x points of the oir spline
        xsplopir=np.append(np.array([0]),10000./np.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]))
        
        ysplir = np.array([0.0,0.26469,0.82925])*self.R_V/3.1
        ysplop0 = -4.22809e-01 + self.R_V*1.00270 + self.R_V**2*2.13572e-04
        ysplop1 = -5.13540e-02 + self.R_V*1.00216 + self.R_V**2*-7.35778e-05
        ysplop2 = 7.00127e-01 + self.R_V*1.00184 + self.R_V**2*-3.32598e-05
        ysplop3 = 1.19456 + self.R_V*1.01707 + self.R_V**2*-5.46959e-03 + self.R_V**3*7.97809e-04 + self.R_V**4*-4.45636e-05
        ysplop = np.array([ysplop0,ysplop1,ysplop2,ysplop3])
        ysplopir = np.append(ysplir,ysplop)
        
        xspluvopir = np.append(xsplopir,xspluv)
        yspluvopir = np.append(ysplopir,yspluv)

        # the spline differes at about the 2nd-3rd decimal place compared to
        # the IDL implementation
        from scipy.interpolate import InterpolatedUnivariateSpline
        uvopirspline = InterpolatedUnivariateSpline(xspluvopir,yspluvopir,k=3)
        curve[iopir] = uvopirspline(x[iopir])
        
        curve=ebv*curve
        self.curve = curve
        
        
        
def avgsmc(R_V = 2.74):
    '''
    These models should just take Av, Beta, flux vector, and wavelength vector, 
    and return a new flux vector..
    
    In the average SMC fit, the only free parameters are Av and Beta'''    
    c1 = -4.959
    c2 =  2.264
    c3 = 0.389
    c4 = 0.461
    gamma=1.05
    x0=4.626
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel
    
def lmc2(R_V = 3.1):   
    c1 = -2.16
    c2 = 1.31
    c3 = 1.92
    c4 = 0.42
    gamma = 1.05
    x0 = 4.626
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel

def avglmc(R_V = 3.2):    
    c1 = -1.28
    c2 = 1.11
    c3 = 2.73
    c4 = 0.64
    gamma = 0.91
    x0 = 4.596  
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel
    
def avgmw(R_V = 3.1,c3 = 3.23,c4 = 0.41):
    '''Take a look at Reichart 2001 for an empirical relation between
    c1, c2, and R_V to remove the degeneracy between them, recommended if you
    dont have amazing enough data to try to constrain all 6 parameters, but 
    amazing enough to try to constrain c3 and c4.
    
    Gamma and x0 can generally be fixed if you do not have enough free paraamters 
    to break the degeneracy.
    '''
    c2 = -0.824 + 4.717/R_V
    c1 = 2.030 - 3.007*c2
    gamma = 0.99
    x0 = 4.596
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel

def powerlawExt(wave,flux,Av=0.0,beta=0.0):
    '''
    wave = wavelength vector in angstroms 
    Assumes a power law of the type F = F_0*(nu/nu_0)^-beta
    
    FOR NOW JUST USE SMC - Update this if you ever want to use this function
    for something other than TestPowerlawExt
    '''
    c = 2.998E10 #cm/s
    nu = c/(wave*1E-8) #hertz, for wave in angstroms
    nu_0 = min(nu) #valid???
    flux_normal = flux*(nu/nu_0)**beta
    
    Rv=2.74
    extmodel=avgsmc(R_V=Rv)
    extmodel.EvalCurve(wave,Av/Rv)
    extmodel.UnreddenFlux(flux_normal)
    return extmodel


def timeDepAvBeta(wave_time_vec,norm_vec,Av_0=0.0,Av_1=0.0,Av_2=0.0,beta_0=0.0,beta_1=0.0,
                    beta_2=0.0,Rv=2.74,c1=-4.959,
                    c2=2.264,c3=0.389,c4=0.461,gamma=1.05,x0=4.626):
    '''Allowing for Av and Beta to change with time
    time_wave_dict: wavelength vector in angstroms [[1,2,3,4,5],[1,2,4,5],[1,3,4,5]]
    time: time vector in seconds [50,250,500]
    norm_vec: vector of normalizations (replacing const) - one for each time
    Av_0: limit at large time
    Av_1: normalization const
    Av_2: decay index
    beta_0: limit at large time
    beta_1: normalization const
    beta_2: decay index
    
    Av(t) = Av_0 + Av_1*t^{Av_2}
    beta(t) = beta_0 + beta_1*t^{beta_2}
    
    f(t) = f_corr(t)*(nu/nu_0)^-beta(t)
    where f_corr(t) is the time-dependent extinction corrected normalized flux
    '''
    assert len(wave) == len(time)
    assert len(norm_vec) == len(time)
    
    # build up x-vector of tuples and y_vector and y_err vector
    
    
    # assume extmodel doesn't change with time
    extmodel=FM(Rv,c1,c2,c3,c4,gamma,x0)
    # for each time
    count = 0 # need?
    flux_out=[]
    for time in time_vec:
        ind = time_vec.index(time) # index
        const = norm_vec[ind]
        wave = np.array(wave_vec[ind]) # vector of wavelengths; should be multi-length
         
        beta = beta_0 + beta_1*(time**(beta_2))
        Av = Av_0 + Av_1*(time**(Av_2))
        
        fluxarr = fluxarr=np.ones(len(wave)) #just an array of ones
        
        extmodel.EvalCurve(wave,Av/Rv)
        extmodel.UnreddenFlux(fluxarr)
        
        extmodel.UnreddenFlux(fluxarr)
        flux_normal = const*extmodel.funred*(10000/wave)**beta
        flux_out.append(flux_normal) # append the resultant array to the list
        
    
def powerlawExtRetFlux(wave,Av=0.0,beta=0.0,const=1.0E3,Rv=2.74,c1=-4.959,
                        c2=2.264,c3=0.389,c4=0.461,gamma=1.05,x0=4.626):
    
    '''
    spectral power law + extinction - return flux
    
    Given a vector of wavelengths, spectral power law slope (beta), 
    extinction in V (Av), FM params (Rv,c1,c2,c3,c4,gamma,x0), and normalization
    constant, return the flux vector. 
     
    wave = wavelength vector in angstroms 
    Assumes a power law of the type F = F_0*(nu/nu_0)^-beta
    
    '''
    fluxarr=np.ones(len(wave)) #just an array of ones
    # c = 2.998E10 #cm/s
    # nu = c/(wave*1E-8) #hertz, for wave in angstroms
    # # nu_0 = const #valid???    
    
    extmodel=FM(Rv,c1,c2,c3,c4,gamma,x0)
    extmodel.EvalCurve(wave,Av/Rv)
    extmodel.UnreddenFlux(fluxarr)
    flux_normal = const*extmodel.funred*(10000/wave)**beta
    # here we're doing 10000/wave to get the constant into a more reasonable
    # range of values for the fit. if we did actual frequency, the constant
    # would be too large and the fit would be too difficult.
    return flux_normal

def SEDFitTest(initial_param='smc'):
    '''Proof of concept of SED fitting, using averaged fluxes from dans 
    lcurve fitting.  
    NOTE THE ERRORS ARE MADE UP HERE. Do not trust wholly.
    '''
    filtlist=[qObs.B,qObs.r,qObs.Rc,qObs.i,qObs.Ic,qObs.z,qObs.J,qObs.H,qObs.Ks]
    z=1.728
    galebv=0.108 
    fluxarr=np.array([4.52318,17.7811,19.9167,38.1636,48.2493,78.5432,145.145,288.604,499.728])
    fluxerrarr=np.array([0.2,0.85,1.0,2.0,2.5,3.9,6.0,11.0,20.0])
    SEDFit(filtlist,fluxarr,fluxerrarr,initial_param=initial_param,z=z,galebv=galebv)

def SEDFitTest2(initial_param='smc'):
    '''another SED fit test, this time starting with magnitudes'''
    z=1.728
    galebv=0.108
    
    # using SMARTS first epoch magnitudes
    filtlist=[qObs.B,qObs.V,qObs.Rc,qObs.Ic,qObs.J,qObs.H,qObs.Ks]
    maglist=[20.07,18.91,18.02,16.96,15.14,14.02,12.94]
    magerrlist=[0.06,0.05,0.03,0.03,0.05,0.06,0.10]
    # Convert the fluxes to magnitudes
    fluxarr, fluxerrarr = maglist2fluxarr(maglist,magerrlist,filtlist)
    SEDFit(filtlist,fluxarr,fluxerrarr,initial_param=initial_param,z=z,galebv=galebv)
    

def _align_SED_times(objblock,sedtimelist,time_thresh=10):
    '''
    Given an object block and a list of desired times, return a dictionary 
    of the following format: 
    {'t = 1216.8 s': {'filtlist': [<MiscBin.qObs.filt instance at 0x6931800>,
                                   <MiscBin.qObs.filt instance at 0x69318a0>,
                                   <MiscBin.qObs.filt instance at 0x6931850>,
                                   <MiscBin.qObs.filt instance at 0x6931878>,
                                   <MiscBin.qObs.filt instance at 0x6931a30>,
                                   <MiscBin.qObs.filt instance at 0x6931828>],
                      'magerrlist': [0.025,
                                     0.0413638318369,
                                     0.0389245320622,
                                     0.0364034889082,
                                     0.117,
                                     0.028],
                      'maglist': [17.337,
                                  12.2781,
                                  14.4305,
                                  13.25495,
                                  19.092,
                                  16.118],
                      'sedtime': 1216.82470999968},
     't = 178.2 s': {'filtlist': [<MiscBin.qObs.filt instance at 0x6931800>,
                                  <MiscBin.qObs.filt instance at 0x69318a0>,
                                  <MiscBin.qObs.filt instance at 0x6931850>,
                                  <MiscBin.qObs.filt instance at 0x6931878>,
                                  <MiscBin.qObs.filt instance at 0x6931828>],
                     'magerrlist': [0.112,
                                    0.0597335757223,
                                    0.0869512692113,
                                    0.0659107210151,
                                    0.063],
                     'maglist': [17.555, 11.8928, 14.7104, 13.1396, 16.233],
                     'sedtime': 178.153240999968}}
    '''
    
    aligndict={}
    for sedtime in sedtimelist:
        timestr = 't = %.1f s' % (sedtime)
        maglist = []
        magerrlist = []
        filtlist = []
        # loop through the obsblocks to collect information 
        for obsblock in objblock.obsdict.itervalues():
            if obsblock.filt != None: # IF WE HAVE A FILTER
                tmidarr=np.array(obsblock.tmidlist)
                # determine the number of observations, if any, within the threshold
                matched_indices = np.nonzero(abs(tmidarr-sedtime)<time_thresh)[0]
                # find the length. if its 1, we found one, if its >1, something is wrong
                if len(matched_indices) == 0:
                    continue
                elif len(matched_indices) != 1:
                    raise Exception("Too many matched indices! Threshold too high?")
                elif len(matched_indices) == 1: # we found one. Build up the lists.
                    index = matched_indices[0]
                    maglist.append(obsblock.maglist[index])
                    magerrlist.append(obsblock.magerrlist[index])
                    filtlist.append(obsblock.filt)
        aligndict.update({timestr:{'sedtime':sedtime,'maglist':maglist,'magerrlist':magerrlist,'filtlist':filtlist}})
    return aligndict
    
def SEDvsTime(initial_param='smc',z=1.728,galebv=0.108,
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedData.dat'):
    '''A function which will take a phot objBlock object from the revamping of
    PhotParse, loop through each time in a given time list and search through 
    each set of observations to find those times that overlap within a given 
    threshhold of that time, and build up an SED for each.'''
    
    
    objblock=PhotParse.PhotParse(directory)
    utburststr = objblock.utburst
    
    time_thresh = 10 # Number of seconds we can be off in time from the reference 

    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    
    aligndict = _align_SED_times(objblock,sedtimelist,time_thresh=time_thresh)
    
    # loop through each time, building up an SED for each
    for timestr, val in aligndict.iteritems():
        maglist=val['maglist'] 
        magerrlist=val['magerrlist']
        filtlist=val['filtlist']
        # convert collected magnitudes to fluxes:
        fluxarr, fluxerrarr = maglist2fluxarr(maglist,magerrlist,filtlist)
        print '**'
        print [filt.name for filt in filtlist]
        print maglist
        print magerrlist
        print fluxarr
        print fluxerrarr
        print '**'
        # perform SED fit
        SEDFit(filtlist,fluxarr,fluxerrarr,initial_param=initial_param,z=z,galebv=galebv,
                fit_list=['const','Av','beta'], timestr=timestr)
        
def maglist2fluxarr(maglist,magerrlist,filtlist):
    '''Given a list of magnitudes, their uncertainties, and filt objects,
    return arrays of fluxes and flux errors.'''
    for filt in filtlist:
        filtcheck(filt)
    assert len(maglist) == len(magerrlist)
    assert len(maglist) == len(filtlist)
    # build up flux array and flux err array
    fluxarr=[]
    fluxerrarr=[]
    count=0
    for mag in maglist:
        magerr=magerrlist[count]
        filt=filtlist[count]
        fluxtup=mag2flux(mag_1=mag,mag_1err=magerr,filt=filt)
        flux=fluxtup[0]
        fluxerr=(fluxtup[1]+fluxtup[2])/2.0 # just take the avg error
        fluxarr.append(flux)
        fluxerrarr.append(fluxerr)
        count+=1  # I forgot this initially, wow!
    fluxarr=np.array(fluxarr)
    fluxerrarr=np.array(fluxerrarr)
    return fluxarr, fluxerrarr
    
def SEDFit(filtlist,fluxarr,fluxerrarr,initial_param='smc',z=0.0,galebv=0.0,
            fit_list = ['const','Av','beta'],timestr='', plot=True):
    '''Fit an SED to a list of fluxes with a FM paramaterization.
    
    Required Inputs:
    * filtlist: list of filt objects (see qObs.filt) 
    * fluxarr: list of fluxes
    * fluxerrarr: list of uncertainties on the fluxes
    
    Optional Inputs:
    * initial_param: set of FM parameters to start with 
        * allowed values: 'smc', 'lmc', 'lmc2', 'mw' (acceptable_initial_param_list)
    
    * fit_list: 
        const - normalization
        beta - power law spectral slope
        Av - Extinction in V
        Rv - scalar specifying the ratio of total to selective extinction
                 R(V) = A(V) / E(B - V).  
        x0 - Centroid of 2200 A bump in microns 
        gamma - Width of 2200 A bump in microns 
        c3 - Strength of the 2200 A bump
        c4 - FUV curvature 
        c2 - Slope of the linear UV extinction component 
        c1 - Intercept of the linear UV extinction component
    ''' 
    
    # error checking
    for filt in filtlist:
        filtcheck(filt)
    assert len(fluxarr) == len(fluxerrarr)
    assert len(fluxarr) == len(filtlist)
    
    # Check if initial parameter list is there
    acceptable_fit_param_list=['const','beta','Av','Rv','x0','gamma','c1','c2','c3','c4']
    for fitparam in fit_list:
        if not fitparam in acceptable_fit_param_list:
            print 'Requested fit parameter of %s is not valid. Please choose from:' % fitparam
            return
    acceptable_initial_param_list=['smc','lmc','lmc2','mw']
    if not initial_param in acceptable_initial_param_list:
        print 'Initial parameter set of %s is not valid. Please choose from:' % (initial_param)
        print acceptable_initial_param_list
        return
    
    rc('font', family='Times New Roman') 
    
    # extract values from the filter list
    wavearr=[]
    wavenamearr=[]
    for filt in filtlist:
        wavearr.append(filt.wave_A)
        wavenamearr.append(filt.name)
    wavearr=np.array(wavearr)
    
    # correct for galactic extinction
    if galebv == 0.0:
        print "\nWARNING: Not correcting for galactic extinction\n"
    
    # Use the FM implementation of the average milkyway extinction to unredden
    # the flux due to the galactic sightlines.
    galAv=galebv*3.1 #Rv=3.1 for mw
    mw=avgmw()
    mw.EvalCurve(wavearr,galebv)
    mw.UnreddenFlux(fluxarr) #need to correct the uncertainties as well
    galcorrectedfluxarr=mw.funred
    mw.UnreddenFlux(fluxerrarr) #since uncertainties also scale with extinction, can do direct correction 
    galcorrectedfluxerrarr=mw.funred
    
    #correct for redshift
    waverestarr=wavearr/(1+z)
    
    # determine initial values
    Av_init = -0.62
    beta_init = -1.45
    const_init = 1000
       
    if initial_param == 'smc':
        Rv_init = 2.74
        c1_init = -4.959
        c2_init = 2.264        
        c3_init = 0.389
        c4_init = 0.461
        gamma_init = 1.05
        x0_init = 4.626

    elif initial_param == 'lmc':
        Rv_init = 3.2
        c1_init = -1.28
        c2_init = 1.11
        c3_init = 2.73
        c4_init = 0.64
        gamma_init = 0.91
        x0_init = 4.596
    
    elif initial_param == 'lmc2':
        Rv_init = 3.1
        c1_init = -2.16
        c2_init = 1.31
        c3_init = 1.92
        c4_init = 0.42
        gamma_init = 1.05
        x0_init = 4.626
    
    elif initial_param == 'mw':
        Rv_init = 3.1
        c3_init = 3.23
        c4_init = 0.41
        c2_init = -0.824 + 4.717/Rv_init
        c1_init = 2.030 - 3.007*c2_init
        gamma_init = 0.99
        x0_init = 4.596

    
    #set parameters
    Av=qFit.Param(Av_init,name='Av')
    beta=qFit.Param(beta_init,name='beta')
    const=qFit.Param(const_init,name='const')
    Rv=qFit.Param(Rv_init,name='Rv')
    c1 = qFit.Param(c1_init,name='c1')
    c2 =  qFit.Param(c2_init,name='c2')
    c3 = qFit.Param(c3_init,name='c3')
    c4 = qFit.Param(c4_init,name='c4')
    gamma= qFit.Param(gamma_init,name='gamma')
    x0=qFit.Param(x0_init,name='x0')
        
    # build up the param list of things we want to fit
    fullparamlist = [Av,beta,const,Rv,c1,c2,c3,c4,gamma,x0] # list of ALL possible parameters, fit or otherwise
    fitparamlist = [] 
    fixparamlist = []
    fix_list = []
    for param in fullparamlist:
        if param.name in fit_list:
            fitparamlist.append(param) # if its one of the params we want to fit, append it
        else:
            fix_list.append(param.name)
            fixparamlist.append(param)
    # error check        
    assert len(fit_list) == len(fitparamlist)
        
    
    def f(x): return powerlawExtRetFlux(x,Av=Av(),beta=beta(),Rv=Rv(),const=const(),
                c1=c1(),c2=c2(),c3=c3(),c4=c4(),gamma=gamma(),x0=x0())
    
    fitdict = qFit.fit(f,fitparamlist,galcorrectedfluxarr,galcorrectedfluxerrarr,waverestarr)
    
    
    if plot:
        #get model array
        w=1000. + np.arange(500)*100
        f = w*0. + 1
        c=powerlawExtRetFlux(w,Av=Av.value,beta=beta.value,Rv=Rv.value,const=const.value,
            c1=c1.value,c2=c2.value,c3=c3.value,c4=c4.value,gamma=gamma.value,x0=x0.value)
        # fig = c.plotModel(show=False,color='green')

        
        fig2=plt.figure()
        ax=fig2.add_axes([0.1,0.1,0.8,0.8])
    
        string = 'chi2 / dof = %.2f / %i' % (fitdict['chi2'],fitdict['dof'])
        fig2.text(0.2,0.2,string)
        textoffset=0.24    
        for string in fitdict['strings']:
            if not string.find('const') != -1:
                fig2.text(0.2,textoffset,string)
                textoffset+=0.04
    
    
        if timestr:
            fig2.text(0.45,0.8,timestr)
    
        string = 'Fixed params (%s)' % (initial_param)
        fig2.text(0.7,0.8,string)
        textoffset=0.76        
        for fixedparam in fixparamlist:
            string = fixedparam.name + ': ' + str(fixedparam.value)
            fig2.text(0.7,textoffset,string)
            textoffset-=0.04
    
        #underplot the model
        ax.plot(w,c) 

    
        # i cant seem to get a scatter plot to appear above a line plot. HMM. 
        # Forget it. No need to have colors on the points.
        # ax.scatter(waverestarr,galcorrectedfluxarr,c=wavearr,cmap='jet',s=30,edgecolors='none')
    
        # plot data
        ax.errorbar(waverestarr,galcorrectedfluxarr,yerr=galcorrectedfluxerrarr,fmt='.')
        # annotate data with the filter names
        for name in wavenamearr:
            ind = wavenamearr.index(name)
            xy=(waverestarr[ind],galcorrectedfluxarr[ind])
            xytext=(waverestarr[ind]*0.95,galcorrectedfluxarr[ind]*1.05)
            ax.annotate(name,xy=xy,xytext=xytext,fontsize='small')

        
        # get the bounds for the other axes, which are to be AB mags and obs wavelength
        # this will ensure the shared axes properly line up
        ymax=galcorrectedfluxarr.max()+1.0*galcorrectedfluxarr.max()
        #ymin = galcorrectedfluxarr.min()-0.5*galcorrectedfluxarr.min()
        ymin = ymax*(1-0.997)
        ylimflux=(ymin,ymax)
        xlimrest=(10000,1000)
        ylimmag0=q.flux2abmag(ylimflux[0])
        ylimmag1=q.flux2abmag(ylimflux[1])
        ylimmag=(ylimmag0,ylimmag1)
        xlimobs0=xlimrest[0]*(1+z)
        xlimobs1=xlimrest[1]*(1+z)
        xlimobs=(xlimobs0,xlimobs1)   
    
        # Initialize the other shared axes and make some of them log axes
        ax.loglog()        
        ax2=ax.twinx()
        ax3=ax.twiny()
        ax3.semilogx()
    
        # Set the limits of the axes based on our definitions above.
        ax.set_xlim(xlimrest)
        ax.set_ylim(ylimflux)    
        ax3.set_xlim(xlimobs)
        ax2.set_ylim(ylimmag)
    
        # Label the 4 axes
        ax.set_ylabel(r'$F_\nu$ (uJy)')
        ax.set_xlabel(r'$\lambda_{\mathrm{eff,rest}}$ ($\AA$)')
        ax2.set_ylabel('AB Mag')
        ax3.set_xlabel(r'$\lambda_{\mathrm{eff}}$ ($\AA$)')

        # Explicitly define which ticks to label on the x axes
        ax.set_xticks([10000,6000,4000,2000,1000])
        ax3.set_xticks([20000,10000,6000,4000,3000])

        # use the function formatter for the ticks to show them as full integers
        # rather then exponential format
        formatter = FuncFormatter(log_10_product)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax3.xaxis.set_major_formatter(formatter)

        fig2.show()

def log_10_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    return '%1i' % (x)

def TestPowerLawExt():
    w=1250. + np.arange(100)*100.
    f = w*0. + 1
    a=powerlawExt(w,f,Av=-0.5,beta=-0.5)
    fig = a.plotModel(show=False,color='violet',loglog=True,xlim=(12000,1000))
    
    b=powerlawExt(w,f,Av=-0.5,beta=-0.7)
    fig = b.plotModel(fig=fig,show=False,color='blue')
    
    c=powerlawExt(w,f,Av=-0.5,beta=-0.9)
    fig = c.plotModel(fig=fig,show=False,color='green')
    
    d=powerlawExt(w,f,Av=-0.5,beta=-1.1)
    fig = d.plotModel(fig=fig,show=False,color='orange')
    
    e=powerlawExt(w,f,Av=-0.5,beta=-1.3)
    fig = e.plotModel(fig=fig,show=False,color='red')
    
    fig.text(0.2,0.54,'SMC law',color='black')
    fig.text(0.2,0.5,'Rv=2.74',color='black')
    fig.text(0.2,0.46,'Av=-0.5',color='black')
    fig.text(0.2,0.42,'beta=-0.5',color='violet')
    fig.text(0.2,0.38,'beta=-0.7',color='blue')
    fig.text(0.2,0.34,'beta=-0.9',color='green')
    fig.text(0.2,0.30,'beta=-1.1',color='orange')
    fig.text(0.2,0.26,'beta=-1.3',color='red')
    
    fig.show()
    
    a=powerlawExt(w,f,Av=-0.1,beta=-0.5)
    fig = a.plotModel(show=False,color='violet',loglog=True,xlim=(12000,1000))
    
    b=powerlawExt(w,f,Av=-0.1,beta=-0.7)
    fig = b.plotModel(fig=fig,show=False,color='blue')
    
    c=powerlawExt(w,f,Av=-0.1,beta=-0.9)
    fig = c.plotModel(fig=fig,show=False,color='green')
    
    d=powerlawExt(w,f,Av=-0.1,beta=-1.1)
    fig = d.plotModel(fig=fig,show=False,color='orange')
    
    e=powerlawExt(w,f,Av=-0.1,beta=-1.3)
    fig = e.plotModel(fig=fig,show=False,color='red')
    
    fig.text(0.2,0.54,'SMC law',color='black')
    fig.text(0.2,0.5,'Rv=2.74',color='black')
    fig.text(0.2,0.46,'Av=-0.1',color='black')
    fig.text(0.2,0.42,'beta=-0.5',color='violet')
    fig.text(0.2,0.38,'beta=-0.7',color='blue')
    fig.text(0.2,0.34,'beta=-0.9',color='green')
    fig.text(0.2,0.30,'beta=-1.1',color='orange')
    fig.text(0.2,0.26,'beta=-1.3',color='red')
    
    fig.show()
    
def Test120119A():
    z=1.72
    wv_angst=np.array([  6470. ,   7865.,  12500. ,  16500. ,  21500. ])
    wv_angst_rest = wv_angst/(1+z)
    
    w=1250. + np.arange(100)*100.
    f = w*0. + 1
    a=lmc2(R_V=3.1)
    a.EvalCurve(w,-0.1)
    a.UnreddenFlux(f)
    fig = a.plotModel(show=False,color='green',loglog=True,xlim=(12000,1000))
    b=avgsmc(R_V=3.1)
    b.EvalCurve(w,-0.1)
    b.UnreddenFlux(f)
    fig=b.plotModel(fig=fig,show=False,color='blue')
    c=avgmw(R_V=3.1)
    c.EvalCurve(w,-0.1)
    c.UnreddenFlux(f)
    fig=c.plotModel(fig=fig,show=False,color='red')
    
    fig.text(0.2,0.5,'E(B-V)=-0.1',color='black')
    fig.text(0.2,0.46,'Rv=3.1',color='black')
    fig.text(0.2,0.42,'SMC',color='blue')
    fig.text(0.2,0.38,'LMC',color='green')
    fig.text(0.2,0.34,'MW',color='red')
    fig.show()
    # Note in the figure above, all extinction models converge in the optical
    # when R_V is set to a given number.
    
    # Now see how a change in Rv manifests with the other parameters set as SMC
    av=-0.31
    rv=4.3
    q=avgsmc(R_V=rv)
    q.EvalCurve(w,av/rv)
    q.UnreddenFlux(f)
        
    rv=3.9
    a=avgsmc(R_V=rv)
    a.EvalCurve(w,av/rv)
    a.UnreddenFlux(f)
    
    rv=3.5
    b=avgsmc(R_V=rv)
    b.EvalCurve(w,av/rv)
    b.UnreddenFlux(f)
    
    rv=3.1
    c=avgsmc(R_V=rv)
    c.EvalCurve(w,av/rv)
    c.UnreddenFlux(f)
    
    rv=2.7
    d=avgsmc(R_V=rv)
    d.EvalCurve(w,av/rv)
    d.UnreddenFlux(f)
    
    rv=2.4
    e=avgsmc(R_V=rv)
    e.EvalCurve(w,av/rv)
    e.UnreddenFlux(f)
    
    fig2=q.plotModel(show=False,color='purple',loglog=True,xlim=(12000,1000))
    fig2=a.plotModel(fig=fig2,show=False,color='blue')
    fig2=b.plotModel(fig=fig2,show=False,color='green')
    fig2=c.plotModel(fig=fig2,show=False,color='yellow')
    fig2=d.plotModel(fig=fig2,show=False,color='orange')
    fig2=e.plotModel(fig=fig2,show=False,color='red')
    fig2.text(0.2,0.54,'SMC law',color='black')
    fig2.text(0.2,0.5,'Av=-0.31',color='black')
    fig2.text(0.2,0.46,'Rv=4.3',color='purple')
    fig2.text(0.2,0.42,'Rv=3.9',color='blue')
    fig2.text(0.2,0.38,'Rv=3.5',color='green')
    fig2.text(0.2,0.34,'Rv=3.1',color='yellow')
    fig2.text(0.2,0.30,'Rv=2.7',color='orange')
    fig2.text(0.2,0.26,'Rv=2.4',color='red')
    fig2.show()
    
    #now see how varying Av for a set Rv will affect things.
    
    rv=2.74 #average smc
    av=-0.3
    q=avgsmc(R_V=rv)
    q.EvalCurve(w,av/rv)
    q.UnreddenFlux(f)
    
    av=-0.25
    a=avgsmc(R_V=rv)
    a.EvalCurve(w,av/rv)
    a.UnreddenFlux(f)
    
    av=-0.20
    b=avgsmc(R_V=rv)
    b.EvalCurve(w,av/rv)
    b.UnreddenFlux(f)
    
    av=-0.15
    c=avgsmc(R_V=rv)
    c.EvalCurve(w,av/rv)
    c.UnreddenFlux(f)
    
    av=-0.10
    d=avgsmc(R_V=rv)
    d.EvalCurve(w,av/rv)
    d.UnreddenFlux(f)
    
    av=-0.05
    e=avgsmc(R_V=rv)
    e.EvalCurve(w,av/rv)
    e.UnreddenFlux(f)
    
    fig3=q.plotModel(show=False,color='purple',loglog=True,xlim=(12000,1000))
    fig3=a.plotModel(fig=fig3,show=False,color='blue')
    fig3=b.plotModel(fig=fig3,show=False,color='green')
    fig3=c.plotModel(fig=fig3,show=False,color='yellow')
    fig3=d.plotModel(fig=fig3,show=False,color='orange')
    fig3=e.plotModel(fig=fig3,show=False,color='red')
    fig3.text(0.2,0.54,'SMC law',color='black')
    fig3.text(0.2,0.5,'Rv=2.74',color='black')
    fig3.text(0.2,0.46,'Av=-3.0',color='purple')
    fig3.text(0.2,0.42,'Av=-2.5',color='blue')
    fig3.text(0.2,0.38,'Av=-2.0',color='green')
    fig3.text(0.2,0.34,'Av=-1.5',color='yellow')
    fig3.text(0.2,0.30,'Av=-1.0',color='orange')
    fig3.text(0.2,0.26,'Av=-0.5',color='red')
    
    fig3.show()
    
    
    return b
    

