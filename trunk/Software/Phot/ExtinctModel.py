import matplotlib.pyplot as plt
import numpy as np
import copy

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
    
    
    def plotModel(self,fig=None,color='blue',show=True):    
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
        
        
        
def avgsmc():
    '''In the average SMC fit, the only free parameters are Av and Beta'''
    R_V = 2.74
    c1 = -4.959
    c2 =  2.264
    c3 = 0.389
    c4 = 0.461
    gamma=1.05
    x0=4.626
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel
    
def avglmc():
    R_V=3.1
    c1=-2.16
    c2=1.31
    c3=1.92
    c4=0.42
    gamma=1.05
    x0=4.626
    fmmodel=FM(R_V,c1,c2,c3,c4,gamma,x0)
    return fmmodel
        
def Test120119A():
    z=1.72
    w=1250. + np.arange(80)*50.
    f = w*0. + 1
    a=avglmc()
    a.EvalCurve(w,-0.1)
    a.UnreddenFlux(f)
    fig = a.plotModel(show=False,color='red')
    b=avgsmc()
    b.EvalCurve(w,-0.1)
    b.UnreddenFlux(f)
    fig=b.plotModel(fig=fig,show=False,color='blue')
    
    
    wv_angst=np.array([  6470. ,   7865.,  12500. ,  16500. ,  21500. ])
    wv_angst_rest = wv_angst/(1+z)
    
    return b
    
