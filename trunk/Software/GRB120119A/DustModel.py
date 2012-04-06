from Modelling.ExtinctModel import SEDtimeSimulFit
from Phot import PhotParse
from Modelling.ExtinctModel import _align_SED_times
from Modelling.ExtinctModel import SEDvsTime
from Modelling.ExtinctModel import _getfitdict
from Modelling.ExtinctModel import DecayingExponential

import numpy as np
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

def ChiSqMap():
    '''
    Loop through a bunch of different parameter values and map out how 
    chisq changes as you change these values.
    '''
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    initial_param='smc'

    time_thresh=10    
    objblock=PhotParse.PhotParse(directory)    
    
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist
    for time in addl_sed_times:
        sedtimelist.append(time)
    
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    objblock=PhotParse.PhotParse(directory)
    
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist
    for time in addl_sed_times:
        sedtimelist.append(time)
    
    objblock=PhotParse.PhotParse(directory)    
    if not sedtimelist: # use this as a default if it is not explicitly defined
        sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist    
    # need to do this to see how many alignments we have    
    aligndict = _align_SED_times(objblock,sedtimelist,time_thresh=time_thresh)
    
    
    # need to do this to see how many alignments we have    
    aligndict = _align_SED_times(objblock,sedtimelist,time_thresh=time_thresh)
    # set up the fit dict            
    fitdict = _getfitdict(initial_param)
    fitdict.pop('const') # get rid of the old const
    fitdict.pop('Av')
    fitdict.pop('beta')
    # add the correct number of normalization constants based on the number of aligned SEDs
    for constnumber in np.arange(len(aligndict)):
        name = 'const' + str(constnumber)
        fitdict.update({name:{'init':10000,'fixed':False}})    
    
    # do a for loop here, fixing Av_1 and beta_1 at different values
    numav=25
    numbeta=25
    steps = numav*numbeta
    Av_1_steps = np.linspace(-2.01,0.01,num=numav) # 50 items from -2 to 0
    beta_1_steps = np.linspace(-1.01,3.01,num=numbeta)
    chisqarr=np.zeros((numav,numbeta))
    count=0
    for avind in np.arange(numav):
        Av_1 = Av_1_steps[avind]
        
        for betaind in np.arange(numbeta):
            count += 1
            print "Doing step %i of %i" % (count,steps)
            
            beta_1 = beta_1_steps[betaind]
            
            fitdict2={
                    'Av_0':{'init':-0.62,'fixed':True},
                    'Av_1':{'init':Av_1,'fixed':True},
                    'Av_2':{'init':1000,'fixed':False},
                    'beta_0':{'init':-1.45,'fixed':True},
                    'beta_1':{'init':beta_1,'fixed':True},            
                    'beta_2':{'init':1000,'fixed':False}
            }
            
            fitdict.update(fitdict2)
            
            try:
                outdict = SEDtimeSimulFit(objblock,sedtimelist,fitdict)
                chisq = outdict['chi2']
            except:
                chisq = np.nan
            
            chisqarr[avind,betaind]=chisq
    return chisqarr
    
    # chisqarr2[np.isnan(chisqarr2)] = -1
    # histeqchisqarr2=ExtinctModel.histeq(chisqarr2)
    # imshow(histeqchisqarr2[0],interpolation='none')
    
    # In [144]: pylab.twinx()
    # Out[144]: <matplotlib.axes.Axes at 0x20ac85b0>
    # 
    # In [145]: pylab.twiny()
    # Out[145]: <matplotlib.axes.Axes at 0x20965b10>
    # 
    # In [149]: pylab.xlim((-1,3))
    # Out[149]: (-1, 3)
    # 
    # In [150]: pylab.ylim((0,-2))
    # Out[150]: (0, -2)
    # 
    # In [167]: pylab.text(0.79,-0.72,'x')
    # Out[167]: <matplotlib.text.Text at 0x20acd850>
    
def SEDtimeSimulFit120119A(initial_param='smc',fixparam='Av', sedtimelist=None,
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    ):
    '''
    time_thresh: Number of seconds we can be off in time from the reference 
    '''
    
    initial_param='smc'

    time_thresh=10    
    objblock=PhotParse.PhotParse(directory)    
    if not sedtimelist: # use this as a default if it is not explicitly defined
        sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist    
    # need to do this to see how many alignments we have    
    aligndict = _align_SED_times(objblock,sedtimelist,time_thresh=time_thresh)

    
    # set up the fit dict            
    fitdict = _getfitdict(initial_param)
    fitdict.pop('const') # get rid of the old const
    fitdict.pop('Av')
    fitdict.pop('beta')
    # add the correct number of normalization constants based on the number of aligned SEDs
    for constnumber in np.arange(len(aligndict)):
        name = 'const' + str(constnumber)
        fitdict.update({name:{'init':10000,'fixed':False}})    
    
    if fixparam == 'beta':
        fitdict2={   
                'Av_0':{'init':-0.62,'fixed':False},
                'Av_1':{'init':-0.3,'fixed':False},
                'Av_2':{'init':300,'fixed':False},
                'beta_0':{'init':-1.45,'fixed':True},
                'beta_1':{'init':0,'fixed':True},            
                'beta_2':{'init':100,'fixed':True}
                }    
    elif fixparam == 'Av':
        fitdict2={
                'Av_0':{'init':-0.62,'fixed':True},
                'Av_1':{'init':0,'fixed':True},
                'Av_2':{'init':100,'fixed':True},
                'beta_0':{'init':-1.30,'fixed':False},
                'beta_1':{'init':-0.2,'fixed':False},            
                'beta_2':{'init':300,'fixed':False}
        }
    elif fixparam == 'both':
        fitdict2={
                'Av_0':{'init':-0.62,'fixed':True},
                'Av_1':{'init':-0.3,'fixed':False},
                'Av_2':{'init':100,'fixed':False},
                'beta_0':{'init':-1.45,'fixed':True},
                'beta_1':{'init':-0.3,'fixed':False},            
                'beta_2':{'init':100,'fixed':False}
        }
    elif fixparam == 'none': #just fit the constants - for testing purposes 
        fitdict2={
                'Av_0':{'init':-0.62,'fixed':True},
                'Av_1':{'init':0,'fixed':True},
                'Av_2':{'init':0,'fixed':True},
                'beta_0':{'init':-1.45,'fixed':True},
                'beta_1':{'init':0,'fixed':True},            
                'beta_2':{'init':0,'fixed':True}
        }
    else: raise ValueError('invalid fixparam')
    fitdict.update(fitdict2)
    
    outdict = SEDtimeSimulFit(objblock,sedtimelist,fitdict)
    
    return outdict


def SEDsimulfit120119Atest(fixparam='beta'):
    # crappy hack way to go about building up the time list if one of the
    # scopes isnt available for each of the SEDs.  Take two timelists, and append
    # the values of one onto the other.  Usually wont have to do this. 
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    objblock=PhotParse.PhotParse(directory)
    
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist
    for time in addl_sed_times:
        sedtimelist.append(time)
    
    
    outdict = SEDtimeSimulFit120119A(sedtimelist=sedtimelist,directory=directory,fixparam=fixparam)
    
    # return outdict 
    for param in outdict['parameters']:
        if param.name == 'Av_1' or param.name == 'beta_1':
            Av1=param.value
        elif param.name == 'Av_2' or param.name == 'beta_2':
            Av2=param.value
        elif param.name == 'Av_0' or param.name ==  'beta_0':
            Av0=param.value
    
    fig=plt.figure()
    ax=fig.add_axes([0.1,0.1,0.8,0.8])        
    ax.semilogx()
    t=np.arange(100000)+1
    # c=BrokenPowerLaw(t,Av0,Av1,Av2)
    c = DecayingExponential(t,Av0,Av1,Av2)
    ax.plot(t,c)
    fig.show()
    return fig


def SEDvsTime120119A():
    # first do fixed Av:
    fig = SEDsimulfit120119Atest(fixparam='Av')
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='grey',fig=fig)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='red',fig=fig)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataSMARTSonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['SMARTS_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='green',fig=fig)
    
    # crappy hack way to go about building up the time list if one of the
    # scopes isnt available for each of the SEDs.  Take two timelists, and append
    # the values of one onto the other.  Usually wont have to do this. 
    # sedtimelist=objblock.obsdict[align_key].tmidlist
    # if second_align_key:
    #     addl_sed_times=objblock.obsdict[align_key].tmidlist
    #     for time in addl_sed_times:
    #         sedtimelist.append(time)
    
    fig.show()
    
    # now fixed beta
    fig = SEDsimulfit120119Atest(fixparam='beta')
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='grey',fig=fig)

    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='red',fig=fig)

    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataSMARTSonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['SMARTS_J'].tmidlist
    fig = SEDvsTime(directory=directory, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='green',fig=fig)
    fig.show()
