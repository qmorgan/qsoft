from Modelling.ExtinctModel import SEDtimeSimulFit
from Phot import PhotParse
from Modelling.ExtinctModel import _align_SED_times
from Modelling.ExtinctModel import SEDvsTime
from Modelling.ExtinctModel import _getfitdict
from Modelling.Functions import DecayingExponential

import numpy as np
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

def ChiSqMap(initial_param='smc'):
    '''
    Loop through a bunch of different parameter values and map out how 
    chisq changes as you change these values.
    '''
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    

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


def BuildDumbInterpolation():
    objblock=PhotParse.PhotParse('/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat')

    
    newobjblock=PhotParse.ObjBlock()
    newobjblock.utburst = objblock.utburst
    newobjblock.redshift = objblock.redshift
    newobjblock.galebv = objblock.galebv
    
    # taking the times from pairitel and using them as the extrapolation times
    ptelj=objblock.obsdict['PAIRITEL_J']
    extraptime=ptelj.tmidlist[0:20]
    take_as_given_list = ['PAIRITEL_J','PAIRITEL_H','PAIRITEL_K',
        'SMARTS_B','SMARTS_V','SMARTS_R','SMARTS_I','SMARTS_J','SMARTS_H','SMARTS_K']
    for obs in take_as_given_list:
        obsblock=objblock.obsdict[obs]
        newobjblock.obsdict.update({obs:obsblock})
    
    interpolate_list = ['PROMPT_I','PROMPT_R']
    for obs in interpolate_list:
        obsblock = objblock.obsdict[obs]
        newobsblock=PhotParse.DumbInterpolation(obsblock,extraptime,fake_error=0.1)
        newobjblock.obsdict.update({obs:newobsblock})
        assert newobsblock != None
    newobjblock.CalculateFlux()
    return newobjblock
    
def SEDtimeSimulFit120119A(initial_param='smc',fixparam='Av', sedtimelist=None,
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat',
    Av_0init=-0.62,
    Av_1init=0,
    Av_2init=300,
    beta_0init=-1.45,
    beta_1init=0,
    beta_2init=300,
    randomize_inits=False,
    unred_latetime=False, 
    test_interpolation=False
    ):
    '''
    time_thresh: Number of seconds we can be off in time from the reference 
    
    set unred_latetime=True and initial param to mw for the assumption of late-time SMC and early time MW
    '''


    time_thresh=10    
    objblock=PhotParse.PhotParse(directory)    
    
    if test_interpolation == True:
        directory='/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat'
        objblock = BuildDumbInterpolation()

        sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist[0:20] #only take the first ptel ones
        addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist[0:3] #add the smarts
        for time in addl_sed_times:
            sedtimelist.append(time)
    
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
        fitdict.update({name:{'init':20000,'fixed':False}})    
    
    # setting up the new parameters 
    fitdict2={   
            'Av_0':{'init':Av_0init,'fixed':False},
            'Av_1':{'init':Av_1init,'fixed':False},
            'Av_2':{'init':Av_2init,'fixed':False},
            'beta_0':{'init':beta_0init,'fixed':False},
            'beta_1':{'init':beta_1init,'fixed':False},            
            'beta_2':{'init':beta_2init,'fixed':False}
            }
    
    # fixing those as desired
    if fixparam == 'beta':
        fixlist=['beta_0','beta_1','beta_2']
    elif fixparam == 'Av':
        fixlist=['Av_0','Av_1','Av_2']
    elif fixparam == 'both':
        fixlist=['Av_0','beta_0']
    elif fixparam == 'all':
        fixlist=['Av_0','Av_1','Av_2','beta_1','beta_2']
    elif fixparam == 'none':
        fixlist=[]
    else: 
        raise ValueError('invalid fixparam')
    
    for key in fixlist:
        if key in fitdict2:
            fitdict2[key]['fixed']=True
        else:
            errmsg = '%s is not a valid fit parameter'
            raise ValueError(errmsg)
    
    # Go through and randomize the initial parameters, if desired, according to
    # The acceptable physical ranges. Crappy coding practice ahoy!
    if randomize_inits:
        import random
        Av_0range=(0,-2.0)
        Av_1range=(0,-2.0)
        Av_2range=(0,1000)
        beta_0range=(-2.0,2.0)
        beta_1range=(-2.0,2.0)
        beta_2range=(0,1000)
        const_range=(10000,30000)
        if not fitdict2['Av_0']['fixed']:
            fitdict2['Av_0']['init'] = random.uniform(Av_0range[0],Av_0range[1])
        if not fitdict2['Av_1']['fixed']:
            fitdict2['Av_1']['init'] = random.uniform(Av_1range[0],Av_1range[1])
        if not fitdict2['Av_2']['fixed']:
            fitdict2['Av_2']['init'] = random.uniform(Av_2range[0],Av_2range[1])
        if not fitdict2['beta_0']['fixed']:
            fitdict2['beta_0']['init'] = random.uniform(beta_0range[0],beta_0range[1])
        if not fitdict2['beta_1']['fixed']:
            fitdict2['beta_1']['init'] = random.uniform(beta_1range[0],beta_1range[1])
        if not fitdict2['beta_2']['fixed']:
            fitdict2['beta_2']['init'] = random.uniform(beta_2range[0],beta_2range[1])
        for key,value in fitdict.iteritems(): # loop through all the constants
            if key[0:5]=='const':
                value['init'] = random.uniform(const_range[0],const_range[1])
    
    fitdict.update(fitdict2)
    
    outdict = SEDtimeSimulFit(objblock,sedtimelist,fitdict,correct_late_time_dust=unred_latetime)
    
    return outdict




def SEDsimulfit120119Atest(fixparam='beta',randomize_inits=False,plot=True,
    test_interpolation=False):
    # crappy hack way to go about building up the time list if one of the
    # scopes isnt available for each of the SEDs.  Take two timelists, and append
    # the values of one onto the other.  Usually wont have to do this. 
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat'
    objblock=PhotParse.PhotParse(directory)
        
    
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist #add the smarts
    for time in addl_sed_times:
        sedtimelist.append(time)
    
    
    outdict = SEDtimeSimulFit120119A(sedtimelist=sedtimelist,directory=directory,
        randomize_inits=randomize_inits,fixparam=fixparam,test_interpolation=test_interpolation)
    
    # return outdict 
    for param in outdict['parameters']:
        if param.name == 'Av_1' or param.name == 'beta_1':
            Av1=param.value
        elif param.name == 'Av_2' or param.name == 'beta_2':
            Av2=param.value
        elif param.name == 'Av_0' or param.name ==  'beta_0':
            Av0=param.value
    
    if plot:
        fig=plt.figure()
        ax=fig.add_axes([0.1,0.1,0.8,0.8])        
        ax.semilogx()
        t=np.arange(100000)+1
        # c=BrokenPowerLaw(t,Av0,Av1,Av2)
        c = DecayingExponential(t,Av0,Av1,Av2)
        ax.plot(t,c)
        fig.show()
        return fig
    else:
        return outdict

def LoopThroughRandomInits(N=1000):
    '''Testing whether the initial conditions change the final fit value
    as found by leastsq, or whether it is robust against the choices. Loop 
    through N iterations and return the outdict for each.  The acceptable 
    ranges are set in SEDtimeSimulFit120119A.
    '''
    outlist=[]
    faillist=[]
    count=0
    while count < N:
        print 'Now doing Count %i of %i' % (count, N)
        try:
            outdict= SEDsimulfit120119Atest(fixparam='both',randomize_inits=True,plot=False)
            outlist.append(outdict)
        except:
            faillist.append(count)
        count += 1
    print "Fit Failed %i times out of %i" % (len(faillist),N)
    
    chisqlist = []
    for value in outlist:
        chisqlist.append(value['chi2'])
    minchisq = min(chisqlist)
    chisqarr=np.array(chisqlist)
    roundchisqarr=np.round(chisqarr,decimals=2)
    roundminchisq=np.round(minchisq,decimals=2)
    lenfound=len(np.nonzero(roundchisqarr == roundminchisq)[0])
    
    print "The minimum chisq value of %.2f was found %i times." % (roundminchisq,lenfound)
    
    return outlist

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
