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

def ChiSqMap(initial_param='smc',Av_init=-0.62,beta_init=-1.45):
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
    fitdict = _getfitdict(initial_param,Av_init=Av_init,beta_init=beta_init)
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


def BuildInterpolation(objblock,extraptime,interpolate_list,take_as_given_list=None,interp_type='dumb'):
    '''Inputs:
    objblock: block of all relevant observations
    extraptime: list of desired times to extrapolate to
    interpolate_list: obsdict names to interpolate
    take_as_given_list: obsdict names to leave alone. if None, just add all names not in interpolate_list
    
    
    
    interpolate_list = ['PAIRITEL_J','PAIRITEL_H','PAIRITEL_K']
    take_as_given_list = ['PROMPT_I','PROMPT_R', 'PROMPT_V', 'PROMPT_B',"Liverpool_z'","Liverpool_r'",
        'SMARTS_B','SMARTS_V','SMARTS_R','SMARTS_I','SMARTS_J','SMARTS_H','SMARTS_K']
    
    '''
    
    
    newobjblock=PhotParse.ObjBlock() #initialize blank objblock
    newobjblock.utburst = objblock.utburst
    newobjblock.redshift = objblock.redshift
    newobjblock.galebv = objblock.galebv
    
    # taking the times from pairitel and using them as the extrapolation times
    #FIXME - instead, loop through the take_as_given list and extrapolate to those times
    # perhaps allow for some non-extrap times to be included in the sedtime list (SMARTS) 
    # As these are better for the late time values.. worried about contamination
    # from nearby star for small telescopes
    
    # perhaps this could be an attribute of the object block instead.. this would be better 
    # than returning two values.  could have a caveat to align_sed_times or further
    # down the road (i.e. in the SEDfit) to not try to do a fit if not more than 3 or 4 points 

    
    # 30 on 128
    
    if take_as_given_list != None:
        for obs in take_as_given_list:
            obsblock=objblock.obsdict[obs]
            newobjblock.obsdict.update({obs:obsblock})
    else: # just default to adding all keys not in the interpolate list
        for obs, obsblock in objblock.obsdict.iteritems():
            if obs not in interpolate_list:
                newobjblock.obsdict.update({obs:obsblock})
    
    newobjblock.sedtimelist = extraptime
    

    for obs in interpolate_list:
        obsblock = objblock.obsdict[obs]
        if interp_type=='dumb':
            newobsblock=PhotParse.DumbInterpolation(obsblock,extraptime,fake_error=0.1)
        elif interp_type=='smart':
            # if obs == 'PROMPT_R': # cheap hack to fix the extrapolation time for R band since it doesnt go back as far
            #     extraptime2 = extraptime_r
            # else:
            #     extraptime2 = extraptime
            newobsblock=PhotParse.SmartInterpolation(obsblock,extraptime,plot=True)
        newobjblock.obsdict.update({obs:newobsblock})
        assert newobsblock != None
    newobjblock.CalculateFlux()
    return newobjblock 
    
def SEDtimeSimulFit120119A(initial_param='smc',fixparam='Av', sedtimelist=None, defaulttimelist='PTELSMARTS',
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT+SMARTS.dat',
    Av_0init=-0.62,
    Av_1init=0,
    Av_2init=300,
    beta_0init=-1.45,
    beta_1init=0,
    beta_2init=300,
    randomize_inits=False,
    unred_latetime=False, 
    interp_type=None,
    plot=False,
    plotchi2=True    
    ):
    '''
    time_thresh: Number of seconds we can be off in time from the reference 
    
    set unred_latetime=True and initial param to mw for the assumption of late-time SMC and early time MW
    '''


    time_thresh=10    
    objblock=PhotParse.PhotParse(directory)    
    
    if interp_type != None:
        directory='/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat'
        objblock = BuildInterpolation(interp_type=interp_type)
        sedtimelist = objblock.sedtimelist #made sedtimelist an attribute in BuildInterpolation
        # FIXME : need to have a smarter way of building the sedtimelist; should be built up in 
        # the smartinterpolation stage, using all the "take as given" frames
        # sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist[0:20] #only take the first ptel ones
        addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist[0:3] #add the smarts
        for time in addl_sed_times:
            sedtimelist.append(time)
    
    # interpolation will give us an sedtimelist. If we didnt do this, we need to choose it..
    if sedtimelist == None:
        if defaulttimelist == 'PAIRITEL': # use this as a default if it is not explicitly defined
            sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist    
        elif defaulttimelist == 'PTELSMARTS':
            sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist[0:20]#only take the first ptel ones
            addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist[0:3] #add the smarts
            for time in addl_sed_times:
                sedtimelist.append(time)
        else:
            raise ValueError('Invalid defaulttimelist')
    # need to do this to see how many alignments we have    
    aligndict = _align_SED_times(objblock,sedtimelist,time_thresh=time_thresh)

    
    # set up the fit dict            
    fitdict = _getfitdict(initial_param,Av_init=Av_0init,beta_init=beta_0init)
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
    elif fixparam == 'all': # for testing purposes
        fixlist=['Av_0','Av_1','Av_2','beta_0','beta_1','beta_2']
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
    ### END RANDOMIZATION 
    
    
    fitdict.update(fitdict2)
    
    outdict = SEDtimeSimulFit(objblock,sedtimelist,fitdict,correct_late_time_dust=unred_latetime)
    
    
    
    #### FIX THIS - from the fitdict you should be able to get thie initial value 
    # from the ones what were fixed, and those that were free
    # return outdict 
    Av1=None
    Av2=None
    Av0=None
    beta1=None
    beta2=None
    beta0=None
    for param in outdict['parameters']:
        if param.name == 'Av_1':
            Av1=param.value
        elif fitdict['Av_1']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            Av1=fitdict['Av_1']['init']
        if param.name == 'Av_2':
            Av2=param.value
        elif fitdict['Av_2']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            Av2=fitdict['Av_2']['init']
        if param.name == 'Av_0':
            Av0=param.value
        elif fitdict['Av_0']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            Av0=fitdict['Av_0']['init']
        if param.name == 'beta_1':
            beta1=param.value
        elif fitdict['beta_1']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            beta1=fitdict['beta_1']['init']
        if param.name == 'beta_2':
            beta2=param.value
        elif fitdict['beta_2']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            beta2=fitdict['beta_2']['init']
        if param.name == 'beta_0':
            beta0=param.value
        elif fitdict['beta_0']['fixed'] == True: # if we fixed it, wouldnt have gone to parameters
            beta0=fitdict['beta_0']['init']
    
    ####
    #####
    # recalculate chi2 as a function of time based on the model and each time point
    chi2list=[]
    for time in sedtimelist:
        corrtime = time/(1.+objblock.redshift)
        # grab the Av and beta values for each time point
        # note we use the corrected time here, but the observer time in the SEDvsTime,
        # because SEDvsTime only uses it as an index
        Av_snap = DecayingExponential(corrtime,Av0,Av1,Av2)
        beta_snap = DecayingExponential(corrtime,beta0,beta1,beta2)
        
        #hack!
        chi2=SEDvsTime(objblock,initial_param='smc',plotsed=False, fitlist=[],
        sedtimelist=[time],retfig=False,fig=None,plotchi2=False, retchi2=True,
        Av_init=Av_snap,beta_init=beta_snap)
        chi2list.append(chi2)
    print chi2list
    print sum(chi2list)
    corrtimelist = np.array(sedtimelist)/(1.+objblock.redshift)    
    # raise Exception
    ####
    ####
    if plot:
        fig=plt.figure()
        t=np.arange(100000)+1
        
        if not fitdict['beta_1']['fixed'] and not fitdict['Av_1']['fixed']:
            if not plotchi2:
                ax1=fig.add_axes([0.1,0.1,0.8,0.4])
                ax2=fig.add_axes([0.1,0.5,0.8,0.4])
                ax1.semilogx()
                ax2.semilogx()
            else:
                ax2=fig.add_axes([0.1,0.6,0.8,0.3])
                ax1=fig.add_axes([0.1,0.3,0.8,0.3])
                ax3=fig.add_axes([0.1,0.1,0.8,0.2])
                ax3.semilogx()
                ax1.semilogx()
                ax2.semilogx()
            
            c = DecayingExponential(t,Av0,Av1,Av2)
            d = DecayingExponential(t,beta0,beta1,beta2)
            ax1.plot(t,c)
            ax2.plot(t,d)
            ax1.set_ylabel(r'$-1*A_v$')
            ax2.set_ylabel(r'$\beta$')
            ax2.set_xlabel(r'$t$ (s, rest frame)')
            
            if plotchi2:
                ax3.scatter(corrtimelist,chi2list)

                ax2.set_xlabel('')
                ax3.set_xlabel(r'$t$ (s, rest frame)')
                ax3.set_ylabel(r'$\chi^2$')
                
                ylim = ax3.get_ylim()
                ax3.set_ylim([0,ylim[1]]) # ensure bottom is 0; cant have chi2 < 0
                xlim=ax3.get_xlim()# set all xlims to be the same
                
                ax1.set_xlim(xlim)
                ax1.set_xticks(ax1.get_xticks()[1:-1])
                ax2.set_xlim(xlim)
                ax2.set_xticks(ax2.get_xticks()[1:-1])
                
                ax1.set_yticks(ax1.get_yticks()[1:-1])
                ax2.set_yticks(ax2.get_yticks()[1:])
                ax3.set_yticks(ax3.get_yticks()[:-1])
                
                string = 'Total chi2 / dof = %.2f / %i' % (outdict['chi2'],outdict['dof'])
                fig.text(0.55,0.2,string)
            
        elif not fitdict['Av_1']['fixed']:
            ax=fig.add_axes([0.1,0.1,0.8,0.8])        
            ax.semilogx()
            # c=BrokenPowerLaw(t,Av0,Av1,Av2)
            c = DecayingExponential(t,Av0,Av1,Av2)
            ax.plot(t,c)
            
        elif not fitdict['beta_1']['fixed']:
            ax=fig.add_axes([0.1,0.1,0.8,0.8])        
            ax.semilogx()
            # c=BrokenPowerLaw(t,Av0,Av1,Av2)
            c = DecayingExponential(t,beta0,beta1,beta2)
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
            outdict= SEDtimeSimulFit120119A(fixparam='both',randomize_inits=True,plot=False)
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
    fig = SEDtimeSimulFit120119A(fixparam='Av',plot=True)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='grey',fig=fig)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='red',fig=fig)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataSMARTSonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['SMARTS_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['beta'],sedtimelist=sedtimelist, retfig=True, color='green',fig=fig)
    
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
    fig = SEDtimeSimulFit120119A(fixparam='beta',plot=True)
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='grey',fig=fig)

    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataPTELPROMPT.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['PAIRITEL_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='red',fig=fig)

    directory = '/Users/amorgan/Data/PAIRITEL/120119A/PTELDustCompare/AlignedDataSMARTSonly.dat'
    objblock=PhotParse.PhotParse(directory)
    sedtimelist=objblock.obsdict['SMARTS_J'].tmidlist
    fig = SEDvsTime(objblock, plotsed=False,fitlist=['Av'],sedtimelist=sedtimelist, retfig=True, color='green',fig=fig)
    fig.show()

    # Doesn't make sense to plot Av and Beta simultaneously with PhotParse because 
    # these are each independent fits, whereas the simulfit ties the allowed values 
    # of Av and Beta together as a function of time.
    # Perhaps instead we could calculate the Chi2 statistic of flux at each time given
    # the model? Would this be a useful comparison?