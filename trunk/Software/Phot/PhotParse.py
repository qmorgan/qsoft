import copy
import matplotlib.pyplot as plt
import numpy as np
from MiscBin.q import maglist2fluxarr
from MiscBin.q import flux2abmag
from MiscBin.q import mag2alpha
from matplotlib import rc
import os
import sys

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
splinedir = os.environ.get("Q_DIR") + '/trunk/Software/Modelling/'
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'


class ObjBlock:
    '''
    Block of ObsBlocks
    '''
    def __init__(self):
        self.obsdict = {}
        self.utburst = None
        self.galebv = None
        self.redshift = None
        
    def updateObj(self,indict):
        if not self.utburst:
            if 'utburst' in indict:
                self.utburst = indict['utburst']
        else:
            if 'utburst' in indict:
                if self.utburst != indict['utburst']:
                    raise ValueError('utburst times do not match!')
        
        if not self.galebv:
            if 'galebv' in indict:
                try:
                    self.galebv = float(indict['galebv'])
                except:
                    self.galebv = None
        else:
            if 'galebv' in indict:
                if self.galebv != float(indict['galebv']):
                    raise ValueError('galebv values do not match!')
                    
        
        if not self.redshift:
            if 'redshift' in indict:
                try:
                    self.redshift = float(indict['redshift'])
                except:
                    self.redshift = None
        else:
            if 'redshift' in indict:
                if self.redshift != float(indict['redshift']):
                    raise ValueError('redshift values do not match!')
        
        name = indict['source'] + '_' + indict['filt']
        if not name in self.obsdict:
            newobs = ObsBlock(indict)
            self.obsdict.update({name:newobs})
        else:
            self.obsdict[name].updateObs(indict)
    
    def CalculateFlux(self):
        '''Loop through all the ObsBlock instances in the ObjBlock and 
        calculate the flux for each, assigning them as attributes of the 
        ObsBlock.
        '''
        from Modelling.ExtinctModel import CorrectFluxForGalExt
        
        for key, ob in self.obsdict.iteritems():
            if ob.filt != None:
                ob.fluxarr,ob.fluxerrarr = \
                    maglist2fluxarr(ob.maglist,ob.magerrlist,ob.filt,singlefilt=True)
                if self.galebv != None:
                    # if we know galebv, do gal corrected flux conversion
                    # do we want to do change this function to handle wave objects? Nah..
                    #array of the same wavelength in angstroms
                    wavearr = np.zeros(len(ob.maglist)) + ob.filt.wave_A 
                    ob.gcfluxarr,ob.gcfluxerrarr = \
                        CorrectFluxForGalExt(self.galebv,wavearr,ob.fluxarr,ob.fluxerrarr)
            else:
                print "No filter for %s, Skipping flux conversion" % (key)
    
    def PlotLC(self):
        # set font
        rc('font', family='Times New Roman')
        
        fig=plt.figure()
        ax=fig.add_axes([0.1,0.1,0.8,0.8])
        ax.loglog()

        for key, ob in self.obsdict.iteritems():
            if ob.filt != None:
                upperinds = np.array(ob.isupperlist)
                detectinds = np.array([not a for a in ob.isupperlist])

                if detectinds.any(): # only plot here if we have at least one detection
                    ax.errorbar(np.array(ob.tmidlist)[detectinds],ob.gcfluxarr[detectinds],yerr=ob.gcfluxerrarr[detectinds], color=ob.color, fmt=ob.marker)
                if upperinds.any(): # only plot here if we have at least one detection
                    ax.errorbar(np.array(ob.tmidlist)[upperinds],ob.gcfluxarr[upperinds],yerr=ob.gcfluxerrarr[upperinds], color=ob.color, fmt='v')
        

        
        old_ylim=ax.get_ylim() # saving for later, as gets modified.. 
        
        ax2=ax.twinx()
        ax3=ax.twiny()
        ax3.loglog()
        
        xobstime=ax.get_xlim()
        xresttime0=xobstime[0]/(1+self.redshift)
        xresttime1=xobstime[1]/(1+self.redshift)
        xrest=(xresttime0,xresttime1)
        ax3.set_xlim(xrest)
        
        
        # duplicate axis for AB mag
        ax.set_ylim(old_ylim)
        ylimflux=ax.get_ylim()
        ylimmag0=flux2abmag(ylimflux[0])
        ylimmag1=flux2abmag(ylimflux[1])
        ylimmag=(ylimmag0,ylimmag1)
        ax2.set_ylim(ylimmag)
        
        # Label the axes
        ax.set_ylabel(r'$F_\nu$ (uJy)')
        zsubscript=str(self.redshift)
        topxlabel = r'$t_{z=%s}$ (s)'  % zsubscript
        ax.set_xlabel(r'$t$ (s)')
        ax2.set_ylabel('AB Mag')
        ax3.set_xlabel(topxlabel)
                
        fig.show()
        
    def SimpleAlphaVSTime(self):
        '''
        Step through every pair of points in a lightcurve to calculate a very
        rough measurement of alpha based on the slope of those two points.
        '''
        # set font
        rc('font', family='Times New Roman')
        
        fig=plt.figure()
        ax=fig.add_axes([0.1,0.1,0.8,0.8])
        ax.semilogx()

        for key, ob in self.obsdict.iteritems():
            upperinds = np.array(ob.isupperlist)
            detectinds = np.array([not a for a in ob.isupperlist])
        
            tmidarr = np.array(ob.tmidlist)[detectinds]
            magarr = np.array(ob.maglist)[detectinds]
            if detectinds.any(): # only plot here if we have at least one detection
                ind = 0
                # loop through each pair of mags 
                while ind < len(magarr) - 1:
                    mag1 = magarr[ind]
                    mag2 = magarr[ind+1]
                    t1 = tmidarr[ind]
                    t2 = tmidarr[ind+1]
                    tmid = (t1 + t2)/2.0
                    alpha = mag2alpha(mag_1=mag1,mag_2=mag2,t_1=t1,t_2=t2)
                    ax.plot(tmid,alpha, color=ob.color, marker=ob.marker)
                    ind += 1
        
        ax.set_ylabel(r'$\alpha$')
        ax.set_xlabel(r'$t_{mid}$ (s)')                
        fig.show()
        
class ObsBlock:
    '''
    Block of Observations of a given Observatory/Filter
    '''
    def __init__(self,indict):
        self.source = indict['source']
        self.filtstr = indict['filt'] 
        self._AssignFilter()
        self._AssignMarker()
        self.maglist=[]
        self.magerrlist=[]
        self.isupperlist=[]
        self.tmidlist=[]
        self.explist=[]
        self.updateObs(indict)
    
    def _AssignMarker(self):
        '''
        Given a source, interpret and assign a marker for plotting purposes
        '''
        if self.source.lower()=='pairitel':
            self.marker='o' # circle
        elif self.source.lower()=='prompt':
            self.marker='d' # diamond
        elif self.source.lower()=='smarts':
            self.marker='s' # square
        elif self.source.lower()=='liverpool':
            self.marker='*' # star
        else:
            self.marker='p' 
            print "unknown source of %s, using default marker" % self.source
    
    def _AssignFilter(self):
        '''
        Given a filtstring, interpret and assign an instance of the filt 
        object'''
        from MiscBin import qObs
        if self.filtstr == 'K' or self.filtstr == 'Ks': 
            self.filt=qObs.Ks
            self.color='#FF9E9E'
        if self.filtstr == 'H': 
            self.filt=qObs.H
            self.color='#FF6969'
        if self.filtstr == 'J': 
            self.filt=qObs.J
            self.color='#FF2929'
        if self.filtstr == "z" or self.filtstr == "z'": 
            self.filt=qObs.z
            self.color='#FF0000'
        if self.filtstr == "i" or self.filtstr == "i'": 
            self.filt=qObs.i
            self.color='#FF8400'
        if self.filtstr == 'I' or self.filtstr == 'Ic': 
            self.color='#FF9D00'
            self.filt=qObs.Ic
        if self.filtstr == "r" or self.filtstr == "r'": 
            self.filt=qObs.r
            self.color='#F5C800'
        if self.filtstr == 'R' or self.filtstr == 'Rc': 
            self.color='#F5C800'
            self.filt=qObs.Rc
        if self.filtstr == 'V': 
            self.color='#9CE805'
            self.filt=qObs.V
        if self.filtstr == "g" or self.filtstr == "g'": 
            self.color='#00D912'
            self.filt=qObs.g
        if self.filtstr == 'B': 
            self.color='#0033FF'
            self.filt=qObs.B
        if self.filtstr == "u" or self.filtstr == "u'": 
            self.color='#9000FF'
            self.filt=qObs.u
        if self.filtstr == 'U': 
            self.filt=qObs.U
            self.color='#9000FF'
        
        if not hasattr(self,'filt'): 
            print '  Could not find appropriate filt object for filtstr %s on source %s; assigning as None' % (self.filtstr,self.source)
            self.filt = None
            self.color='#DDDDDD'
        
    def updateObs(self,indict):
        # update flux/mag values
        self.maglist.append(float(indict['mag']))
        self.magerrlist.append(float(indict['emag']))        
        if 'lim' in indict:
            isupperchar = str(indict['lim']).lower()[0] # will return 'n' if None
            if isupperchar == 'n': self.isupperlist.append(False)
            elif isupperchar == 'y' or isupperchar == 'x': self.isupperlist.append(True) 
            else: raise ValueError('Cannot parse whether upper limit or not!')
        else:
            self.isupperlist.append(False)
            
        # TODO: DO CONVERSIONS
        # find tmid and exp values and perform conversions
        # take first character of inunit/expunit to determine conversion
        # Currently convert everything to seconds
        if indict['inunit'][0] == 'd': inmultfactor = 24*3600.
        elif indict['inunit'][0] == 'h': inmultfactor = 3600.
        elif indict['inunit'][0] == 'm': inmultfactor = 60.
        elif indict['inunit'][0] == 's': inmultfactor = 1.0
        else: raise ValueError('Cannot parse determine mult factor inunit!')
        
        if indict['expunit'][0] == 'd': expmultfactor = 24*3600.
        elif indict['expunit'][0] == 'h': expmultfactor = 3600.
        elif indict['expunit'][0] == 'm': expmultfactor = 60.
        elif indict['expunit'][0] == 's': expmultfactor = 1.0
        else: raise ValueError('Cannot parse determine mult factor for exp!')
        
        
        if 'tmid' in indict and 'exp' in indict: 
        #if tmid is in there, just use that for tmid
            tmid = float(indict['tmid'])*inmultfactor
            exp = float(indict['exp'])*expmultfactor
        elif 'tstart' in indict and 'tend' in indict: 
            # take midtime based on start/end
            tstart = float(indict['tstart'])
            tend = float(indict['tend'])
            tmid = ((tstart+tend)/2.0)*inmultfactor
            if 'exp' not in indict:
                # calculate exposure time if not explicit in indict
                # no conversion necessary since already converted inunit
                exp = (tend-tstart)
            else:
                exp = float(indict['exp'])*expmultfactor
        elif 'tstart' in indict and 'exp' in indict:
            # take midtime based on start+exptime/2
            tstart = float(indict['tstart'])*inmultfactor
            exp = float(indict['exp'])*expmultfactor
            tmid = tstart+exp/2.0
        else:
            errmsg = "Cannot determine tmid and/or exp for %s filter %s" % (self.source,self.filtstr)
            raise ValueError(errmsg)    
        
        self.tmidlist.append(tmid)   
        self.explist.append(exp)                 
    


def SmartInterpolation(obsblock,desired_time_array,errestimate='spline',plot=False):
    '''Will take an obsblock (lightcurve) and desired array of times to 
    interpolate to. 
    
    It will use qSpline optimizing over a number of nodes using GCV methodology
    the model uncertainty given by the spline fit (grey in Figure X) 
    is not the estimated uncertainty on that point had a photometric measurement been made there.  
     instead, this must be added in quadrature with some estimate of what the 
     uncertainty on a measurement would be at that time.  To achieve this, we 
     take the set of all PAIRITEL photometric uncertainties with time and fit
     another spline fit for interpolating to achieve an instrumental uncertainty estimate.
     The final uncertainty on each interpolated point is given by these two added in quadrature.
     '''
    from Modelling.qSpline import qSpline
    from Modelling.qSpline import qSplinePlot
    
    allowed_error_estimates = ['simple','spline']
    if errestimate not in allowed_error_estimates:
        raise ValueError('Please Specify an allowed Error Estimate type')

    logtarray = np.log10(desired_time_array)
    xoutvals = np.array(logtarray)
        
    # for interpolation to work, the tmidlist must be in increasing order
    detectinds = np.array([not a for a in obsblock.isupperlist])
    
    timevals = np.array(obsblock.tmidlist)[detectinds] # get rid of the upper limits 
    yvals = np.array(obsblock.maglist)[detectinds]
    yerrvals = np.array(obsblock.magerrlist)[detectinds]
    
    xvals = np.log10(timevals)
    ylab = r"$m_%s$" % obsblock.filtstr
    newyarr, spline_model_errarr = qSpline(xvals,yvals,yerrvals,xoutvals,plot=False)
    if plot:
        fig=plt.figure(figsize=(8,8))
        ax1=fig.add_axes([0.1,0.4,0.8,0.5])
        qSplinePlot(xvals,yvals,yerrvals,fig=fig,ax_index=0, inverse_y=True,xlabel=r'$t$(s)',ylabel=ylab) #repeat for plot
    newylist = list(newyarr)
    
    # NOW, ESTIMATE THE AVERAGE OBSERVATIONAL ERROR AS A FUNCTION OF TIME
    # TO ADD IN QUADRATURE TO THE MODEL ERROR     
    print list(spline_model_errarr)
    
    if errestimate == 'simple':
        insterrestimate = np.average(yerrvals)
        
    if errestimate == 'spline':
        yerr_errvals = yerrvals*0.1 # Assume 10% error on the errors??
        ylab = r"$m_%s$ err" % obsblock.filtstr
        insterrestimate, error_on_error = qSpline(xvals,yerrvals,yerr_errvals,xoutvals,plot=False)
        if plot:
            ax2=fig.add_axes([0.1,0.1,0.8,0.3])
            qSplinePlot(xvals,yerrvals,yerr_errvals,fig=fig,ax_index=1,xlabel=r'$t$ (s)',ylabel=ylab)
            
    
    if plot:
        #touching up the labels. May just adjust to make margins slightly larger rather than removing ticks.
        ax1lim = ax1.get_ylim()
        ax1adjust = (max(ax1lim)-min(ax1lim))*0.05 #range times a small number
        
        ax2lim = ax2.get_ylim()
        ax2adjust = (max(ax2lim)-min(ax2lim))*0.05
        
        ax1.set_ylim(ax1.get_ylim()[0]+ax1adjust,ax1.get_ylim()[1]-ax1adjust) # backwards due to inverse mag scale
        ax2.set_ylim(0,ax2.get_ylim()[1]+ax2adjust) # bottom axis is zero, cant have - error
        ax1.set_xticks(ax1.get_xticks()[1:-1]) # removing edge xticks for middle plot
        ax1.set_yticks(ax1.get_yticks()[:-1])
        ax2.set_yticks(ax2.get_yticks()[:-1])
        
        filepath = storepath + 'spline' + obsblock.source + '_' + obsblock.filtstr + '.png'
        fig.savefig(filepath)
        
    print insterrestimate       #add uncertainties in quadrature
    spline_model_errlist = list(np.sqrt(spline_model_errarr**2 + insterrestimate**2)) 
    print spline_model_errlist 
    
    obsblock.explist = None
    obsblock.fluxarr=None
    obsblock.fluxerrarr=None
    obsblock.gcfluxerrarr=None
    obsblock.gcfluxarr=None
    obsblock.isupperlist=list(np.ones(len(newylist))==1) # since we forced all upper limits out
    
    obsblock.tmidlist=10**logtarray
    obsblock.maglist=newylist
    obsblock.magerrlist=spline_model_errlist
    return obsblock
    
def DumbInterpolation(obsblock,desired_time_array,fake_error=0.0):
    '''
    if promptr is an obsblock, it has the following attributes:
    In [36]: promptr.
    promptr.color         promptr.fluxarr       promptr.isupperlist   promptr.source
    promptr.explist       promptr.fluxerrarr    promptr.magerrlist    promptr.tmidlist
    promptr.filt          promptr.gcfluxarr     promptr.maglist       promptr.updateObs
    promptr.filtstr       promptr.gcfluxerrarr  promptr.marker
    
    '''
    # for interpolation to work, the tmidlist must be in increasing order
    detectinds = np.array([not a for a in obsblock.isupperlist])
    timevals = np.array(obsblock.tmidlist)[detectinds] # get rid of the upper limits 
    magvals = np.array(obsblock.maglist)[detectinds]
    
    xarray=np.log10(timevals) #log of times is our x array
    assert np.all(np.diff(xarray) > 0)
    yarray = magvals
    logtarray = np.log10(desired_time_array)
    newmaglist = np.interp(logtarray,xarray,yarray,left=0,right=0)
    
    print newmaglist
    logtarray = logtarray[np.nonzero(newmaglist)] # have to have this first before i change newmaglist!
    newmaglist = newmaglist[np.nonzero(newmaglist)] # getting rid of the out of bounds interpolations

    # Figure out some better way to do this..
    errorlist=np.ones(len(newmaglist))
    errorlist*=fake_error
    
    obsblock.explist = None
    obsblock.fluxarr=None
    obsblock.fluxerrarr=None
    obsblock.gcfluxerrarr=None
    obsblock.gcfluxarr=None
    obsblock.isupperlist=list(np.ones(len(newmaglist))==1) # since we forced all upper limits out
    
    obsblock.tmidlist=10**logtarray
    obsblock.maglist=newmaglist
    obsblock.magerrlist=errorlist
    return obsblock
        
def PhotParse(filename,verbose=False):
    object_block = ObjBlock()
    f=file(filename)
    wholefile=f.read() # read in the whole file as a string
    f.close()
    
    # split the string into blocks where there are two line breaks
    strblocks = wholefile.split('\n\n')
    
    headblock=strblocks[0]
    bodyblocks=strblocks[1:]
    
    # set the default for keydict, which will be updated for each block
    default_keydict={'inunit':'sec',
            'expunit':'sec',
            'filt':'unknown',
            'source':'unkown',
            'utburst':'unknown',
            'galebv':'unknown',
            'redshift':'unknown'}
    
    parseable_names=['tmid','tstart','tend','exp','mag','emag','filt','lim']
    name_replace_dict={'filter':'filt',
                        'exptime':'exp',
                        'exposure':'exp',
                        'magnitude':'mag',
                        'limit':'lim',
                        'ulim':'lim',
                        'upper':'lim',
                        'isupper':'lim',
                        't_mid':'tmid',
                        'tstop':'tend'
                        }
    
    # the header should only contain keydict stuff such as: 
    # @inunit=sec
    # @expunit=sec
    # @utburst=04:04:30.21
    # these default values will be used unless otherwise specified in
    # further text blocks.
    for head in headblock.split('\n'):
        if head[0] == '@': #special delimiter denoting default param
            key, val = head[1:].split('=')
            if key not in default_keydict:
                print ' I do not know how to handle key %s, skipping' % (key)
            else:
                default_keydict.update({key:val})
    
    # now loop through the bodyblocks
    for body in bodyblocks:
        if verbose: print "*Moving on to next source...*"
        bodylines = body.split('\n')
        keydict = copy.copy(default_keydict) # set as default but update if @s are defined
        for bodyline in bodylines: # loop through each line in the bodyblock
            if not bodyline: #if blank line, skip
                continue
            if bodyline[0] == '#': #if comment line, skip
                continue 
            #### Grab the HEADERS of the bodyblock to get the default values
            if bodyline[0] == '@': #special delimiter denoting default param
                key, val = bodyline[1:].split('=')
                # replace names as necessary
                if key in name_replace_dict:
                    key = name_replace_dict[key]
                if key not in keydict:
                    print ' I do not know how to handle key %s, skipping' % (key)
                else:
                    keydict.update({key:val})
                    if key == 'source' and verbose:
                        print ' Now reading in data for source %s' % (val)
                    
            #### Grab the FORMAT of the remaining lines
            elif bodyline[0] == '%': #special delimiter of the format line
                fmtlist = bodyline[1:].split()
                #convert to lowercase
                fmtlist = [fmt.lower() for fmt in fmtlist]
                # replace the names
                for name, replacename in name_replace_dict.iteritems():
                    fmtlist = [fmt.replace(name,replacename) for fmt in fmtlist]
                if verbose:
                    # rename the format list if necessary
                    for fmtname in fmtlist:
                        if fmtname not in parseable_names:
                            print "  Cannot parse %s, skipping this column" % fmtname
                        
            
            #### Otherwise it is a line of numbers; split and parse and update relevant ObsBlock
            else:
                datalist = bodyline.strip().split()
                current_data_dict={}
                for fmt in fmtlist:
                    if fmt not in parseable_names:
                        pass
                    else:
                        try:
                            current_val = datalist[fmtlist.index(fmt)]
                        # if the format list is longer than the value list,
                        # we assume its blank lines at the end and assign None
                        except(IndexError): 
                            current_val = None   
                        current_data_dict.update({fmt:current_val})
                
                keydict.update(current_data_dict) # override any defaults with current
                if not 'filt' in keydict and not 'source' in keydict:
                    raise Exeption('Dont have both filt and source!')
                
                # if keydict['source']=='PAIRITEL':
                #     raise Exception
                object_block.updateObj(keydict)
    
    object_block.CalculateFlux()            
    
    return object_block