from Modelling.ExtinctModel import SEDtimeSimulFit
from Phot import PhotParse
from Modelling.ExtinctModel import _align_SED_times
from Modelling.ExtinctModel import SEDvsTime
from Modelling.ExtinctModel import _getfitdict
from Modelling.Functions import DecayingExponentialbeta
from Modelling.Functions import DecayingExponentialAv
import os
import sys
import numpy as np
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from Modelling import qFit
import copy
from Modelling import ExtinctModel
from GRB120119A import DustModel

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'
paperdir = os.environ.get("Q_DIR") + '/Papers/GRB120119A/'
figuresdir = paperdir + 'Figures/'
tablesdir = paperdir + 'Tables/'


def _build_extraptime(objblock):
    fullextraptime=False
    justIRextraptime=True
    
    ### BUILD UP A REASONABLE EXTRAPOLTION TIME CODE
    ### BUILD UP AN SED AT EACH POINT IN TIME
    promptI=objblock.obsdict['PROMPT_I']
    promptR=objblock.obsdict['PROMPT_R']
    # promptV=objblock.obsdict['PROMPT_V']
    livz=objblock.obsdict["Liverpool_z'"]
    livi=objblock.obsdict["Liverpool_i'"]
    livr=objblock.obsdict["Liverpool_r'"]
    promptB=objblock.obsdict['PROMPT_B']
    smarts=objblock.obsdict['SMARTS_J']
    print livi.tmidlist
    # raise Exception
    extraptime=promptI.tmidlist[1:]### SORT THIS
    for time in promptR.tmidlist:
        extraptime.append(time)
    # for time in promptV.tmidlist[0:8]:
    #     extraptime.append(time)
    # for time in promptB.tmidlist:
    #     extraptime.append(time)
    # for time in livz.tmidlist:
    #     extraptime.append(time)
    # for time in livr.tmidlist:
    #     extraptime.append(time)
    extraptime.sort()
    print extraptime
    # fullextraptime has promptI, R, V, B, liverpool z, and liverpool r 
    if justIRextraptime: #also smarts maybe
        extraptime = [62.208000000000006,
        # 78.624,
        78.624,
        # 93.312,
        95.04,
        # 118.368,
        118.368,
        # 146.88,
        146.88,
        # 176.256,
        177.12,
        # 202.17600000000002,
        204.768,
        # 234.144,
        239.328,
        # 288.576,
        288.576,
        # 338.688,
        338.688,
        # 388.79999999999995,
        388.79999999999995,
        # 458.784,
        458.784,
        # 550.3679999999999,
        551.232,
        # # 1209.6000000000001, # SUP WITH THIS OUTLIER
        # 1209.6000000000001,
        # 1391.04,
        1391.04,
        # 1570.7520000000002,
        1571.6160000000002,
        # 1752.192,
        1753.0559999999998,
        # 1975.968,
        1977.6960000000001,
        2291.328,
        2337.984,
        2320.776,  #SMARTS
        3523.392,
        3549.312,
        # 4173.9839999999995,
        4178.304,
        # 4877.28,
        4879.872,
        5640.192,
        # # 5641.0560000000005,
        # 6163.776#, #SMARTS
        # # 6494.688,
        # # 6498.144,
        # # 7383.744,
        # # 7437.312,
        # # 8589.024,
        # # 8616.672,
        # # 9838.368,
        # # 9908.352,
        # # 10678.788, #SMARTS
        # # 11163.743999999999,
        # # 11219.039999999999,
        # # 12500.352,
        # # 12553.920000000002,
        # # 14442.624,
        # # 14552.351999999999,
        # # 16687.296000000002,
        # # 16818.624,
        # # 88782.912
    ]
    if fullextraptime:
        extraptime=[62.208000000000006,
        78.624,
        93.312,
        118.368,
        146.88,
        177.12,
        193.98000000000002,
        202.17600000000002,
        239.328,
        # 241.05599999999998,
        276.84,
        288.576,
        # 290.30400000000003,
        338.688,
        # 338.688,
        360.48,
        388.79999999999995,
        # 388.79999999999995,
        444.12,
        458.784,
        # 458.784,
        527.76,
        # 550.3679999999999,
        551.232,
        611.4,
        694.8,
        872.76,
        894.6,
        916.5,
        961.8599999999999,
        984.1200000000001,
        1005.9599999999999,
        1058.1599999999999,
        # # 1209.6000000000001, SUP WITH THIS OUTLIER?
        # # 1209.6000000000001,
        # # 1210.26,
        1373.4599999999998,
        1391.04,
        # 1391.04,
        1463.616,
        1541.1,
        1570.7520000000002,
        # 1571.6160000000002,
        1658.0159999999998,
        1673.8799999999999,
        1753.0559999999998,
        # 1753.0559999999998,
        # 1834.44,
        1838.592,
        1966.26,
        1977.6960000000001,
        # 1977.6960000000001,
        2128.2599999999998,
        2196.288,
        2260.2000000000003,
        2276.64,
        2320.776,  #SMARTS
        # 2337.984,
        # 2455.56,
        # 2587.44,
        # 2719.6800000000003,
        # 2892.54,
        # 3024.42,
        # 3157.32,
        # 3506.1119999999996,
        # 3523.392,
        # 3547.584,
        # 4161.024,
        # 4178.304,
        # # 4182.624,
        # 4879.872,
        # 4889.376,
        # 5444.928000000001,
        # 5577.984,
        # 5641.0560000000005,
        # 6163.776, #SMARTS
        # 6498.144,
        # 7437.312,
        # 8589.024,
        # 9838.368#,
        # 10678.788, #SMARTS
        # # 11219.039999999999,
        # # 12553.920000000002,
        # # 14552.351999999999,
        # # 16818.624
        ]
    return extraptime
     
    
#### INTERPOLATION PLOTS
def _interpplot(initial_param = 'DPgrb120119A',
    late_time_av = -0.963, 
    late_time_beta=-1.0,
    Av_1init = -1.7,
    beta_1init= 2.0,
    Av_2init = 70,
    beta_2init = 70,
    sedvstimeylimdict={"betaonly":None,"Avonly":None,"bothAv":None,"bothbeta":None},
    retfigfortesting=False
    ):
    # sedvstimeylimdict dictionary of ylimits for sedvstime
    # late_time_av=-0.86
    # late_time_beta=-0.90
    

    zoomtime = 2000 #seconds
    addl_sed_times=None
    testrandom=False
    
    
    objblock_original=PhotParse.PhotParse('/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat')
    objblock_original_2 = copy.deepcopy(objblock_original) #not sure why i have do this..
    
    # take_as_given_list = ['PROMPT_I','PROMPT_R', 'PROMPT_V', 'PROMPT_B',"Liverpool_z'","Liverpool_r'",
        # 'SMARTS_B','SMARTS_V','SMARTS_R','SMARTS_I','SMARTS_J','SMARTS_H','SMARTS_K']
    interpolate_list = ['PAIRITEL_J','PAIRITEL_H','PAIRITEL_K']
    # ptelj=objblock.obsdict['PAIRITEL_J']
    # extraptime=ptelj.tmidlist[0:20]
    # extraptime_r=ptelj.tmidlist[1:20] # special hack for when extrapolating prompt R band
    extraptime=_build_extraptime(objblock_original)
    objblock_interpolated=DustModel.BuildInterpolation(objblock_original,extraptime,interpolate_list,
        take_as_given_list=None,interp_type='smart',plot=True)
    cmd = "mv " + storepath + 'splinePAIRITEL* '+ figuresdir
    os.system(cmd)
    print 'New interpolation plots moved to paper directory.'
    # 
    # do a zoomed in interpolation plot
    if zoomtime:
        zoomobjblock=DustModel.BuildInterpolation(objblock_original_2,extraptime,interpolate_list,
            take_as_given_list=None,interp_type='smart',plot=True,plotzoom=zoomtime)
        cmd = "mv " + storepath + 'splinePAIRITEL* '+ figuresdir
        os.system(cmd)
        print 'New interpolation plots moved to paper directory.'

    # Now do the SEDvsTime Plot
    sedtimelist = extraptime
    # addl_sed_times=objblock.obsdict['SMARTS_J'].tmidlist[0:3] #add the smarts
    if addl_sed_times: # shouldnt have to do this; define via extraptimes above
        for time in addl_sed_times:
            sedtimelist.append(time)
    
    if testrandom:
        outlist_rand = DustModel.LoopThroughRandomInits(objblock=objblock_interpolated,
            sedtimelist=sedtimelist,fixparam='both',N=100)
        return outlist_rand
    
    fig = None
    # fig = DustModel.SEDtimeSimulFit120119A(objblock=objblock_interpolated,sedtimelist=sedtimelist,fixparam='Av',plot=True,Av_0init=late_time_av,Av_1init=0,Av_2init=100)
    
    # Make SEDvsTime for fixed Av
    SEDvsTime(objblock_interpolated,sedtimelist=sedtimelist,plotsed=False,
        fitlist=['beta'],plotchi2=True,
        Av_init=late_time_av,beta_init=late_time_beta,fig=fig,initial_param=initial_param,
        retfig = False, retchi2=False, color='grey',
        time_thresh=5,fixylimbeta=sedvstimeylimdict['betaonly'])
     
    cmd = "mv " + storepath + 'SEDvsTime.png '+ figuresdir + 'SEDvsTime_fixedAv_' + initial_param + '.png'
    os.system(cmd)
    print 'New SEDvstime fixed Av plots moved to paper directory.'
    
    # Make g for fixed beta
    SEDvsTime(objblock_interpolated,sedtimelist=sedtimelist,plotsed=False,
        fitlist=['Av'],plotchi2=True,
        Av_init=late_time_av,beta_init=late_time_beta,fig=fig,initial_param=initial_param,
        retfig = False, retchi2=False, color='grey',
        time_thresh=5,fixylimAv=sedvstimeylimdict['Avonly'])
     
    cmd = "mv " + storepath + 'SEDvsTime.png '+ figuresdir + 'SEDvsTime_fixedbeta_' + initial_param + '.png'
    os.system(cmd)
    print 'New SEDvstime fixed Av plots moved to paper directory.'
    
    # Make SEDvsTime for free Av and Beta

    SEDvsTime(objblock_interpolated,sedtimelist=sedtimelist,plotsed=False,
        fitlist=['Av','beta'],plotchi2=True,
        Av_init=late_time_av,beta_init=late_time_beta,initial_param=initial_param,
        retfig = False, retchi2=False, fig=None, color='grey',
        time_thresh=5,fixylimbeta=sedvstimeylimdict['bothbeta'],fixylimAv=sedvstimeylimdict['bothAv'])
        
    cmd = "mv " + storepath + 'SEDvsTime.png '+ figuresdir + 'SEDvsTime_freeAv_' + initial_param + '.png'
    os.system(cmd)
    print 'New SEDvstime plots moved to paper directory.'
    
    
    # Now do the SEDtimeSimulfit plot
    # Av_0init=-0.57,
    # Av_1init=0,
    # Av_2init=100,
    # beta_0init=-1.57,
    # beta_1init=0,
    # beta_2init=100,
    fitdict=DustModel.SEDtimeSimulFit120119A(objblock=objblock_interpolated,
        sedtimelist=sedtimelist,fixparam='both', #fixparam both is Av0 and beta0
        initial_param=initial_param,
        time_thresh=5,
        plot=True,plotchi2=False,retfig=False,
        Av_0init=late_time_av,beta_0init=late_time_beta,
        Av_1init=Av_1init,beta_1init=beta_1init,
        Av_2init=Av_2init,beta_2init=beta_2init,
        randomize_inits=False,
        unred_latetime=False,
        plot_every_model=True)
        
    if retfigfortesting:
        fig=DustModel.SEDtimeSimulFit120119A(objblock=objblock_interpolated,
            sedtimelist=sedtimelist,fixparam='both',
            initial_param=initial_param,
            time_thresh=5,
            plot=True,plotchi2=False,retfig=True,
            Av_0init=late_time_av,beta_0init=late_time_beta,
            Av_1init=Av_1init,beta_1init=beta_1init,
            Av_2init=Av_2init,beta_2init=beta_2init,
            randomize_inits=False,
            unred_latetime=False)
        return fig
        
    cmd = "mv " + storepath + 'SEDtimesimulfit.png '+ figuresdir + 'SEDtimesimulfit_' + initial_param + '.png'
    os.system(cmd)
    
    
    # Marginialization plot 
    qFit.plot_marg_from_fitdict(fitdict,('Av_1','beta_1'))
    cmd = "mv " + storepath + 'marginalization.png '+ figuresdir + 'SEDtimesimulfit_' + initial_param + '_marg.png'
    os.system(cmd)
    print 'New SEDsimulfit plots moved to paper directory.'
    return fitdict
 
def _lateSED(initial_param = 'smc',incl_xray=False):
    from Modelling import ExtinctModel
    if incl_xray:
        addl_str = '_xrt'
        addl = 'xray'
    else:
        addl_str = ''
        addl = ''
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119A_SED.dat'
    a=ExtinctModel.DefaultSEDFit(directory,initial_param=initial_param,fitlist=['Av','beta'],
        plotmarg=True, incl_xray=incl_xray)
    
    cmd = "mv " + storepath + 'SED' + addl + '.png '+ figuresdir + 'lateSED_' + initial_param + addl_str +  '.png'
    os.system(cmd)
    cmd = "mv " + storepath + 'marginalization.png '+ figuresdir + 'SEDmarg_' + initial_param + addl_str + '.png'
    os.system(cmd)
    print 'New SED plots moved to paper directory.'

def _lateSEDx2(initial_param = 'smc'):
    _lateSED(initial_param=initial_param,incl_xray=False)
    _lateSED(initial_param=initial_param,incl_xray=True)
    

def _lightcurves():
    objblock_original=PhotParse.PhotParse('/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat')
    objblock_xrt = PhotParse.PhotParse('/Users/amorgan/Data/PAIRITEL/120119A/Combined/xrt/xrtlcbin_test.dat')
    
    objblock_original.PlotLC(show=True,save=True,legend=True,residualscale=False,
        xlimits=(3e1,1e5),ylimits=(1e0,3e4))
    cmd = "mv " + storepath + 'LC_GRB120119A.png '+ figuresdir
    os.system(cmd)
    
    objblock_xrt.PlotXRTlc(show=True,save=True,legend=True,
        obslist=['BAT_unknown','XRT_unknown','PAIRITEL_K'],
        xlimits=(1e0,1e6),ylimits=None)

    cmd = "mv " + storepath + 'LCxrt_GRB120119A.png '+ figuresdir
    os.system(cmd)
    print 'New lightcurve plots moved to paper directory.'

def _make_proposal_figure():
    modeldict={
        'DPgrb120119Axrt':{
            'dust_model':'DPgrb120119Axrt',
            'xrt_incl':True,
            'late_time_av':-1.00, #NEED TO RECHECK THESE VALUES
            'late_time_beta':-0.92, #NEED TO RECHECK THESE VALUES
            'Av_1init':-1.7,
            'beta_1init':2.0,
            'Av_2init':70,
            'beta_2init':70,
            'sedvstimeylimdict':{"betaonly":(-1.9,-0.8),"Avonly":(0.95,1.79),"bothAv":(0,6),"bothbeta":(-2.5,4.5)}
            }
        }

    for paramlist in modeldict.itervalues():
        myfig = _interpplot(initial_param = paramlist['dust_model'],
            late_time_av = paramlist['late_time_av'], 
            late_time_beta= paramlist['late_time_beta'],
            Av_1init = paramlist['Av_1init'],
            beta_1init= paramlist['beta_1init'],
            Av_2init = paramlist['Av_2init'],
            beta_2init = paramlist['beta_2init'],
            sedvstimeylimdict=paramlist['sedvstimeylimdict'],
            retfigfortesting=True
            ) 
            
        
        ax1=myfig.get_axes()[0]
        ax2=myfig.get_axes()[1]
        ax2.set_ylabel(r'$\beta$',size=20)
        ax1.set_xlabel('$t$ (s, rest frame)',size=20)
        ax1.set_ylabel('$A_V$',size=20)
        ax2.set_xticks([])
        # ax1.set_ylim(0.95,1.79)
        
        ax1.set_xlim(10,1e4)
        ax2.set_xlim(10,1e4)
        
        # beta_1: -0.24 +/- 0.07
        # beta_2: 799.48 +/- 289.05
        # Av_1: -0.72 +/- 0.10
        # Av_2: 53.62 +/- 13.48
        myfig.text(0.6,0.28,r'$\Delta A_V = 0.72 \pm 0.1$',size=20)
        myfig.text(0.6,0.35,r'$\tau_{A_V} = 53 \pm 13$ s',size=20)
        myfig.text(0.6,0.55,r'$\Delta \beta = -0.24 \pm 0.07$',size=20)
        myfig.text(0.6,0.62,r'$\tau_{\beta} = 800 \pm 290$ s',size=20)
        
        return myfig

def _make_colorchange_table():
    contentlist=[]
    
    modeldict={
        # 'smc':{
            # 'dust_model':'smc',
            # 'xrt_incl':False,
            # 'late_time_av':-0.70, 
            # 'late_time_beta':-1.21,
            # 'Av_1init':-1.7,
            # 'beta_1init':2.0,
            # 'Av_2init':70,
            # 'beta_2init':70
            # },
        'smcx':{
            'dust_model':'smc',
            'xrt_incl':True,
            'late_time_av':-0.85, 
            'late_time_beta':-0.90,
            'Av_1init':-1.7,
            'beta_1init':2.0,
            'Av_2init':70,
            'beta_2init':70,
            'sedvstimeylimdict':{"betaonly":(-1.8,-0.8),"Avonly":(0.7,1.5),"bothAv":(0,2.5),"bothbeta":(-2.0,0.9)}
            },
        'lmc2':{
            'dust_model':'lmc2',
            'xrt_incl':True,
            'late_time_av':-1.15, 
            'late_time_beta':-0.94,
            'Av_1init':-1.7,
            'beta_1init':2.0,
            'Av_2init':70,
            'beta_2init':70,
            'sedvstimeylimdict':{"betaonly":(-1.7,-0.7),"Avonly":(0.9,1.8),"bothAv":(0,3),"bothbeta":(-3,0.9)}
            },
        'DPgrb120119A':{
            'dust_model':'DPgrb120119A',
            'xrt_incl':False,
            'late_time_av':-0.81, # NEED TO RECHECK THESE VALUES 
            'late_time_beta':-1.26, # NEED TO RECHECK THESE VALUES
            'Av_1init':-1.7,
            'beta_1init':2.0,
            'Av_2init':70,
            'beta_2init':70,
            'sedvstimeylimdict':{"betaonly":(-2.0,-1.0),"Avonly":(0.8,1.5),"bothAv":(0,6.0),"bothbeta":(-2.5,4.5)}
            },
        'DPgrb120119Axrt':{
            'dust_model':'DPgrb120119Axrt',
            'xrt_incl':True,
            'late_time_av':-1.00, #NEED TO RECHECK THESE VALUES
            'late_time_beta':-0.92, #NEED TO RECHECK THESE VALUES
            'Av_1init':-1.7,
            'beta_1init':2.0,
            'Av_2init':70,
            'beta_2init':70,
            'sedvstimeylimdict':{"betaonly":(-1.9,-0.8),"Avonly":(0.9,2.0),"bothAv":(0,6),"bothbeta":(-2.5,4.5)}
            }
        }
            
    for paramlist in modeldict.itervalues():
        fitdict = _interpplot(initial_param = paramlist['dust_model'],
            late_time_av = paramlist['late_time_av'], 
            late_time_beta= paramlist['late_time_beta'],
            Av_1init = paramlist['Av_1init'],
            beta_1init= paramlist['beta_1init'],
            Av_2init = paramlist['Av_2init'],
            beta_2init = paramlist['beta_2init'],
            sedvstimeylimdict=paramlist['sedvstimeylimdict'],
            ) 
        
        beta_1str = 'error'
        Av_1str = 'error'
        taubetastr = 'error'
        tauAvstr = 'error'
            
        for outstr in fitdict['strings']: #assume beta_0 and Av_0 are fixed and thus not in strin list
            if 'beta_1:' in outstr:
                beta_1str = outstr.lstrip('beta_1:')
            if 'beta_2:' in outstr:
                taubetastr = outstr.lstrip('beta_2:')
            if 'Av_1:' in outstr:
                if outstr.find('Av_1: -') == 0:
                    outstr=outstr.replace('Av_1: -','Av_1: ') # since Av is negative what it should be in this code
                else:
                    outstr=outstr.replace('Av_1: ','Av_1: -') #replacing in case the fit actually shows the best fit Av IS negative                    
                Av_1str = outstr.lstrip('Av_1:')
            if 'Av_2:' in outstr:
                tauAvstr = outstr.lstrip('Av_2:')
        
        Av_0str = str(paramlist['late_time_av']) + ' (fixed)'
        beta_0str = str(paramlist['late_time_beta']) + ' (fixed)'
        
        newstring = '%s & %s & $%s$ & $%s$ & %s & $%s$ & $%s$ & %.1f / %i \\\\' % (paramlist['dust_model'],Av_0str,Av_1str.replace("+/-","\pm"),tauAvstr.replace("+/-","\pm"),beta_0str,beta_1str.replace("+/-","\pm"),taubetastr.replace("+/-","\pm"),fitdict['chi2'],fitdict['dof'])
        contentlist.append(newstring)
    

    header='''
\\begin{deluxetable}{llllllll}
\\tablecaption{Results of Color Change Fits}
\\tabletypesize{\scriptsize}
\\tablewidth{0pt}
\\tablehead{
\\colhead{Dust} & \\colhead{$A_{V,0}$} & \\colhead{$A_{V,1}$} & \\colhead{$\\tau_{A_V}$} & \\colhead{$\\beta_{0}$} & \\colhead{$\\beta_{1}$} & \\colhead{$\\tau_{\\beta}$} & \\colhead{$\chi^2$ $/$ dof} \\\\
\\colhead{Model}           & \\colhead{(mag)}    & \\colhead{(mag)}        & \\colhead{(s)}    & \\colhead{} & \\colhead{} & \\colhead{(s)} & \\colhead{}}
\\startdata

'''

    footer= '''
\\enddata
\\tablecomments{Results of color change model fits to the interpolated early-time SEDs of GRB~110119A.}
\\label{tab:ccfits}
\\end{deluxetable}
    '''
    
    
    content = ''
    for contentline in contentlist:
        content += ''' %s
''' % contentline
    
    tabletext = header + content + footer
    tabletext = tabletext.replace("DPGRB120119AXRT","FMX")
    tabletext = tabletext.replace("DPGRB120119A","FM")
    try:
        filename = tablesdir + 'colorchangetable.tex'
        f = open(filename,'w')
        f.write(tabletext)
        f.close()
        
        print ''
        print "Wrote colorchange table to tables directory"
    except:
        print "FAILED TO WRITE COLOR CHANGE TABLE" 
    

def _make_SED_table():
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119A_SED.dat'
    
    
    paramlist=['mw','lmc','lmc2','smc','DPgrb120119Axrt','DPgrb120119A']
    inclxrt=[True,True,True,True,True,False]
    
    contentlist = []
        
    for param in paramlist:
        
        if inclxrt[paramlist.index(param)]:
            incllist = [False,True]
        else:
            incllist = [False]
        for incl in incllist:    # repeat the fit if including XRT
            fitdict=ExtinctModel.DefaultSEDFit(directory,initial_param=param,fitlist=['Av','beta'],
            plotmarg=False, plot=False, incl_xray=incl)
            for outstr in fitdict['strings']:
                if 'beta:' in outstr:
                    betastr = outstr.lstrip('beta:')
                if 'Av:' in outstr:
                    if outstr.find('Av: -') == 0:
                        outstr=outstr.replace('Av: -','Av: ') # since Av is negative what it should be in this code
                    else:
                        outstr=outstr.replace('Av: ','Av: -') #replacing in case the fit actually shows the best fit Av IS negative                    
                    Avstr = outstr.lstrip('Av:')    
            
            
            if incl:
                inclXstr = 'Y'
            else:
                inclXstr = 'N'

            newstring = '%s & %s & $%s$ & $%s$ & %.1f / %i \\\\' % (param.upper(),inclXstr,betastr.replace("+/-","\pm"),Avstr.replace("+/-","\pm"), fitdict['chi2'],fitdict['dof'])
            contentlist.append(newstring)


    content = ''
    for contentline in contentlist:
        content += ''' %s
''' % contentline
        

    header='''
\\begin{deluxetable}{lllll}
\\tablecaption{Results of Extinction Fits}
\\tablewidth{0pt}
\\tablehead{
\\colhead{Dust} & \\colhead{+XRT?} & \\colhead{$\\beta$} & \\colhead{$A_V$}  & \\colhead{$\chi^2$ $/$ dof} \\\\
\\colhead{Model}           & \\colhead{}    & \\colhead{}        & \\colhead{(mag)}    & \\colhead{}}
\\startdata

'''

    footer= '''
\\enddata
\\tablecomments{Results of standard dust model fits to the SED of GRB~110119A from the contemporaneous SMARTS data at $t=38.7$ minutes after the burst. [Converged on negative extinction for MW - what is standard practice here? Say just $<0.02$? Leave as-is? Rewire modelling to disallow negative Av?]}
\\label{tab:extfits}
\\end{deluxetable}
    '''
    
    tabletext = header + content + footer
    tabletext = tabletext.replace("DPgrb120119Axrt","FMX")
    tabletext = tabletext.replace("DPgrb120119A","FM")
    try:
        filename = tablesdir + 'dusttable.tex'
        f = open(filename,'w')
        f.write(tabletext)
        f.close()
        
        print ''
        print "Wrote dust table to tables directory"
    except:
        print "FAILED TO WRITE TABLE" 
