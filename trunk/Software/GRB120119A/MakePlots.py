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
        # 78.624,
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
        # 2291.328,
        # 2337.984,
        # 2320.776,  #SMARTS
        # 3523.392,
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
        # 78.624,
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
        # 1373.4599999999998,
        # 1391.04,
        # # 1391.04,
        # 1463.616,
        # 1541.1,
        # 1570.7520000000002,
        # # 1571.6160000000002,
        # 1658.0159999999998,
        # 1673.8799999999999,
        # 1753.0559999999998,
        # # 1753.0559999999998,
        # # 1834.44,
        # 1838.592,
        # 1966.26,
        # 1977.6960000000001,
        # # 1977.6960000000001,
        # 2128.2599999999998,
        # 2196.288,
        # 2260.2000000000003,
        # 2276.64,
        # 2320.776,  #SMARTS
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
def _interpplot():
    # late_time_av=-0.86
    # late_time_beta=-0.90
    
    late_time_av = -0.963 # for DPgrb120119A
    late_time_beta=-1.0
    Av_1init = -1.7
    beta_1init= 2.0
    Av_2init = 70
    beta_2init = 70
    zoomtime = 2000 #seconds
    from GRB120119A import DustModel
    addl_sed_times=None
    testrandom=False
    initial_param = 'DPgrb120119A'
    
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
        time_thresh=5)
     
    cmd = "mv " + storepath + 'SEDvsTime.png '+ figuresdir + 'SEDvsTime_fixedAv_' + initial_param + '.png'
    os.system(cmd)
    print 'New SEDvstime fixed Av plots moved to paper directory.'
    
    # Make SEDvsTime for fixed beta
    SEDvsTime(objblock_interpolated,sedtimelist=sedtimelist,plotsed=False,
        fitlist=['Av'],plotchi2=True,
        Av_init=late_time_av,beta_init=late_time_beta,fig=fig,initial_param=initial_param,
        retfig = False, retchi2=False, color='grey',
        time_thresh=5)
     
    cmd = "mv " + storepath + 'SEDvsTime.png '+ figuresdir + 'SEDvsTime_fixedbeta_' + initial_param + '.png'
    os.system(cmd)
    print 'New SEDvstime fixed Av plots moved to paper directory.'
    
    # Make SEDvsTime for free Av and Beta

    SEDvsTime(objblock_interpolated,sedtimelist=sedtimelist,plotsed=False,
        fitlist=['Av','beta'],plotchi2=False,
        Av_init=late_time_av,beta_init=late_time_beta,initial_param=initial_param,
        retfig = False, retchi2=False, fig=None, color='grey',
        time_thresh=5)
        
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
        sedtimelist=sedtimelist,fixparam='both',
        initial_param=initial_param,
        time_thresh=5,
        plot=True,plotchi2=True,retfig=False,
        Av_0init=late_time_av,beta_0init=late_time_beta,
        Av_1init=Av_1init,beta_1init=beta_1init,
        Av_2init=Av_2init,beta_2init=beta_2init,
        randomize_inits=False,
        unred_latetime=False)
        

    cmd = "mv " + storepath + 'SEDtimesimulfit_' + initial_param + '.png '+ figuresdir
    os.system(cmd)
    
    
    # Marginialization plot 
    qFit.plot_marg_from_fitdict(fitdict,('Av_1','beta_1'))
    cmd = "mv " + storepath + 'marginalization.png '+ figuresdir + 'SEDtimesimulfit_marg.png'
    os.system(cmd)
    print 'New SEDsimulfit plots moved to paper directory.'
    return fitdict
 
def _lateSED():
    from Modelling import ExtinctModel
    initial_param = 'smc'
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119A_SED.dat'
    a=ExtinctModel.DefaultSEDFit(directory,initial_param=initial_param,fitlist=['Av','beta'],plotmarg=True)
    cmd = "mv " + storepath + 'SED.png '+ figuresdir + 'lateSED.png'
    os.system(cmd)
    cmd = "mv " + storepath + 'marginalization.png '+ figuresdir + 'SEDmarg.png'
    os.system(cmd)
    print 'New SED plots moved to paper directory.'
 