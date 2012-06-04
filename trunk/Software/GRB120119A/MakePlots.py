from Modelling.ExtinctModel import SEDtimeSimulFit
from Phot import PhotParse
from Modelling.ExtinctModel import _align_SED_times
from Modelling.ExtinctModel import SEDvsTime
from Modelling.ExtinctModel import _getfitdict
from Modelling.Functions import DecayingExponential
import os
import sys
import numpy as np
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt


if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'
paperdir = os.environ.get("Q_DIR") + '/Papers/GRB120119A/'
figuresdir = paperdir + 'Figures/'
tablesdir = paperdir + 'Tables/'

#### INTERPOLATION PLOTS
def _interpplot():
    from GRB120119A import DustModel
    
    objblock=PhotParse.PhotParse('/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119Afinal.dat')
    take_as_given_list = ['PROMPT_I','PROMPT_R', 'PROMPT_V', 'PROMPT_B',"Liverpool_z'","Liverpool_r'",
        'SMARTS_B','SMARTS_V','SMARTS_R','SMARTS_I','SMARTS_J','SMARTS_H','SMARTS_K']
    interpolate_list = ['PAIRITEL_J','PAIRITEL_H','PAIRITEL_K']
    # ptelj=objblock.obsdict['PAIRITEL_J']
    # extraptime=ptelj.tmidlist[0:20]
    # extraptime_r=ptelj.tmidlist[1:20] # special hack for when extrapolating prompt R band
    
    ### BUILD UP A REASONABLE EXTRAPOLTION TIME CODE
    ### BUILD UP AN SED AT EACH POINT IN TIME
    promptR=objblock.obsdict['PROMPT_R']
    promptV=objblock.obsdict['PROMPT_V']
    livz=objblock.obsdict["Liverpool_z'"]
    livr=objblock.obsdict["Liverpool_r'"]
    promptB=objblock.obsdict['PROMPT_B']
    extraptime=promptR.tmidlist[:-10]### SORT THIS
    for time in promptV.tmidlist[0:5]:
        extraptime.append(time)
    for time in livz.tmidlist:
        extraptime.append(time)
    for time in livr.tmidlist:
        extraptime.append(time)
    
    
    objblock=DustModel.BuildInterpolation(objblock,extraptime,interpolate_list,take_as_given_list=None,interp_type='smart')
    raise Exception
    cmd = "mv " + storepath + 'splinePAIRITEL* '+ figuresdir
    os.system(cmd)
    print 'New interpolation plots moved to paper directory.'
    
def _lateSED():
    from Modelling import ExtinctModel
    
    directory = '/Users/amorgan/Data/PAIRITEL/120119A/Combined/120119A_SED.dat'
    a=ExtinctModel.DefaultSEDFit(directory,plotmarg=True)
    cmd = "mv " + storepath + 'SED.png '+ figuresdir + 'lateSED.png'
    os.system(cmd)
    cmd = "mv " + storepath + 'marginalization.png '+ figuresdir + 'SEDmarg.png'
    os.system(cmd)
    print 'New SED plots moved to paper directory.'
    
def _sedvtime():
    from Modelling import ExtinctModel
    ExtinctModel.testSEDvsTime(Av_init=-0.57)