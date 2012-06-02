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
    DustModel.BuildInterpolation(interp_type='smart')
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
    