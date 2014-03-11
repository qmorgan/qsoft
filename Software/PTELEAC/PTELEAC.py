'''Wrapper for the PTEL Early Afterglow Catalog final products'''

import matplotlib
import os
import numpy
from pylab import *
from numpy import array as arr
import cosmocalc
import pandas as pd
from Phot import q_phot
from MiscBin import qPickle
from glob import glob

meta_dict = {'061126':{
                        'z':None,
                        'zlim':1.1588
                        },
            '080310': {
                        'z':2.4266,
                        'zlim':None
                        },
            '080330':{
                       'z': 1.51,
                       'zlim':None
                       },
            '090618':{
                       'z':0.54,
                       'zlim':None
                       },
            '080319C':{
                        'z':1.95,
                        'zlim':None
                        },
            '090530':{
                       'z':None,
                       'zlim':1.6
                       },
            '080607':{
                       'z':3.036,
                       'zlim':None
                       },
            '051109A':{
                        'z':2.346,
                        'zlim':None
                        },
            '070208':{
                       'z':1.165,
                       'zlim':None
                       },
            '071025':{
                       'z':4.8,
                       'zlim':None
                       }
        }

extinction_grb_list = ['080330','071025','090618','070208','080310','080319C','090709A','051109A','090530','061126','080607']
# average from S & F = Schlafly & Finkbeiner 2011 (ApJ 737, 103) at a nearby calibration star
extinction_list = ['0.0133','0.0635','0.0753','0.0140','0.0336','0.0223','0.0742','0.1553','0.0198','0.1568','0.0191']

df_z= pd.DataFrame(meta_dict).T

dfext=pd.DataFrame({'galebv':extinction_list},index=extinction_grb_list)

df_grb = pd.merge(df_z,dfext,how='outer',left_index=True,right_index=True)


def get_dict(burstlist):
    #Burstlist = update_z_all_GRBs(GRB_list)
    GRB_list =  ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '070208', '071025']
    z_list_bad = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8, 0, 0]
    Burstlist = makedict(GRB_list,outputtext=True)
    return Burstlist
    

def makedict(df_grb, outputtext=False,
    very_good_pickle_path='/Volumes/MyPassport/PTELBACKUP2/picklefiles/very_good_pickles/'):
    ''' make dict for all_bursts() '''
    all_GRB_dict = {}
    
    GRB_list = list(df_grb.index)
    
    for index, GRB in enumerate(GRB_list):
        globstr = very_good_pickle_path + GRB + '*'
        pathstr = glob(globstr)[0]
        print pathstr
        result = qPickle.load(pathstr)
        if outputtext:
            q_phot.textoutput(result,name=GRB)
        GRB_dict = {GRB:result}
        all_GRB_dict.update(GRB_dict)
    return all_GRB_dict