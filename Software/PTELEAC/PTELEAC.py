'''Wrapper for the PTEL Early Afterglow Catalog final products'''

import matplotlib
import os
import numpy
from pylab import *
from numpy import array as arr
import cosmocalc
import pandas as pd
from Phot import q_phot
from Phot import PhotParse
from MiscBin import qPickle
import glob

meta_dict = {'061126':{
                        'z':None,
                        'zlim':1.1588,
                        'z':2.346,
                        'zlim':None,
                        'xmin':100,
                        'xmax':2e3,
                        'ymin':100,
                        'ymax':3e4
                        },
            '080310': {
                        'z':2.4266,
                        'zlim':None,
                        'z':2.346,
                        'zlim':None,
                        'xmin':1.5e3,
                        'xmax':3e3,
                        'ymin':6e2,
                        'ymax':2e3
                        },
            '080330':{
                       'z': 1.51,
                       'zlim':None,
                       'z':2.346,
                       'zlim':None,
                       'xmin':1000,
                       'xmax':2e4,
                       'ymin':40,
                       'ymax':2e3
                       },
            '090618':{
                       'z':0.54,
                       'zlim':None,
                       'z':2.346,
                       'zlim':None,
                       'xmin':100,
                       'xmax':3e3,
                       'ymin':1000,
                       'ymax':3e4
                       },
            '080319C':{
                        'z':1.95,
                        'zlim':None,
                        'z':2.346,
                        'zlim':None,
                        'xmin':300,
                        'xmax':2e3,
                        'ymin':100,
                        'ymax':5e3
                        },
            '090530':{
                       'z':None,
                       'zlim':1.6,
                       'z':2.346,
                       'zlim':None,
                       'xmin':100,
                       'xmax':5e3,
                       'ymin':50,
                       'ymax':8e3
                       },
            '080607':{
                       'z':3.036,
                       'zlim':None,
                       'z':2.346,
                       'zlim':None,
                       'xmin':100,
                       'xmax':5e3,
                       'ymin':30,
                       'ymax':3e4
                       },
            '051109A':{
                        'z':2.346,
                        'zlim':None,
                        'xmin':100,
                        'xmax':5e3,
                        'ymin':100,
                        'ymax':1e4
                        },
            '070208':{
                       'z':1.165,
                       'zlim':None,
                       'z':2.346,
                       'zlim':None,
                       'xmin':300,
                       'xmax':2e4,
                       'ymin':50,
                       'ymax':1e3
                       },
            '071025':{
                       'z':4.8,
                       'zlim':None,
                       'z':2.346,
                       'zlim':None,
                       'xmin':100,
                       'xmax':5e3,
                       'ymin':40,
                       'ymax':2e4
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
        globresult = glob.glob(globstr)
        if len (globresult) == 0:
            print "Can't Find any files in {}. External HD plugged in?".format(globresult)
            raise(ValueError)
        pathstr = globresult[0]
        if not os.path.exists(pathstr):
            print "The path {} doesn't appear to exist. Wrong path?".format(pathstr)
        print pathstr
        result = qPickle.load(pathstr)
        galebv=df_grb.loc[GRB]['galebv']
        redshift=df_grb.loc[GRB]['z']
        if outputtext:
            q_phot.textoutput(result,name=GRB,galebv=galebv,redshift=redshift)
        GRB_dict = {GRB:result}
        all_GRB_dict.update(GRB_dict)
    return all_GRB_dict
    
# MAKE PARSEABLE PHOTOMETRY FILES
def make_outdict():
    outdict = makedict(df_grb, outputtext=True)
    return outdict
    
# TRANSPORT TABLES 
#todo

# load tables
catdir = "/Users/amorgan/qsoft/store/ptelcat/"
phottables=glob.glob(catdir+"*data.txt")

def obj_loop():
    for table in phottables:
        phot_obj=PhotParse.PhotParse(table)
        
        GRB_key = phot_obj.name.split(',')[-1] # eg 071025
        GRB_df = df_grb.loc['GRB_key']
        
        phot_obj.PlotLC(show=False,legend=False,xlimits=(100,20000))
        phot_obj.WriteCondensedTable()
    
