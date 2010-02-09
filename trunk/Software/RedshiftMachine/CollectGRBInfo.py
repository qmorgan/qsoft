#!/usr/bin/env python
# encoding: utf-8
"""
CollectGRBInfo.py
Author: Adam N Morgan
Created: Aug 24, 2009
Last Updated: Aug 24, 2009

This program collects all the parsed info from ParseSwiftCat.py, 
ParseGCNNotice.py, and ParseNatCat.py into a single place and has the option
to output a .arff file for machine learning with Weka.
"""

import sys
import os
import time
from RedshiftMachine import ParseSwiftCat
from RedshiftMachine import LoadGCN
from AutoRedux import Signal
from MiscBin.q import object2dict
from MiscBin import qPickle
import pylab
import scipy
import numpy

#from MiscBin.q import Standardize

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def_bat_natcats = [storepath+'bat_catalog_07061275.fits',storepath+'bat_catalog_current.fits']
def_xrt_natcats = [storepath+'xrt_catalog_090831.fits']


def whatis(keyword):
    '''Look up what a particular keyword represents in the collected dictionary.
    
    example values in the helpdict is for '060908' except Beta, which is from 081118
    
    '''

    helpdict = {
    'A': {'definition':'Alpha','type':'double','source':'NatBat Spectra','sample':-0.77270000000000005},
    'A1': {'definition':'Alpha Lower limit','type':'double','source':'NatBat Spectra','sample':-1.0469999999999999},
    'A2': {'definition':'Alpha Upper limit','type':'double','source':'NatBat Spectra','sample':-0.48049999999999998},
    'B': {'definition':'Beta','type':'double','source':'NatBat Spectra','sample':-2.1383000000000001},
    'B1': {'definition':'Beta Lower limit','type':'double','source':'NatBat Spectra','sample':-2.3252999999999999},
    'B2': {'definition':'Alpha Upper limit','type':'double','source':'NatBat Spectra','sample':-1.9612000000000001},
    'CHI2': {'definition':'BAT Spectral Fit Chi-squared','type':'double','source':'NatBat Spectra','sample':36.090000000000003},
    'CHI2_PC': {'definition':'default definition','type':'string','source':'NatBat','sample':'31.90/44'},
    'CHI2_PC_LATE': {'definition':'default definition','type':'string','source':'NatBat','sample':'17.49/17'},
    'CHI2_WT': {'definition':'default definition','type':'string','source':'NatBat','sample':'17.85/20'},
    'DPK_O_CTS': {'definition':'Peak Countrate over Counts Uncertainty (s^-1)','type':'string','source':'NatBat','sample':0.0084310400000000008},
    'DRT45': {'definition':'rT_0.45 Uncertainty (s)','type':'string','source':'NatBat Timing','sample':0.10100000000000001},
    'DT50': {'definition':'T50 Uncertainty (s)','type':'string','source':'NatBat Timing','sample':0.22600000000000001},
    'DT90': {'definition':'T90 Uncertainty (s)','type':'string','source':'NatBat Timing','sample':0.16900000000000001},
    'DT_MAX_SNR': {'definition':'Duration of the time window containing the maximum s/n detection (s)','type':'string','source':'NatBat','sample':11.220000000000001},
    'EISO': {'definition':'Bayes Isotropic equivalent Energy [erg] (Approximate Bolometric Fluence [erg/cm^2] in the 1-10^4 keV band if no Redshift)','type':'double','source':'NatBat Spectra','sample':6.9099999999999996e+52},
    'EISO1': {'definition':'Bayes E_iso lower limit','type':'double','source':'NatBat','sample':5.6600000000000002e+52},
    'EISO2': {'definition':'Bayes E_iso upper limit','type':'double','source':'NatBat','sample':1.0900000000000001e+53},
    'EP': {'definition':'Bayes E_peak','type':'double','source':'NatBat Spectra','sample':161.47},
    'EP0': {'definition':'E_peak','type':'double','source':'NatBat Spectra','sample':145.96340000000001},
    'EP01': {'definition':'E_peak lower limit','type':'double','source':'NatBat Spectra','sample':112.1662},
    'EP02': {'definition':'E_peak upper limit','type':'double','source':'NatBat Spectra','sample':253.4169},
    'EP1': {'definition':'Bayes E_peak lower limit ','type':'double','source':'NatBat Spectra','sample':122.63},
    'EP2': {'definition':'Bayes E_peak upper limit','type':'double','source':'NatBat Spectra','sample':334.51999999999998},
    'FL': {'definition':'Energy Fluence (15-350 keV) [erg/cm^2]','type':'double','source':'NatBat Spectra','sample':4.1281000000000003e-06},
    'FL1': {'definition':'Energy Fluence (15-350 keV) lower limit','type':'double','source':'NatBat Spectra','sample':3.7743999999999998e-06},
    'FL2': {'definition':'Energy Fluence (15-350 keV) upper limit','type':'double','source':'NatBat Spectra','sample':4.5635999999999998e-06},
    'FLX_PC': {'definition':'Average Flux from XRT PC Data (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':6.1600000000000007e-05},
    'FLX_PC1': {'definition':'Avg Flux from PC Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':5.5609999999999998e-05},
    'FLX_PC2': {'definition':'Avg Flux from PC Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':6.8150000000000003e-05},
    'FLX_PC_LATE': {'definition':'Average Flux from Late Time (>10ks) XRT PC Data (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':8.8419999999999994e-06},
    'FLX_PC_LATE1': {'definition':'Avg Flux from Late Time PC Data Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':5.5670000000000002e-06},
    'FLX_PC_LATE2': {'definition':'Avg Flux from Late Time PC Data Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':1.429e-05},
    'FLX_WT': {'definition':'Average Flux from XRT WT Data (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':0.0060889999999999998},
    'FLX_WT1': {'definition':'Avg Flux from WT Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':0.0052449999999999997},
    'FLX_WT2': {'definition':'Avg Flux from WT Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','sample':0.007038},
    'GAM_PC': {'definition':'Gamma inferred from XRT PC Data','type':'double','source':'NatXRT','sample':1.8776999999999999},
    'GAM_PC1': {'definition':'PC Gamma lower Limit','type':'double','source':'NatXRT','sample':1.7541},
    'GAM_PC2': {'definition':'PC Gamma upper limit','type':'double','source':'NatXRT','sample':2.0095999999999998},
    'GAM_PC_LATE': {'definition':'Gamma inferred from late time (>10ks) XRT PC Data','type':'double','source':'NatXRT','sample':1.9380999999999999},
    'GAM_PC_LATE1': {'definition':'Late time PC Gamma lower Limit','type':'double','source':'NatXRT','sample':1.2849999999999999},
    'GAM_PC_LATE2': {'definition':'Late time PC Gamma upper Limit','type':'double','source':'NatXRT','sample':2.6629999999999998},
    'GAM_WT': {'definition':'Gamma inferred from XRT WT Data','type':'double','source':'NatXRT','sample':2.3073999999999999},
    'GAM_WT1': {'definition':'WT Gamma lower Limit','type':'double','source':'NatXRT','sample':2.0787},
    'GAM_WT2': {'definition':'WT Gamma upper Limit','type':'double','source':'NatXRT','sample':2.5615999999999999},
    'MAX_SNR': {'definition':'Maximum BAT S/N Ratio','type':'double','source':'NatBat','sample':47.399999999999999},
    'MODEL': {'definition':'BAT Spectral Fit Model Used for Frequentist approach: Either Simple Power Law Model (PLM), powerlaw times an exponetial cutoff (PLEXP), or smoothly connected broken powerlaw - the Band GRB Model (GRBM)','type':'string','source':'NatBat Spectra','sample':'PLEXP'},
    'NH_GAL': {'definition':'Galactic Neutral Hydrogen Column at GRB Location (cm^-2)','type':'double','source':'NatXRT','sample':0.023400000000000001},
    'NH_PC': {'definition':'Excess N_H Column inferred from XRT PC Data (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.017999999999999999},
    'NH_PC1': {'definition':'Excess PC N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':-0.0080000000000000002},
    'NH_PC2': {'definition':'Excess PC N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.050000000000000003},
    'NH_PC_LATE': {'definition':'Excess N_H Column inferred from late time (>10ks) XRT PC Data (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.057000000000000002},
    'NH_PC_LATE1': {'definition':'Excess late time PC N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':-0.128},
    'NH_PC_LATE2': {'definition':'Excess late time PC N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.27500000000000002},
    'NH_WT': {'definition':'Excess N_H Column inferred from XRT WT Data (10^22 cm^-2)','type':'double','source':'NatXRT Spectra','sample':0.052999999999999999},
    'NH_WT1': {'definition':'Excess WT N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.0070000000000000001},
    'NH_WT2': {'definition':'Excess WT N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','sample':0.106},
    'NISO': {'definition':'Isotropic equivalent Photons [ph] (approximate bolometric photon fluence [ph/cm^-2] in observer frame 1-10^4 keV band if no redshift)','type':'double','source':'NatBat','sample':3.4699999999999998e+59},
    'NISO1': {'definition':'default definition','type':'double','source':'NatBat','sample':2.61e+59},
    'NISO2': {'definition':'default definition','type':'double','source':'NatBat','sample':5.0600000000000004e+59},
    'NU': {'definition':'BAT Spectral Fit D.O.F.','type':'double','source':'NatBat Spectra','sample':54.0},
    'OBS': {'definition':'Swift Trigger Number','type':'string','source':'NatBat','sample':'00228581'},
    'PK_O_CTS': {'definition':'Ratio of the peak rate Rate_p (in a time bin of width 0.01 dt_(S/N)) over the total source counts (cts).  Used to approximately relate the burst fluences to peak fluxes.','type':'double','source':'NatBat','sample':0.114758},
    'RT45': {'definition':'BAT rT0.45 (Defined in Reichart et al. 2001)','type':'double','source':'NatBat Timing','sample':5.2800000000000002},
    'T0': {'definition':'Lower BAT Spectral Time Region','type':'string','source':'NatBat Spectra','sample':'-9.615'},
    'T1': {'definition':'Upper BAT Spectral Time Region','type':'string','source':'NatBat Spectra','sample':'13.815'},
    'T50': {'definition':'BAT T50 Difference between 75th and 25th percentile time of total counts relative to start of burst interval','type':'double','source':'NatBat Timing','sample':7.2599999999999998},
    'T90': {'definition':'BAT T90 Difference between 95th and 5th percentile time of total counts relative to start of burst interval - loosely consistant with Swift T90; Highly dependent on burst start & stop times','type':'double','source':'NatBat Timing','sample':18.48},
    'T_PC': {'definition':'XRT Photon Counting Mode Time Region (ks post burst)','type':'string','source':'NatXRT Spectra','sample':'0.156101-1087.47'},
    'T_PC_LATE': {'definition':'XRT Late Photon Counting Mode Time Region (at least 10ks post burst)','type':'string','source':'NatXRT Spectra','sample':'10.0000-1087.4741'},
    'T_WT': {'definition':'XRT WindowTiming Mode Region Time Region (ks post burst)','type':'string','source':'NatXRT Spectra','sample':'0.0802510-5.72213'},
    'UT': {'definition':'UT Time of Burst','type':'string','source':'NatBat','sample':'20060908_085722.340000'},
    'Z': {'definition':'Redshift (from Nat)','type':'double','source':'NatBat Spectra','sample':2.4300000000000002},
    'b_mag_isupper': {'definition':'Is the UVOT B Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'bat_bkg_dur': {'definition':'GCN BAT Background duration (seconds)','type':'double','source':'Swift-BAT GRB Position','sample':8.0},
    'bat_bkg_inten': {'definition':'GCN BAT Background intensity (counts)','type':'double','source':'Swift-BAT GRB Position','sample':41524.0},
    'bat_dec': {'definition':'Declination of burst, as determined by Swift BAT','type':'double','source':'Swift-BAT GRB Position','sample':0.37},
    'bat_img_peak': {'definition':'GCN GRB BAT Image Peak (image_counts)','type':'double','source':'Swift-BAT GRB Position','sample':240.0},
    'bat_inten': {'definition':'GCN GRB BAT Intensity (counts)','type':'double','source':'Swift-BAT GRB Position','sample':6439.0},
    'bat_is_rate_trig': {'definition':'Was the event a Rate trigger (yes) or image trigger (no)','type':'string','source':'Swift-BAT GRB Position','sample':'yes'},
    'bat_pos_err': {'definition':'Uncertainty in BAT position (arcmin)','type':'double','source':'Swift-BAT GRB Position','sample':3.0},
    'bat_image_signif': {'definition':'Significane of BAT Image Trigger (sigma)','type':'double','source':'Swift-BAT GRB Position','sample':7.19},
    'bat_ra': {'definition':'RA of burst, as determined by Swift BAT','type':'double','source':'Swift-BAT GRB Position','sample':31.835000000000001},
    'bat_rate_signif': {'definition':'Significance of BAT Rate Trigger (sigma)','type':'double','source':'Swift-BAT GRB Position','sample':21.02},
    'bat_trig_ind': {'definition':'BAT Trigger Index?','type':'double','source':'Swift-BAT GRB Position','sample':147.0},
    'bat_trig_ind_range': {'definition':'BAT Trigger Energy Range string (keV)','type':'string','source':'Swift-BAT GRB Position','sample':' 25-100 keV'},
    'bat_trigger_dur': {'definition':'Duration of BAT Trigger (seconds)','type':'double','source':'Swift-BAT GRB Position','sample':1.024},
    'burst_time_str': {'definition':'Burst time in HH:MM:SS format, read directly from Swift Catalog (less precise than grb_time_str)','type':'string','source':'SwiftCat','sample':'08:57:22'},
    'fluence': {'definition':'BAT fluence (15-150 keV) [10^-7 erg/cm^2]','type':'double','source':'SwiftCat','sample':28.0},
    'fluence_str': {'definition':'BAT Fluence String, read from Swift catalog','type':'string','source':'SwiftCat','sample':'28.00'},
    'grb_date_doy': {'definition':'Day of Year of GRB','type':'integer','source':'Swift-BAT GRB Position','sample':251},
    'grb_date_str': {'definition':'UT String Date in YY/MM/DD format','type':'string','source':'Swift-BAT GRB Position','sample':'06/09/08'},
    'grb_date_tjd': {'definition':'Truncated Julian Date of GRB','type':'integer','source':'Swift-BAT GRB Position','sample':13986},
    'grb_time_sod': {'definition':'Seconds of Day (SOD) of the Burst trigger','type':'string','source':'Swift-BAT GRB Position','sample':32242.34},
    'grb_time_str': {'definition':'UT Burst time string in HH:MM:SS.ss format','type':'string','source':'Swift-BAT GRB Position','sample':'08:57:22.34'},
    'm2_mag_isupper': {'definition':'Is the UVOT M2 Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'yes'},
    'moon_dist': {'definition':'Burst distance from Moon at time of XRT Position Determination','type':'double','source':'Swift-XRT Position','sample':37.619999999999997},
    'moon_illum': {'definition':'Percent illumination of the moon at time of XRT Position Determination','type':'double','Swift-XRT Position':'SwiftCat','sample':99.0},
    'out_dir': {'definition':'Directory where the HTML output for this event was created','type':'string','source':'SwiftCat','sample':'/Users/amorgan/Public/TestDir//228581'},
    'peakflux': {'definition':'BAT 1-sec peak photon flux (15-150 keV) [ph/cm^2/s]','type':'double','source':'SwiftCat','sample':3.0299999999999998},
    'peakflux_str': {'definition':'BAT 1-sec peak photon flux string','type':'string','source':'SwiftCat','sample':'3.03'},
    'reg_path': {'definition':'Path to the Region file created for this event','type':'string','source':'SwiftCat','sample':'/Users/amorgan/q_soft//store/sw228581.reg'},
    'sun_dist': {'definition':'Burst distance from Sun at time of XRT Position Determination','type':'double','':'Swift-XRT Position','sample':134.47999999999999},
    't90': {'definition':'Swift t90 value, parsed from t90_str','type':'double','source':'SwiftCat','sample':19.300000000000001},
    't90_str': {'definition':'BAT T90 string, read from Swift Catalog','type':'string','source':'SwiftCat','sample':'19.300'},
    'triggerid_str': {'definition':'Swift TriggerID string, read from Swift Catalog','type':'string','source':'SwiftCat','sample':'228581'},
    'u_mag_isupper': {'definition':'Is the UVOT U Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'uvot_dec': {'definition':'Declination of burst, as determined by Swift UVOT afterglow ','type':'double','source':'Swift-UVOT Position','sample':0.34200000000000003},
    'uvot_list': {'definition':'List of initial UVOT magnitudes and upper limits (all but V)','type':'string','source':'SwiftCat','sample':'B=18.41|U=17.01|UVW1=18.61|UVM2>18.57|UVW2>19.07|White=15.06'},
    'uvot_pos_err': {'definition':'Uncertainty in UVOT Position (arcsec)','type':'double','source':'Swift-UVOT Position','sample':0.40000000000000002},
    'uvot_ra': {'definition':'RA of burst, as determined by Swift UVOT afterglow ','type':'double','source':'Swift-UVOT Position','sample':31.8264},
    'v_mag_isupper': {'definition':'Is the UVOT V Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'v_mag_str': {'definition':'Initial UVOT V Magnitude String','type':'string','source':'SwiftCat','sample':'V=16.85'},
    'w1_mag_isupper': {'definition':'Is the UVOT W1 Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'w2_mag_isupper': {'definition':'Is the UVOT W2 Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'yes'},
    'wh_mag_isupper': {'definition':'Is the UVOT White Magnitude an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'xrt_amplifier': {'definition':'XRT Amplifier?','type':'integer','source':'Swift-XRT Position','sample':2},
    'xrt_column': {'definition':'Excess NH Column (10^21 cm^-2) inferred from XRT Data - from Swift Catalog','type':'double','source':'SwiftCat','sample':0.80000000000000004},
    'xrt_column_str': {'definition':'String of XRT Column density from Swift Cat','type':'string','source':'SwiftCat','sample':'0.80'},
    'xrt_dec': {'definition':'Declination of burst, as determined by Swift XRT afterglow - from GCN','type':'double','source':'Swift-XRT Position','sample':0.34160000000000001},
    'xrt_dec_str': {'definition':'XRT Dec String','type':'string','source':'SwiftCat','sample':'00:20:28.8'},
    'xrt_inten': {'definition':'XRT Flux as measured at the time of position determination','type':'double','source':'Swift-XRT Position','sample':4.9800000000000004e-10},
    'xrt_pos': {'definition':'XRT Position tuple (RA,dec) as determined by swift XRT afterglow; from Swift Catalog','type':'tuple','source':'SwiftCat','sample':(31.825875, 0.34133333333333332)},
    'xrt_pos_err': {'definition':'Uncertainty in XRT position (arcsec)','type':'double','source':'GCN','sample':6.0999999999999996},
    'xrt_ra': {'definition':'Right Ascention of burst, as determined by Swift XRT afterglow - from GCN','type':'double','source':'Swift-XRT Position','sample':31.825500000000002},
    'xrt_ra_str': {'definition':'XRT RA String in HH:MM:SS.ss format','type':'string','source':'SwiftCat','sample':'02:07:18.21'},
    'xrt_signif': {'definition':'Significance of the XRT detection (sigma)','type':'double','source':'Swift-XRT Position','sample':5.0},
    'xrt_tam0': {'definition':'XRT TAM 0?','type':'double','source':'Swift-XRT Position','sample':327.63},
    'xrt_tam1': {'definition':'XRT TAM 1?','type':'double','source':'Swift-XRT Position','sample':237.21000000000001},
    'xrt_tam2': {'definition':'XRT TAM 2?','type':'double','source':'Swift-XRT Position','sample':261.25999999999999},
    'xrt_tam3': {'definition':'XRT TAM 3?','type':'double','source':'Swift-XRT Position','sample':243.41999999999999},
    'xrt_waveform': {'definition':'XRT Waveform?','type':'integer','source':'SwiftCat','sample':134},
    'z': {'definition':'Redshift','type':'double','source':'SwiftCat','sample':2.4300000000000002},
    'z_class': {'definition':'Redshift class (high_z,medium_z,low_z) derived from actual redshift value','type':'string','source':'Derived','sample':'medium_z'},
    'z_isupper': {'definition':'Was the redshift value reported an upper limit?','type':'string','source':'SwiftCat','sample':'no'},
    'z_str': {'definition':'Redshift String','type':'string','source':'SwiftCat','sample':'2.43 (Gemini-North: absorption)'}}

    if keyword in helpdict:
        return helpdict[keyword]
    else: 
        print 'Keyword unknown.  Check keyword or tell Adam to update his dictionary'

def LoadDB(name, clobber=False):
    ### LOAD or CREATE PICKLE STORAGE FILE 
    # Attempt to load pickle file
    pklpath = storepath+'DB_'+str(name)+'.pkl'
    loadeddb = qPickle.load(pklpath)
    # If couldn't load, or clobber == True, create a new instance of the class
    if clobber or not loadeddb:
        # Create new instance of db Notice
        loadeddb = GRBdb(name)
        try:
            # Extract values from GCN Notice
            # loadeddb.extract_values()
            # loadeddb.get_positions()
            if loadeddb.successful_load:
                # Save new Pickle file
                qPickle.save(loadeddb,pklpath,clobber=True)
            else:
                print 'Could not succesfully load db.'
                return
        except:
            print "Could not Extract Values for db."
    return loadeddb

def SaveDB(loadeddb):
    # If the attributes of the loaded GCN have changed since loading, use 
    # this to save the new version of the GCN.
    pklpath = storepath+'DB_'+str(loaddeddb.name)+'.pkl'
    qPickle.save(loaddeddb,pklpath,clobber=True)

class GRBdb:
    '''Instance of a grb database'''
    def __init__(self,name):
        
        self.date_created = time.ctime()
        try: 
            self.dict = self.collect()
            self.name = str(name)
            if len(self.name) > 20:
                raise ValueError('Length of db name too long')
            self.successful_load = True
        except:
            self.successful_load = False
        
        
    def collect(self,incl_nat=True,incl_fc=False,incl_reg=True,make_html=True,\
                html_path='/Users/amorgan/Public/TestDir/'):
        '''A wrapper around all the parsers written to collect GRB information
        into a single dictionary, with GRB phone numbers as the keys.  Each key 
        has several attributes, and as of 01/10/10 the GRBs with the most
        attributes (131, using incl_nat=True, incl_fc=False, incl_reg=True,
        make_html=True (at least 3 more attributes, Beta and two limits, also exist
        but may not be in the 131.  At least 134 total attributes.)) are as follows: 
    
    
        080413B
        060908 - new redshift (1.8836) from Fynbo et al- redo Nat Values!
        080916A
        081008
        090618
        081222
        080319C
        090424
        
        Nat's Norm@15keV (Bat Spectra) not in catalog?
        Energy Fluence (15-150 keV) [erg/cm^2]
        Photon Fluence (15-150 keV) [ph/cm^2]	
        Peak Energy Flux (15-350 keV) [erg/cm^2/s]	
        Peak Energy Flux (15-150 keV) [erg/cm^2/s]	
        Peak Photon Flux (15-150 keV) [ph/cm^2/s]	
        Peak Photon Flux (50-300 keV) [ph/cm^2/s]	
        Eiso [erg]	
    
        rT_0.90	 
        rT_0.50
        T_av	 
        T_max	 
        T_rise	 
        T_fall	 
        Cts	 
        Rate_pk	 (Though Rate_pk/Cts is present)
        Band
        '''
        print '\nNow loading Swift Online Catalog Entries'
        swiftcatdict = ParseSwiftCat.parseswiftcat(storepath+'grb_table_current.txt')
        if incl_nat:
            from RedshiftMachine import ParseNatCat
            print "Now loading Nat's Catalog"
            natcatdict = ParseNatCat.load_natcats(def_bat_natcats,def_xrt_natcats)
        collected_dict = {}
        failed_gcn_grbs = []
        failed_nat_grbs = []
        failed_finding_charts = []
        for grb_str,catdict in swiftcatdict.iteritems():
            # If the swiftcat has a TRIGGERID associated with it, grab the trigger id
            # For now, only collect if it has an associated (redshift)/ Triggerid
            if 'triggerid_str' in catdict:
                fc_path=None
                reg_path=None
            
                trigid_str = catdict['triggerid_str']
                print '\nNow collecting GCN entries for trigger %s, GRB %s' % (trigid_str, grb_str)
                try:
                    triggerid=int(trigid_str)
                    loaded_gcn = LoadGCN.LoadGCN(triggerid)
                    loaded_gcn.extract_values()
                    source_name = 'Swift_%s-GRB%s' % (trigid_str, grb_str)
                    if incl_reg:
                        try:
                            reg_path = Signal._incl_reg(loaded_gcn)
                        except:
                            print 'cannot load reg'
                    if incl_fc:
                        try:
                            fc_path = Signal._incl_fc(loaded_gcn,src_name=source_name)
                        except:
                            failed_finding_charts.append(source_name)
                    if make_html:
                        try:
                            html_inst = Signal.make_grb_html(loaded_gcn, html_path=html_path, reg_path=reg_path, fc_path=fc_path)
                            web_dict =  object2dict(html_inst,include=['out_dir','fc_path','reg_path'])
                            catdict.update(web_dict)
                        except:
                            print 'cannot make html'

                    catdict.update(loaded_gcn.pdict)
                
                
                    # if make_finding_charts:
                    #     try:
                    #         from AutoRedux import qImage
                    #         source_name = 'GRB' + grb_str
                    #         qImage.MakeFindingChart(ra=loaded_gcn.pdict['xrt_ra'],dec=loaded_gcn.pdict['xrt_dec'],uncertainty=loaded_gcn.pdict['xrt_pos_err'],src_name=source_name,pos_label='XRT',survey='dss2red')
                    #     except: 
                    #         failed_finding_charts.append('GRB'+grb_str+' ('+trigid_str+')')
                except:
                    print "Cannot load GCN for trigger %s for GRB %s" % (trigid_str,grb_str)
                    failed_gcn_grbs.append('GRB'+grb_str+' ('+trigid_str+')')
            
                print "Now collecting Nat's Catalog entries for GRB %s" % (grb_str)
                GRBgrb_str = 'GRB'+grb_str
                try:
                    catdict.update(natcatdict[GRBgrb_str])
                except:
                    try:
                        # Try stripping the trailing A off the name and see if its in nat's cat
                        catdict.update(natcatdict[GRBgrb_str.strip('A')])
                    except:
                        try:
                            # Try ADDING the trailing A onto the name and see if it's in natcat
                            catdict.update(natcatdict[GRBgrb_str+'A'])
                        except:
                            print "Cannot load Nat's entries for GRB %s" % (grb_str)
                            failed_nat_grbs.append(GRBgrb_str)
            
                subdict = {grb_str:catdict}
                collected_dict.update(subdict)
        
        self.failed_gcn_grbs = failed_gcn_grbs
        self.failed_nat_grbs = failed_nat_grbs
        self.failed_finding_charts = failed_finding_charts
        self.length = len(collected_dict)
        # WOULD BE NICE IF THE 'COLLECTED DICT' Was an OBJECT so these could be 
        # ATTRIBUTES.  and we could have the functions be attributes too.  Bah.
        # this would be nice.  So I made it so.  02/06/10
            
        print ''
        print len(collected_dict), ' entries in the collected dictionary'
        print 'GRBs failed to gather from GCN: ', failed_gcn_grbs   
        print "GRBs failed to gather from Nat's Catalogue: ", failed_nat_grbs  
        print "GRBs failed to obtain finding charts: ", failed_finding_charts           
        return collected_dict

    def update_class(self):
        '''
        Given conditions, assign and update the redshift class
        '''
        for grb in self.dict:
            if 'z' in self.dict[grb]:
                if self.dict[grb]['z'] > 5.0:
                    self.dict[grb]['z_class'] = 'high_z'
                else:
                    self.dict[grb]['z_class'] = 'low_z'
                
            if 'RT45' in self.dict[grb] and 'FL' in self.dict[grb]:
                flosqrtrt45 = self.dict[grb]['FL']*(self.dict[grb]['RT45']**-0.5)
                self.dict[grb]['FLoverSQRTRT'] = flosqrtrt45  # this should be a measure of S/N
        
            if 'v_mag_isupper' in self.dict[grb]:
                vmagbinarystr = self.dict[grb]['v_mag_isupper']
                if vmagbinarystr == 'no':
                    self.dict[grb]['v_mag_isupper_binary'] = 0
                else:
                    self.dict[grb]['v_mag_isupper_binary'] = 1
        
    def log_update_class(self,keylist):
        for grb in self.dict:
            for key in keylist:
                if key in self.dict[grb]:
                    newname = 'log_'+key
                    try:
                        self.dict[grb][newname] = numpy.log(self.dict[grb][key])
                    except:
                        print 'Cannot take the log of GRB %s %s' % (grb,key)
                    
        # Note that this will make potentially interesting negative numbers into NAN
        # Some better way to deal with negative numbers before taking their log maybe? 
      
      
    def norm_update_class(self,keylist):
        pass
        # 0) Get list of values only for Z?
        # 1) Make array
        # 2) Determine and remove outliers?
        # 3) Normalize (see below)
        # 4) Create new dictionary items? 
        #
        # arr = scipy.array(mylist)
        # return (arr-arr.mean())/(arr.std())
      
    def compare_z(self):
        for i in iter(self.dict):
            try:
                print '\nGRB',i,'sw:',self.dict[i]['z'],'nat:',self.dict[i]['Z']
            except:
                pass


    def ret_list(self,x,y,z=[]):
        '''Returns a list of all the values and index names for a specific key'''
        xlist = []
        ylist = []
        zlist = []
        ilist = []
    
        remove_short = True
    
        if z:
            for i in iter(self.dict):
                if x in self.dict[i] and y in self.dict[i] and z in self.dict[i]:
                    try:
                        # REMOVE SHORT BURSTS FROM LIST.  PUT THIS ELSEWHERE
                        if remove_short and self.dict[i]['T90'] < 2.0:
                            print 'removing %s with T90 %f' % (i, self.dict[i]['T90'])
                            continue
                        x_val=float(self.dict[i][x])
                        y_val=float(self.dict[i][y])
                        z_val=float(self.dict[i][z])
                        xlist.append(x_val)
                        ylist.append(y_val)  
                        zlist.append(z_val)
                        ilist.append(i)
                    except:
                        pass
            return((xlist,ylist,ilist,zlist))
        else:
            for i in iter(self.dict):
                if x in self.dict[i] and y in self.dict[i]:
                    try:
                        # REMOVE SHORT BURSTS FROM LIST.  PUT THIS ELSEWHERE
                        if remove_short and self.dict[i]['T90'] < 2.0:
                            print 'removing %s with T90 %f' % (i, self.dict[i]['T90'])
                            continue
                        x_val=float(self.dict[i][x])
                        y_val=float(self.dict[i][y])
                        xlist.append(x_val)
                        ylist.append(y_val)  
                        ilist.append(i)
                    except:
                        pass
            return((xlist,ylist,ilist))
    
    
    def grbplot(self,x,y,logx=False,logy=False):
        list_tup = ret_list(self.dict,x,y)
        xlist = list_tup[0]
        ylist = list_tup[1]
        if not logx and not logy:
            pylab.plot(xlist,ylist,'ro')
        if logx and not logy:
            pylab.semilogx(xlist,ylist,'ro')
        if logy and not logx:
            pylab.semilogy(xlist,ylist,'ro')
        if logy and logx:
            pylab.loglog(xlist,ylist,'ro')
        pylab.ylabel(y)
        pylab.xlabel(x)


    def grbannotesubplot(self,y_keys=['NH_PC','NH_PC','NH_WT','NH_WT'],x_keys=['NH_PC','NH_WT','NH_PC','NH_WT'],z_keys=['Z','Z','Z','Z'],logx=True,logy=True):
        '''Bah
    
        '''
        from Plotting.Annote import AnnotatedSubPlot
    
        remove_short = True
    
        if len(x_keys) != len(y_keys):
            raise(ValueError('Len of key lists do not match'))
    
        ind = 0
        xlist=[]
        ylist=[]
        zlist=[]
        annotelist=[]
    
        while ind < len(x_keys):
            x=x_keys[ind]
            y=y_keys[ind]
            if z_keys:
                z=z_keys[ind]
                list_tup = ret_list(self.dict,x,y,z)
            else:
                list_tup = ret_list(self.dict,x,y)
        
            currentgrbname = list_tup[2]
            annotelist.append(currentgrbname)
        
        
            xlist.append(list_tup[0])
            ylist.append(list_tup[1])
            if z_keys:
                zlist.append(list_tup[3])
            ind += 1 
        
        # Maybe instead of plotting a loglog plot, create a new parameter which
        # is the log of the value and plot a linear scatter.  This will make 
        # the plotting faster and allow for more flexibility 
        AnnotatedSubPlot(xlist,ylist,annotelist,zlist=zlist,ynames=y_keys,xnames=x_keys,znames=z_keys,logx=False,logy=False)

    def plotallvall(self,keylist,zval=None,remove_redundant=True,single_save=True):
        '''Plot all listed keywords against each other.  If zval is specified, use
        it as the third 'color' dimension.  If remove_redundant, do not plot e.g.
        a vs b AND b vs a (nor a vs a). Single_save makes one plot for each and 
        saves them as a text file. 
    
        '''
        xkeys=[]
        ykeys=[]
        zkeys=[]
        #convert keylist [a,b,c] into xkeys = [a,b,c,a,b,c,a,b,c], ykeys = [a,a,a,b,b,b,c,c,c]
        for key in keylist:
            for yek in keylist:
                ykeys.append(key)
                xkeys.append(yek)
                if zval:
                    zkeys.append(zval)
        if remove_redundant:
            # only save the plots that haven't been plotted already, or that arent X vs X.  
            # the following loop determines which items to strike from the list of keys. 
            # e.g. 
            # [a,b,c,a,b,c,a,b,c]
            # [a,a,a,b,b,b,c,c,c]
            # loop returns 0,3,4,6,7,8 so the lists would be reduced to 
            # [b,c,c]
            # [a,a,b]
            # set redundant values to -1
            q=len(keylist)
            i=0
            while i < q:
                j=0
                while j < i+1:
                    val= j + i*q 
                    print val
                    xkeys[val] = -1
                    ykeys[val] = -1
                    if zval:
                        zkeys[val] = -1
                    j+=1
                i+=1
            # remove redundant values - nice way to remove all values of a certain type from a list
            xkeys=filter(lambda aa: aa != -1, xkeys)
            ykeys=filter(lambda aa: aa != -1, ykeys)
            zkeys=filter(lambda aa: aa != -1, zkeys)
    
        if single_save:
            ind = 0
            while ind < len(xkeys):
                xkey = xkeys[ind]
                ykey = ykeys[ind]
                if zval:
                    zkey = zkeys[ind]
                    z_keys = [zkey]
                else:
                    z_keys=None
                grbannotesubplot(self.dict,x_keys=[xkey],y_keys=[ykey],z_keys=z_keys)
                figname = self.name
                figname += '_'+xkey+'_vs_'+ykey
                if z_keys:
                    figname += '_vs_' + zkey
                figname += '.png'
                figoutdir = storepath + 'figures/'+figname
                pylab.savefig(figoutdir)
                pylab.close()
                ind += 1 
        else:
            grbannotesubplot(self.dict,x_keys=xkeys,y_keys=ykeys,z_keys=zkeys)

    def testplotall(self):
        keylist = ['EP', 'EP0', 'FL', 'NH_PC', 'NH_WT', 'NH_PC_LATE', 'PK_O_CTS', 'T50', 'T90', 'RT45', 'Z', 'MAX_SNR', 'DT_MAX_SNR']
        self.plotallvall(keylist=keylist,zval='Z')
        #['n','n','y','n','n','n','n']
        pass
    
    def printall(self,keywordlist,suppress = True):
        '''Print all the keywords in the keyword list in the collected dictionary
        If suppress, then only print them out if all keywords are present for a 
        particular GRB
        '''
        print keywordlist
        for i in iter(self.dict):
            printstr = 'GRB' + i + ' '
            doprint = True
            for keyword in keywordlist:
                try:
                    printstr += str(self.dict[i][keyword])
                    printstr += ' '
                except:
                    printstr += '??? '
                    if suppress: doprint = False
            if doprint: print printstr
    

        
def createarff(outdict,keylist=['T90','FL','peakflux','NH_PC_LATE','wh_mag_isupper','v_mag_isupper'],\
                    attributeclass='z_class',classlist=['high_z','low_z']):
    # BAT Specific: T90 Duration, Fluence, 1-sec Peak Photon Flux
    # XRT Specific: Location, Column Density (NH)
    # UVOT Specific: V Magnitude, Other Filter Magnitudes
    
    # Open file
    arffpath = storepath+'redshiftmachine.arff'
    f=open(arffpath,'w')
    
    # Create .arff header
    f.write('% 1. Title: Redshift Predictor for Swift GRBs\n')
    f.write('% \n')
    f.write('% 2. Sources:\n')
    f.write('%     (a) Creator: Adam N. Morgan\n')
    f.write('%     (b) Data From: http://swift.gsfc.nasa.gov/docs/swift/archive/grb_table.html/\n')
    f.write('%     (c) Date: '+time.asctime()+'\n')
    f.write('% \n')
    f.write('% 3. This file was created automatically. \n')
    f.write('%    CHECK THE ATTRIBUTES before running Weka. \n')
    f.write('% \n')
    f.write('@RELATION swift_redshift\n')
    f.write('\n')
    
    # Create .arff attributes section 
    for keyitem in keylist:
        # If the key for the first item in the dictonary is not a string, assume it is a numeric quantity
        if type(outdict[outdict.keys()[0]][keyitem]).__name__ != 'str':
            keystring = ('@ATTRIBUTE %s NUMERIC\n') % keyitem
        else:
            # WARNING: MIGHT NOT BE YES OR NO - MORE OPTIONS COULD BE PRESENT
            f.write('% !CHECK ME:\n')
            keystring = ('@ATTRIBUTE %s {yes, no}\n') % keyitem
        f.write(keystring)
    classsubstr = ''
    for classitem in classlist:
        classsubstr += classitem
        if len(classlist) - classlist.index(classitem) != 1:
            classsubstr += ', ' 
    classstring = ('@ATTRIBUTE class {%s}\n') % classsubstr
    f.write(classstring)
    
    # Create .arff data section
    
    f.write('\n')
    f.write('@DATA\n')
    for entry in outdict.keys():
        # Output each entry according that appears in keylist.  If it doesn't
        # appear, output a single '?' as required by the .arff standard
        datastring = ''
        for keyitem in keylist:
            if outdict[entry].has_key(keyitem):
                datastring += str(outdict[entry][keyitem])
            else:
                datastring += '?'
            datastring += ','
        
        #only write the line if the attribute class exists!
        if attributeclass in outdict[entry]:
            datastring += outdict[entry][attributeclass]
            datastring += '\n'
            f.write(datastring)
    
    f.close()

if __name__ == '__main__':
    collect()
    sys.exit(0)     