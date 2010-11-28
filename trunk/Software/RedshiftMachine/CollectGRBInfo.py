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
from MiscBin.q import RemoveNaN
from MiscBin.q import object2dict
from MiscBin.q import where
from MiscBin.q import dec2sex
from MiscBin import qPickle
import pylab
import scipy
import numpy
from Plotting.ColorScatter import ColorScatter
from Plotting.Annote import AnnotatedSubPlot
from Phot import extinction


#from MiscBin.q import Standardize

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def_bat_natcats = [loadpath+'bat_catalog_apj.671.656.fits',\
                    loadpath+'bat_catalog_apj.711.495.fits',\
                    loadpath+'bat_catalog_100706.fits']
def_xrt_natcats = [loadpath+'xrt_catalog_100706.fits']

def whatis(keyword):
    '''Look up what a particular keyword represents in the collected dictionary.
    
    example values in the helpdict is for '060908' except Beta, which is from 081118
    
    '''

    helpdict = {
    'A': {'definition':'Alpha: One of the 4 parameters in the Band-function fit to the BAT gamma-ray spectrum.  Alpha is the power-law index *before* the peak in the fit and is typically clustered around -1','type':'double','source':'NatBat Spectra','speed':'processed','sample':-0.77270000000000005},
    'A1': {'definition':'Alpha Lower limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':-1.0469999999999999},
    'A2': {'definition':'Alpha Upper limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':-0.48049999999999998},
    'B': {'definition':'Beta: One of the 4 parameters in the Band-function fit to the BAT gamma-ray spectrum.  Beta is the power-law index *after* the peak and has a typical value of around -3 (and is required to be less than -2 for the energy to remain finite)','type':'double','source':'NatBat Spectra','speed':'processed','sample':-2.1383000000000001},
    'B1': {'definition':'Beta Lower limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':-2.3252999999999999},
    'B2': {'definition':'Beta Upper limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':-1.9612000000000001},
    'CHI2': {'definition':'Chi-squared residuals for the fit to the BAT gamma-ray spectrum','type':'double','source':'NatBat Spectra','speed':'processed','sample':36.090000000000003},
    'CHI2_PC': {'definition':'Chi-squared residuals for the XRT Spectral Fit using PC (Photon-counting Mode) data - see T_PC for time range','type':'string','source':'NatXRT','speed':'processed','sample':'31.90/44'},
    'CHI2_PC_LATE': {'definition':'XRT Spectral Fit Chi-squared residuals from Late Time (>10ks) XRT PC (Photon Counting Mode) Data - see T_PC_LATE for time range','type':'string','source':'NatXRT','speed':'late_processed','sample':'17.49/17'},
    'CHI2_WT': {'definition':'XRT Spectral Fit Chi-squared residuals from WT (Windowed Timing Mode) XRT Data - Typically the earliest data; see T_WT for time range','type':'string','source':'NatXRT','speed':'processed','sample':'17.85/20'},
    'DPK_O_CTS': {'definition':'Peak Countrate over Counts Uncertainty (s^-1)','type':'string','source':'NatBat','speed':'processed','sample':0.0084310400000000008},
    'DRT45': {'definition':'rT_0.45 Uncertainty (s)','type':'string','source':'NatBat Timing','speed':'processed','sample':0.10100000000000001},
    'DT50': {'definition':'T50 Uncertainty (s)','type':'string','source':'NatBat Timing','speed':'processed','sample':0.22600000000000001},
    'DT90': {'definition':'T90 Uncertainty (s)','type':'string','source':'NatBat Timing','speed':'processed','sample':0.16900000000000001},
    'DT_MAX_SNR': {'definition':'Duration of the time window containing the maximum s/n detection (s)','type':'string','source':'NatBat','speed':'processed','sample':11.220000000000001},
    'EISO': {'definition':'Bayes Isotropic equivalent Energy [erg]  from integrated BAT Gamma-ray light curve - USES REDSHIFT. (Approximate Bolometric Fluence [erg/cm^2] in the 1-10^4 keV band if no Redshift)','type':'double','source':'NatBat Spectra','speed':'processed','sample':6.9099999999999996e+52},
    'EISO1': {'definition':'Bayes E_iso lower limit','type':'double','source':'NatBat','speed':'processed','sample':5.6600000000000002e+52},
    'EISO2': {'definition':'Bayes E_iso upper limit','type':'double','source':'NatBat','speed':'processed','sample':1.0900000000000001e+53},
    'EP': {'definition':'Bayes E_peak: One of the 4 parameters in the Band-function fit to the BAT gamma-ray spectrum.  E_Peak is the energy (frequency) at which most of the energy is emitted.','type':'double','source':'NatBat Spectra - found with a model fit to the spectrum, using bayesian techniques. More trustworthy than EP0.','speed':'processed','sample':161.47},
    'EP0': {'definition':'E_peak: One of the 4 parameters in the Band-function fit to the BAT gamma-ray spectrum.  E_Peak is the energy (frequency) at which most of the energy is emitted.','type':'double','source':'NatBat Spectra - found with a model fit to the spectrum, using frequentist techniques. Less trustworthy than EP.','speed':'processed','sample':145.96340000000001},
    'EP01': {'definition':'E_peak lower limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':112.1662},
    'EP02': {'definition':'E_peak upper limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':253.4169},
    'EP1': {'definition':'Bayes E_peak lower limit ','type':'double','source':'NatBat Spectra','speed':'processed','sample':122.63},
    'EP2': {'definition':'Bayes E_peak upper limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':334.51999999999998},
    'FL': {'definition':'Energy Fluence (15-350 keV) [erg/cm^2] - mathematically, the integral of the BAT gamma-ray light curve over time','type':'double','source':'NatBat Spectra','speed':'processed','sample':4.1281000000000003e-06},
    'FL1': {'definition':'Energy Fluence (15-350 keV) lower limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':3.7743999999999998e-06},
    'FL2': {'definition':'Energy Fluence (15-350 keV) upper limit','type':'double','source':'NatBat Spectra','speed':'processed','sample':4.5635999999999998e-06},
    'FLX_PC': {'definition':'Average Flux from XRT PC Data (mJy @ 1 keV) using XRT PC (Photon-counting Mode) data - see T_PC for time range','type':'double','source':'NatXRT','speed':'processed','sample':6.1600000000000007e-05},
    'FLX_PC1': {'definition':'Avg Flux from PC Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'processed','sample':5.5609999999999998e-05},
    'FLX_PC2': {'definition':'Avg Flux from PC Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'processed','sample':6.8150000000000003e-05},
    'FLX_PC_LATE': {'definition':'Average Flux from Late Time (>10ks) XRT PC (Photon Counting Mode) Data (mJy @ 1 keV) - see T_PC_LATE for time range','type':'double','source':'NatXRT','speed':'late_processed','sample':8.8419999999999994e-06},
    'FLX_PC_LATE1': {'definition':'Avg Flux from Late Time PC Data Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'late_processed','sample':5.5670000000000002e-06},
    'FLX_PC_LATE2': {'definition':'Avg Flux from Late Time PC Data Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'late_processed','sample':1.429e-05},
    'FLX_WT': {'definition':'Average Flux (mJy @ 1 keV) from XRT WT Data (Windowed Timing Mode) - Typically the earliest data; see T_WT for time range','type':'double','source':'NatXRT','speed':'processed','sample':0.0060889999999999998},
    'FLX_WT1': {'definition':'Avg Flux from WT Lower Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'processed','sample':0.0052449999999999997},
    'FLX_WT2': {'definition':'Avg Flux from WT Upper Limit (mJy @ 1 keV)','type':'double','source':'NatXRT','speed':'processed','sample':0.007038},
    'GAM_PC': {'definition':'Gamma inferred from XRT PC (Photon-counting Mode) data - see T_PC for time range','type':'double','source':'NatXRT','speed':'processed','sample':1.8776999999999999},
    'GAM_PC1': {'definition':'PC Gamma lower Limit','type':'double','source':'NatXRT','speed':'processed','sample':1.7541},
    'GAM_PC2': {'definition':'PC Gamma upper limit','type':'double','source':'NatXRT','speed':'processed','sample':2.0095999999999998},
    'GAM_PC_LATE': {'definition':'Gamma inferred from late time (>10ks) XRT PC (Photon Counting Mode) Data - see T_PC_LATE for time range','type':'double','source':'NatXRT','speed':'late_processed','sample':1.9380999999999999},
    'GAM_PC_LATE1': {'definition':'Late time PC Gamma lower Limit','type':'double','source':'NatXRT','speed':'late_processed','sample':1.2849999999999999},
    'GAM_PC_LATE2': {'definition':'Late time PC Gamma upper Limit','type':'double','source':'NatXRT','speed':'late_processed','sample':2.6629999999999998},
    'GAM_WT': {'definition':'Gamma inferred from XRT WT Data (Windowed Timing Mode) - Typically the earliest data; see T_WT for time range','type':'double','source':'NatXRT','speed':'processed','sample':2.3073999999999999},
    'GAM_WT1': {'definition':'WT Gamma lower Limit','type':'double','source':'NatXRT','speed':'processed','sample':2.0787},
    'GAM_WT2': {'definition':'WT Gamma upper Limit','type':'double','source':'NatXRT','speed':'processed','sample':2.5615999999999999},
    'MAX_SNR': {'definition':'Maximum BAT S/N Ratio','type':'double','source':'NatBat','speed':'processed','sample':47.399999999999999},
    'MODEL': {'definition':'BAT Spectral Fit Model Used for Frequentist approach: Either Simple Power Law Model (PLM), powerlaw times an exponetial cutoff (PLEXP), or smoothly connected broken powerlaw - the Band GRB Model (GRBM)','type':'string','source':'NatBat Spectra','speed':'processed','sample':'PLEXP'},
    'NH_GAL': {'definition':'Galactic Neutral Hydrogen Column at GRB Location (cm^-2)','type':'double','source':'NatXRT','speed':'processed','sample':0.023400000000000001},
    'NH_PC': {'definition':'Excess N_H Column inferred from XRT PC Data (10^22 cm^-2) using PC (Photon-counting Mode) data - see T_PC for time range','type':'double','source':'NatXRT','speed':'processed','sample':0.017999999999999999},
    'NH_PC1': {'definition':'Excess PC N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'processed','sample':-0.0080000000000000002},
    'NH_PC2': {'definition':'Excess PC N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'processed','sample':0.050000000000000003},
    'NH_PC_LATE': {'definition':'Excess N_H Column inferred from late time (>10ks) XRT PC (Photon Counting Mode) Data (10^22 cm^-2) - see T_PC_LATE for time range','type':'double','source':'NatXRT','speed':'late_processed','sample':0.057000000000000002},
    'NH_PC_LATE1': {'definition':'Excess late time PC N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'late_processed','sample':-0.128},
    'NH_PC_LATE2': {'definition':'Excess late time PC N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'late_processed','sample':0.27500000000000002},
    'NH_WT': {'definition':'Excess N_H Column (10^22 cm^-2) inferred from XRT WT Data (Windowed Timing Mode) - Typically the earliest data; see T_WT for time range','type':'double','source':'NatXRT Spectra','speed':'processed','sample':0.052999999999999999},
    'NH_WT1': {'definition':'Excess WT N_H Column Lower Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'processed','sample':0.0070000000000000001},
    'NH_WT2': {'definition':'Excess WT N_H Column Upper Limit (10^22 cm^-2)','type':'double','source':'NatXRT','speed':'processed','sample':0.106},
    'NISO': {'definition':'Isotropic equivalent Photons [ph] from integrated BAT Gamma-ray light curve - USES REDSHIFT.  (approximate bolometric photon fluence [ph/cm^-2] in observer frame 1-10^4 keV band if no redshift)','type':'double','source':'NatBat','speed':'processed','sample':3.4699999999999998e+59},
    'NISO1': {'definition':'NISO lower limit','type':'double','source':'NatBat','speed':'processed','sample':2.61e+59},
    'NISO2': {'definition':'NISO upper limit','type':'double','source':'NatBat','speed':'processed','sample':5.0600000000000004e+59},
    'NU': {'definition':'BAT Spectral Fit Degrees Of Freedom','type':'double','source':'NatBat Spectra','speed':'processed','sample':54.0},
    'OBS': {'definition':'Swift Trigger Number','type':'string','source':'NatBat','speed':'processed','sample':'00228581'},
    'PK_O_CTS': {'definition':'Ratio of the peak rate Rate_p (in a time bin of width 0.01 dt_(S/N)) over the total source counts (cts).  Used to approximately relate the burst fluences to peak fluxes.','type':'double','source':'NatBat','speed':'processed','sample':0.114758},
    'RT45': {'definition':'BAT rT0.45 (Defined in Reichart et al. 2001) - yet another way of defining the duration of the initial burst of gamma-rays','type':'double','source':'NatBat Timing','speed':'processed','sample':5.2800000000000002},
    'T0': {'definition':'Lower BAT Spectral Time Region','type':'string','source':'NatBat Spectra','speed':'processed','sample':'-9.615'},
    'T1': {'definition':'Upper BAT Spectral Time Region','type':'string','source':'NatBat Spectra','speed':'processed','sample':'13.815'},
    'T50': {'definition':'BAT T50 Difference between 75th and 25th percentile time of total counts above background level relative to start of burst interval','type':'double','source':'NatBat Timing','speed':'processed','sample':7.2599999999999998},
    'T90': {'definition':'BAT T90 Difference between 95th and 5th percentile time of total counts above background level relative to start of burst interval - loosely consistant with Swift T90; Highly dependent on burst start & stop times','type':'double','source':'NatBat Timing','speed':'processed','sample':18.48},
    'T_PC': {'definition':'XRT Photon Counting Mode Time Region (ks post burst)','type':'string','source':'NatXRT Spectra','speed':'processed','sample':'0.156101-1087.47'},
    'T_PC_LATE': {'definition':'XRT Late Photon Counting Mode Time Region (at least 10ks post burst)','type':'string','source':'NatXRT Spectra','speed':'late_processed','sample':'10.0000-1087.4741'},
    'T_WT': {'definition':'XRT WindowTiming Mode Region Time Region (ks post burst)','type':'string','source':'NatXRT Spectra','speed':'processed','sample':'0.0802510-5.72213'},
    'UT': {'definition':'UT Time of Burst','type':'string','source':'NatBat','speed':'processed','sample':'20060908_085722.340000'},
    'Z': {'definition':'Redshift (from Nat)','type':'double','source':'NatBat Spectra','speed':'na','sample':2.4300000000000002},
    'b_mag_isupper': {'definition':'Is the UVOT B Magnitude an upper limit?','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'bat_bkg_dur': {'definition':'GCN BAT Background duration (seconds)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':8.0},
    'bat_bkg_inten': {'definition':'GCN BAT Background intensity (counts)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':41524.0},
    'bat_dec': {'definition':'Declination of burst, as determined by Swift BAT','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':0.37},
    'bat_img_peak': {'definition':'GCN GRB BAT Image Peak (image_counts)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':240.0},
    'bat_inten': {'definition':'GCN GRB BAT Intensity (counts)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':6439.0},
    'bat_is_rate_trig': {'definition':'Was the event a Rate trigger (yes) or image trigger (no)','type':'string','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':'yes'},
    'bat_pos_err': {'definition':'Uncertainty in BAT position (arcmin)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':3.0},
    'bat_image_signif': {'definition':'Significane of BAT Image Trigger (sigma)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':7.19},
    'bat_ra': {'definition':'RA of burst, as determined by Swift BAT','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':31.835000000000001},
    'bat_rate_signif': {'definition':'Significance of BAT Rate Trigger (sigma)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':21.02},
    'bat_trig_ind': {'definition':'BAT Trigger Index?','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':147.0},
    'bat_trig_ind_range': {'definition':'BAT Trigger Energy Range string (keV)','type':'string','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':' 25-100 keV'},
    'bat_trigger_dur': {'definition':'Duration of BAT Trigger (seconds)','type':'double','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':1.024},
    'best_ra': {'definition':'Right Ascention of burst, as determined - from the most accurate GCN available','type':'double','source':'GCN Notices','speed':'nfi_prompt','sample':31.825500000000002},
    'best_dec': {'definition':'Declination of burst, as determined  - from the most accurate GCN available','type':'double','source':'GCN Notices','speed':'nfi_prompt','sample':0.34160000000000001},
    'best_pos_err': {'definition':'Uncertainty in Best position (arcsec)','type':'double','source':'GCN','speed':'nfi_prompt','sample':6.0999999999999996},
    'best_pos_type': {'definition':'Source of best position','type':'string','source':'GCN','speed':'nfi_prompt','sample':'XRT'},
    'burst_time_str': {'definition':'Burst time in HH:MM:SS format, read directly from Swift Catalog (less precise than grb_time_str)','type':'string','source':'SwiftCat','speed':'bat_prompt','sample':'08:57:22'},
    'fluence': {'definition':'BAT fluence (15-150 keV) [10^-7 erg/cm^2] (see FL for full definition) - trust less than FL?','type':'double','source':'SwiftCat','speed':'processed','sample':28.0},
    'fluence_str': {'definition':'BAT Fluence String, read from Swift catalog','type':'string','source':'SwiftCat','speed':'processed','sample':'28.00'},
    'grb_date_doy': {'definition':'Day of Year of GRB','type':'integer','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':251},
    'grb_date_str': {'definition':'UT String Date in YY/MM/DD format','type':'string','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':'06/09/08'},
    'grb_date_tjd': {'definition':'Truncated Julian Date of GRB','type':'integer','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':13986},
    'grb_time_sod': {'definition':'Seconds of Day (SOD) of the Burst trigger','type':'string','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':32242.34},
    'grb_time_str': {'definition':'UT Burst time string in HH:MM:SS.ss format','type':'string','source':'Swift-BAT GRB Position','speed':'bat_prompt','sample':'08:57:22.34'},
    'm2_mag_isupper': {'definition':'Is the UVOT M2 Magnitude an upper limit?','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'yes'},
    'moon_dist': {'definition':'Burst distance from Moon at time of XRT Position Determination','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':37.619999999999997},
    'moon_illum': {'definition':'Percent illumination of the moon at time of XRT Position Determination','type':'double','Swift-XRT Position':'SwiftCat','speed':'nfi_prompt','sample':99.0},
    'out_dir': {'definition':'Directory where the HTML output for this event was created','type':'string','source':'SwiftCat','speed':'na','sample':'/Users/amorgan/Public/TestDir//228581'},
    'peakflux': {'definition':'1-second integrated peak photon flux of the BAT spectrum (15-150 keV) [ph/cm^2/s]','type':'double','source':'SwiftCat','speed':'processed','sample':3.0299999999999998},
    'peakflux_str': {'definition':'BAT 1-sec peak photon flux string','type':'string','source':'SwiftCat','speed':'processed','sample':'3.03'},
    'reg_path': {'definition':'Path to the Region file created for this event','type':'string','source':'SwiftCat','speed':'na','sample':'/Users/amorgan/q_soft//store/sw228581.reg'},
    'sun_dist': {'definition':'Burst distance from Sun at time of XRT Position Determination','type':'double','':'Swift-XRT Position','speed':'nfi_prompt','sample':134.47999999999999},
    't90': {'definition':'Swift t90 value, parsed from t90_str (see T90 for full definition)','type':'double','source':'SwiftCat','speed':'processed','sample':19.300000000000001},
    't90_str': {'definition':'BAT T90 string, read from Swift Catalog','type':'string','source':'SwiftCat','speed':'processed','sample':'19.300'},
    'triggerid_str': {'definition':'Swift TriggerID string, read from Swift Catalog','type':'string','source':'SwiftCat','speed':'na','sample':'228581'},
    'u_mag_isupper': {'definition':'Is the UVOT U Magnitude an upper limit?','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'uvot_dec': {'definition':'Declination of burst, as determined by Swift UVOT afterglow ','type':'double','source':'Swift-UVOT Position','speed':'nfi_prompt','sample':0.34200000000000003},
    'uvot_list': {'definition':'List of initial UVOT magnitudes and upper limits (all but V)','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'B=18.41|U=17.01|UVW1=18.61|UVM2>18.57|UVW2>19.07|White=15.06'},
    'uvot_pos_err': {'definition':'Uncertainty in UVOT Position (arcsec)','type':'double','source':'Swift-UVOT Position','speed':'nfi_prompt','sample':0.40000000000000002},
    'uvot_ra': {'definition':'RA of burst, as determined by Swift UVOT afterglow ','type':'double','source':'Swift-UVOT Position','speed':'nfi_prompt','sample':31.8264},
    'uvot_time_delta': {'definition':'time since burst of first uvot observation','type':'double','source':'SwiftCat','speed':'nfi_prompt','sample':116.2},
    'v_mag_isupper': {'definition':'Is the UVOT V Magnitude an upper limit? - V is our second best chance of getting a detection; typically second on-sky (sometimes first) and the reddest of the filters','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'uvot_detection': {'definition':'Is there a detection in the UVOT? Combines information from wh_mag_isupper and v_mag_isupper','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'v_mag_str': {'definition':'Initial UVOT V Magnitude String','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'V=16.85'},
    'w1_mag_isupper': {'definition':'Is the UVOT W1 Magnitude an upper limit?','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'w2_mag_isupper': {'definition':'Is the UVOT W2 Magnitude an upper limit?','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'yes'},
    'wh_mag_isupper': {'definition':'Is the UVOT White Magnitude an upper limit? - White mode is our best chance of getting a detection, typically first on-sky (but not always taken, particularly in earlier observations), and widest filter','type':'string','source':'SwiftCat','speed':'nfi_prompt','sample':'no'},
    'xrt_amplifier': {'definition':'XRT Amplifier?','type':'integer','source':'Swift-XRT Position','speed':'nfi_prompt','sample':2},
    'xrt_column': {'definition':'Excess NH Column (10^21 cm^-2) inferred from XRT Data - from Swift Catalog','type':'double','source':'SwiftCat','speed':'processed','sample':0.80000000000000004},
    'xrt_column_str': {'definition':'String of XRT Column density from Swift Cat','type':'string','source':'SwiftCat','speed':'processed','sample':'0.80'},
    'xrt_dec': {'definition':'Declination of burst, as determined by Swift XRT afterglow - from GCN','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':0.34160000000000001},
    'xrt_dec_str': {'definition':'XRT Dec String','type':'string','source':'SwiftCat','speed':'na','sample':'00:20:28.8'},
    'xrt_inten': {'definition':'XRT Flux as measured at the time of position determination','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':4.9800000000000004e-10},
    'xrt_pos': {'definition':'XRT Position tuple (RA,dec) as determined by swift XRT afterglow; from Swift Catalog','type':'tuple','source':'SwiftCat','speed':'nfi_prompt','sample':(31.825875, 0.34133333333333332)},
    'xrt_pos_err': {'definition':'Uncertainty in XRT position (arcsec)','type':'double','source':'GCN','speed':'nfi_prompt','sample':6.0999999999999996},
    'xrt_ra': {'definition':'Right Ascention of burst, as determined by Swift XRT afterglow - from GCN','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':31.825500000000002},
    'xrt_ra_str': {'definition':'XRT RA String in HH:MM:SS.ss format','type':'string','source':'SwiftCat','speed':'na','sample':'02:07:18.21'},
    'xrt_signif': {'definition':'Significance of the XRT detection (sigma)','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':5.0},
    'xrt_tam0': {'definition':'XRT TAM 0?','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':327.63},
    'xrt_tam1': {'definition':'XRT TAM 1?','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':237.21000000000001},
    'xrt_tam2': {'definition':'XRT TAM 2?','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':261.25999999999999},
    'xrt_tam3': {'definition':'XRT TAM 3?','type':'double','source':'Swift-XRT Position','speed':'nfi_prompt','sample':243.41999999999999},
    'xrt_time_delta': {'definition':'time since burst of first xrt observation','type':'double','source':'SwiftCat','speed':'nfi_prompt','sample':86.2},
    'xrt_waveform': {'definition':'XRT Waveform?','type':'integer','source':'Swift-XRT Position','speed':'nfi_prompt','sample':134},
    'z': {'definition':'Redshift','type':'double','source':'SwiftCat','speed':'na','sample':2.4300000000000002},
    'z_class': {'definition':'Redshift class (high_z,medium_z,low_z) derived from actual redshift value','type':'string','source':'Derived','speed':'na','sample':'medium_z'},
    'z_isupper': {'definition':'Was the redshift value reported an upper limit?','type':'string','source':'SwiftCat','speed':'na','sample':'no'},
    'z_str': {'definition':'Redshift String','type':'string','source':'SwiftCat','speed':'na','sample':'2.43 (Gemini-North: absorption)'}}

    if keyword in helpdict:
        return helpdict[keyword]
    elif keyword == 'all':
        return helpdict
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
    pklpath = storepath+'DB_'+str(getattr(loadeddb,'name'))+'.pkl'
    qPickle.save(loadeddb,pklpath,clobber=True)

class GRBdb:
    '''Instance of a grb database'''
    def __init__(self,name,incl_nat=True,incl_fc=False,incl_reg=True,
                make_html=True,html_path='/Users/amorgan/Public/TestDir/'):
        
        self.date_created = time.ctime()
        try: 
            self.name = str(name)
            self.incl_fc = incl_fc
            self.incl_nat = incl_nat
            self.incl_reg = incl_reg
            self.make_html = make_html
            self.html_path = html_path
            self.class_updated = False
            self.dict = self.collect()
        
            if len(self.name) > 20:
                raise ValueError('Length of db name too long')
            self.successful_load = True
        except:
            self.successful_load = False
        
        
    def collect(self):
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
        
        incl_fc = self.incl_fc
        incl_nat = self.incl_nat
        incl_reg = self.incl_reg
        make_html = self.make_html
        html_path = self.html_path
        
        print '\nNow loading Swift Online Catalog Entries'
        swiftcatdict = ParseSwiftCat.parseswiftcat(loadpath+'grb_table_current.txt')
        if incl_nat:
            from RedshiftMachine import ParseNatCat
            print "Now loading Nat's Catalog"
            natcatdict = ParseNatCat.load_natcats(def_bat_natcats,def_xrt_natcats)
        collected_dict = {}
        failed_gcn_grbs = []
        failed_nat_grbs = []
        failed_gcn_grbs_ids = []  # includes just the grb_str e.g. '090313'
        failed_nat_grbs_ids = []  # includes just the grb_str e.g. '090313'
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
                    loaded_gcn.get_positions()
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
                    failed_gcn_grbs_ids.append(grb_str)
                    
            
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
                            failed_nat_grbs_ids.append(grb_str)
            
                subdict = {grb_str:catdict}
                collected_dict.update(subdict)
        
        failed_nat_grbs.sort()
        failed_gcn_grbs.sort()
        failed_finding_charts.sort()
        self.failed_gcn_grbs = failed_gcn_grbs
        self.failed_nat_grbs = failed_nat_grbs
        self.failed_gcn_grbs_ids = failed_gcn_grbs_ids
        self.failed_nat_grbs_ids = failed_nat_grbs_ids
        self.failed_finding_charts = failed_finding_charts
        self.length = len(collected_dict)
        # Can potentially update later to cull out the failed_gcn_grb_ids, etc
        
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
                
            if 'T90' in self.dict[grb] and 'FL' in self.dict[grb]:
                flosqrtrt90 = self.dict[grb]['FL']*(self.dict[grb]['T90']**-0.5)
                self.dict[grb]['FL_over_SQRT_T90'] = flosqrtrt90  # this should be a measure of S/N
            
            # Make new attribute checking if there's ANY uvot detection
            UVOT_Detection = 'no'
            if 'wh_mag_isupper' in self.dict[grb]:
                if self.dict[grb]['wh_mag_isupper'] == 'no':
                    UVOT_Detection = 'yes'
                
            elif 'v_mag_isupper' in self.dict[grb]:
                if self.dict[grb]['v_mag_isupper'] == 'no':
                    UVOT_Detection = 'yes'                
            self.dict[grb]['uvot_detection'] = UVOT_Detection
            
            # Make new attribute grabbing the extinction column.  
            # Might want to change the observation epoch.
            if 'best_ra' in self.dict[grb] and 'best_dec' in self.dict[grb]:
                best_position = dec2sex((self.dict[grb]['best_ra'],self.dict[grb]['best_dec']))
                ext_list = extinction.extinction(lon=best_position[0],\
                    lat=best_position[1],system_in='Equatorial',\
                    system_out='Galactic',obs_epoch="2005.0")
                gal_EB_V = ext_list[0]
            self.dict[grb]['gal_EB_V'] = gal_EB_V
            self.class_updated = True
            
      
    def MakeNomArr(self,key):
        '''Same as MakeAttrArr, but for nominal values (i.e. only create array and subarray;
        we can't calculate a mean or stdev for these values.) 
        '''
        arr = numpy.array(map(lambda x:x[key] if key in x else numpy.nan, self.dict.itervalues()))
        namearr = numpy.array(map(lambda x:x[0] if key in x[1] else numpy.nan, self.dict.iteritems()))
        subarr = arr[where(arr,'nan','!=')]
        keydict = {'array':arr,'names':namearr,'subarr':subarr,'type':'nominal'}
        setattr(self,key,keydict)
        
    def MakeBinArr(self,key,truval):
        zeros = numpy.zeros(self.length) # Create array of zeros
        nans = numpy.nonzero(getattr(self,key)['array'] == 'nan')
        inds = numpy.nonzero(getattr(self,key)['array'] == truval) #get indices
        namearr = getattr(self,key)['names']
        zeros[inds] = 1.0 #convert locations of indices to 1.0
        zeros[nans] = numpy.nan
        subarr = RemoveNaN(zeros)
        newkey = key + '_binary' #make new key name
        keydict = {'array':zeros,'names':namearr,'subarr':subarr,'type':'binary'}
        setattr(self,newkey,keydict)
    
    def MakeAttrArr(self,key,poserrkey=None,negerrkey=None,DeltaErr=True):
        '''For the key, create a numpy array of all the values
        Good for getting the mean, std dev, etc.  And for plotting!
        
        The errors expected are to be DeltaErr's, that is, distance above and below the mean.
        i.e value = 5.0, negerr = 0.2, poserr = 0.3
        Sometimes however they will be reported as High val, low val, medium val.
        i.e. value = 5.0, negerr = 4.8, poserr = 5.3
        In this case, set DeltaErr = False
        
        Leading to plotting with error bars: 
        In [428]: db.MakeAttrArr('T90',negerrkey='DT90',poserrkey='DT90')
        In [429]: db.MakeAttrArr('FL',negerrkey='FL1',poserrkey='FL2',DeltaErr=False)
        In [430]: pylab.errorbar(y=db.FL['array'],x=db.T90['array'],xerr=[db.T90['negerrarr'],db.T90['poserrarr']],yerr=[db.FL['negerrarr'],db.FL['poserrarr']],fmt='ro')
        
        
        Dictonary with keywords:
            array - full array for all the entries - use for plotting!
            poserrarr - positive error in the values of the array
            negerrarr - negative error in the values of the array
            subarray - only the non-NaN entries - do not use for indexing or plotting!
            mean - calculated from subarray
            median - calculated from subarray
            std - calculated from subarray
        '''
        arr = numpy.array(map(lambda x:x[key] if key in x else numpy.nan, self.dict.itervalues()))
        namearr = numpy.array(map(lambda x:x[0] if key in x[1] else numpy.nan, self.dict.iteritems()))
        subarr = RemoveNaN(arr)
        mean = subarr.mean()
        median = numpy.median(subarr)
        std = subarr.std()
        keydict = {'array':arr,'subarr':subarr,'names':namearr,'mean':mean,'median':median,'std':std,'type':'numeric'}
        if poserrkey:
            errarr = numpy.array(map(lambda x:x[poserrkey] if poserrkey in x else numpy.nan, self.dict.itervalues()))
            if not DeltaErr:
                errarr = errarr - arr
            keydict['poserrarr'] = errarr
        if negerrkey:
            errarr = numpy.array(map(lambda x:x[negerrkey] if negerrkey in x else numpy.nan, self.dict.itervalues()))
            if not DeltaErr:
                errarr = arr - errarr
            keydict['negerrarr'] = errarr
        # Errors expected are Figure out whether the     
        
        setattr(self,key,keydict)
   
    def log_update_class(self,keylist):
        '''Create offset if there are negative values before taking the logarithm?
        
        Note that this will make potentially interesting negative numbers into NAN
        Some better way to deal with negative numbers before taking their log maybe?
        '''
        for grb in self.dict:
            for key in keylist:
                if key in self.dict[grb]:
                    newname = 'log_'+key
                    try:
                        self.dict[grb][newname] = numpy.log(self.dict[grb][key])
                        self.MakeAttrArr(newname)
                    except:
                        print 'Cannot take the log of GRB %s %s' % (grb,key)
        
    
    def norm_update_class(self,keylist):
        '''
        # 0) Get list of values only for Z?
        # 1) Make array
        # 2) Determine and remove outliers?
        # 3) Normalize (see below)
        # 4) Create new dictionary items? 
        #
        # arr = scipy.array(mylist)
        # return (arr-arr.mean())/(arr.std())
        
        # In the code below I use the entire db to calc the mean/std.  Maybe only use
        # The values for which I have a redshift?
        
        '''
        for grb in self.dict:
            for key in keylist:
                if not hasattr(self,key):
                    self.MakeAttrArr(key)
                # Grab dict of value array
                tmp_dict = getattr(self,key)
                # Check to make sure there's the right # of entries in the temp arr
                if len(tmp_dict['array']) == self.length:
                    # THESE ARE IN NESTED FOR LOOPS  - Remove? - Did this.  Made function MakeAttrArr
                    if key in self.dict[grb]:
                        newname = 'norm_'+key
                        try:
                            self.dict[grb][newname] = (self.dict[grb][key]-tmp_dict['mean'])/tmp_dict['std']
                        except:
                            print 'Cannot normalize %s for GRB %s' % (key,grb)
                
                else:
                    print 'Not enough etries in array for key %s - is this the right keyword? Did the database change?' % (key)
        # Make an attribute array for every key we looped thru
        for key in keylist:
            newkey = 'norm_'+key
            self.MakeAttrArr(newkey)
                    
    def compare_z(self):
        for i in iter(self.dict):
            try:
                print '\nGRB',i,'sw:',self.dict[i]['z'],'nat:',self.dict[i]['Z']
            except:
                pass


    def ret_list(self,x,y,z=[]):
        '''Returns a list of all the values and index names for a specific key
        
        THIS SHOULD BE DEPRECIATED.  Use the arrays created by MakeAllAttr instead.
        '''
        
        print 'WARNING: ret_list is depreciated.'
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
    

    def grbplot(self,x_key,y_key,z_key=None,logx=False,logy=False,yjitter=0.0,\
        xjitter=0.0,discrete=0,marker='o'):
        '''Plot two keys against each other, with an optional third key as 
        the colorbar parameter.  Specify xjitter or yjitter to add a bit 
        of scatter to the plot for visualization reasons.  This replaces the
        old grbplot function which read from the dictionaries, and instead
        uses the arrays created by MakeAllAttr.
        '''
        # list_tup = self.ret_list(x,y)
        xlist = getattr(self,x_key)['array']
        ylist = getattr(self,y_key)['array']
        zlist = None
        if z_key:
            zlist = getattr(self,z_key)['array']
        if not logx and not logy:
            ColorScatter(xlist,ylist,zlist,yjitter=yjitter,xjitter=xjitter,\
                discrete=discrete,marker=marker)
        if logx and not logy:
            pylab.semilogx()
        if logy and not logx:
            pylab.semilogy()
        if logy and logx:
            pylab.loglog()
        pylab.ylabel(y_key)
        pylab.xlabel(x_key)


    def grbannotesubplot(self,\
        x_keys=['NH_PC','NH_PC','NH_WT','NH_WT'],\
        y_keys=['NH_PC','NH_WT','NH_PC','NH_WT'],\
        z_keys=['Z','Z','Z','Z'],\
        logx=False,logy=False):
        '''Create an annotated sub plot of the GRB parameters specified in the
        keys.
    
        '''
        
        
        xlist = [getattr(self,key)['array'] for key in x_keys]
        ylist = [getattr(self,key)['array'] for key in y_keys]
        annotelist = [getattr(self,key)['names'] for key in y_keys]
        zlist = [getattr(self,key)['array'] for key in z_keys]
        
        AnnotatedSubPlot(xlist,ylist,annotelist,zlist=zlist,ynames=y_keys,\
            xnames=x_keys,znames=z_keys,logx=logx,logy=logy)
                
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
                else:
                    zkey=None
                self.grbplot(x_key=xkey,y_key=ykey,z_key=zkey)
                figname = self.name
                figname += '_'+xkey+'_vs_'+ykey
                if zkey:
                    figname += '_vs_' + zkey
                figname += '.png'
                figoutdir = storepath + 'figures/'+figname
                pylab.savefig(figoutdir)
                pylab.close()
                ind += 1 
        else:
            self.grbannotesubplot(x_keys=xkeys,y_keys=ykeys,z_keys=zkeys)


    def test_log_update(self,plot=True,hist=True):
        # Remove the shorts before removing outliers so as to not bias the sample
        self.removeShort(remove_no_redshift=False)
        
        if not self.class_updated:
            self.update_class()
        self.MakeAllAttr()
        
        
        keys_to_log = ['gal_EB_V','uvot_time_delta','xrt_signif', 'bat_rate_signif', 'bat_image_signif', 'EP', 'EP0', 'FL', 'NH_PC', 'NH_WT', 'NH_PC_LATE', 'PK_O_CTS', 'T90', 'RT45', 'MAX_SNR', 'DT_MAX_SNR','peakflux','bat_inten','xrt_column','FL_over_SQRT_T90']
        self.log_update_class(keys_to_log)
        keys_to_norm = ['log_xrt_signif', 'log_bat_rate_signif', 'log_bat_image_signif','log_EP', 'log_EP0', 'log_FL', 'log_NH_PC', 'log_NH_WT', 'log_NH_PC_LATE', 'log_PK_O_CTS', 'log_T90', 'log_RT45', 'log_MAX_SNR', 'log_DT_MAX_SNR', 'log_peakflux', 'log_bat_inten', 'log_xrt_column','log_FL_over_SQRT_T90']
        self.norm_update_class(keys_to_norm)
        keys_to_plot = ['Z', 'log_uvot_time_delta','norm_log_xrt_signif', 'norm_log_bat_rate_signif', 'norm_log_bat_image_signif', 'norm_log_EP', 'norm_log_EP0', 'norm_log_FL', 'norm_log_NH_PC', 'norm_log_NH_WT', 'norm_log_NH_PC_LATE', 'norm_log_PK_O_CTS', 'norm_log_T90', 'norm_log_RT45', 'norm_log_MAX_SNR', 'norm_log_DT_MAX_SNR', 'norm_log_peakflux', 'norm_log_bat_inten', 'norm_log_xrt_column','norm_log_FL_over_SQRT_T90','v_mag_isupper_binary','wh_mag_isupper_binary','bat_is_rate_trig_binary']
        keys_to_hist = [ 'log_uvot_time_delta','norm_log_xrt_signif', 'norm_log_bat_rate_signif', 'norm_log_bat_image_signif', 'norm_log_EP', 'norm_log_EP0', 'norm_log_FL', 'norm_log_NH_PC', 'norm_log_NH_WT', 'norm_log_NH_PC_LATE', 'norm_log_PK_O_CTS', 'norm_log_T90', 'norm_log_RT45', 'norm_log_MAX_SNR', 'norm_log_DT_MAX_SNR', 'norm_log_peakflux', 'norm_log_bat_inten', 'norm_log_xrt_column','norm_log_FL_over_SQRT_T90']
        
        remove_outliers = True
        if remove_outliers:
            for attr in keys_to_log:
                self.removeOutliers(attr,threshold=0.32)
        
        if plot:
            self.plotallvall(keylist=keys_to_plot,zval='Z')
        if hist:
            self.DistHist(keylist=keys_to_hist)
    
    def MakeAllAttr(self):
        self.MakeAttrArr('A',negerrkey='A1',poserrkey='A2',DeltaErr=False)
        self.MakeAttrArr('B',negerrkey='B1',poserrkey='B2',DeltaErr=False)
        self.MakeAttrArr('CHI2')
        self.MakeNomArr('CHI2_PC')
        self.MakeNomArr('CHI2_PC_LATE')
        self.MakeNomArr('CHI2_WT')
        self.MakeAttrArr('DT_MAX_SNR')
        self.MakeAttrArr('EISO',negerrkey='EISO1',poserrkey='EISO2',DeltaErr=False)
        self.MakeAttrArr('EP',negerrkey='EP1',poserrkey='EP2',DeltaErr=False)
        self.MakeAttrArr('EP0',negerrkey='EP01',poserrkey='EP02',DeltaErr=False)
        self.MakeAttrArr('FL',negerrkey='FL1',poserrkey='FL2',DeltaErr=False)
        self.MakeAttrArr('FLX_PC',negerrkey='FLX_PC1',poserrkey='FLX_PC2',DeltaErr=False)
        self.MakeAttrArr('FLX_PC_LATE',negerrkey='FLX_PC_LATE1',poserrkey='FLX_PC_LATE2',DeltaErr=False)
        self.MakeAttrArr('FLX_WT',negerrkey='FLX_WT1',poserrkey='FLX_WT2',DeltaErr=False)
        self.MakeAttrArr('GAM_PC',negerrkey='GAM_PC1',poserrkey='GAM_PC2',DeltaErr=False)
        self.MakeAttrArr('GAM_PC_LATE',negerrkey='GAM_PC_LATE1',poserrkey='GAM_PC_LATE2',DeltaErr=False)
        self.MakeAttrArr('GAM_WT',negerrkey='GAM_WT1',poserrkey='GAM_WT2',DeltaErr=False)
        self.MakeAttrArr('MAX_SNR')
        self.MakeNomArr('MODEL')
        self.MakeAttrArr('NH_GAL')
        self.MakeAttrArr('NH_PC',negerrkey='NH_PC1',poserrkey='NH_PC2',DeltaErr=False)
        self.MakeAttrArr('NH_PC_LATE',negerrkey='NH_PC_LATE1',poserrkey='NH_PC_LATE2',DeltaErr=False)
        self.MakeAttrArr('NH_WT',negerrkey='NH_WT1',poserrkey='NH_WT2',DeltaErr=False)
        self.MakeAttrArr('NISO',negerrkey='NISO1',poserrkey='NISO2',DeltaErr=False)
        self.MakeAttrArr('NU')
        self.MakeNomArr('OBS')
        self.MakeAttrArr('PK_O_CTS',negerrkey='DPK_O_CTS',poserrkey='DPK_O_CTS',DeltaErr=True)
        self.MakeAttrArr('RT45',negerrkey='DRT45',poserrkey='DRT45',DeltaErr=True)
        self.MakeAttrArr('T50',negerrkey='DT50',poserrkey='DT50',DeltaErr=True)
        self.MakeAttrArr('T90',negerrkey='DT90',poserrkey='DT90',DeltaErr=True)
        self.MakeNomArr('T0')
        self.MakeNomArr('T1')
        self.MakeNomArr('UT')
        self.MakeNomArr('T_PC')
        self.MakeNomArr('T_PC_LATE')
        self.MakeNomArr('T_WT')
        self.MakeAttrArr('Z')
        
        # From the BAT Position GCN
        self.MakeAttrArr('bat_bkg_dur')
        self.MakeAttrArr('bat_bkg_inten')
        self.MakeAttrArr('bat_img_peak')
        self.MakeAttrArr('bat_inten')
        self.MakeNomArr('bat_is_rate_trig')
        self.MakeAttrArr('bat_image_signif')
        self.MakeAttrArr('bat_rate_signif')
        self.MakeAttrArr('bat_trig_ind')
        self.MakeNomArr('bat_trig_ind_range')
        self.MakeAttrArr('bat_trigger_dur')
        
        self.MakeAttrArr('grb_date_doy')
        self.MakeNomArr('grb_date_str')
        self.MakeAttrArr('grb_date_tjd')
        self.MakeAttrArr('grb_time_sod')
        self.MakeNomArr('grb_time_str')
        
        # From the XRT Position GCN
        self.MakeAttrArr('xrt_amplifier')
        self.MakeAttrArr('xrt_inten')
        self.MakeAttrArr('xrt_signif')
        self.MakeAttrArr('xrt_tam0')
        self.MakeAttrArr('xrt_tam1')
        self.MakeAttrArr('xrt_tam2')
        self.MakeAttrArr('xrt_tam3')
        self.MakeAttrArr('xrt_waveform')
        self.MakeAttrArr('moon_dist')
        self.MakeAttrArr('moon_illum')
        self.MakeAttrArr('sun_dist')
        
        # From SwiftCAT
        self.MakeAttrArr('fluence')
        self.MakeAttrArr('peakflux')
        self.MakeAttrArr('t90')
        self.MakeNomArr('v_mag_isupper')
        self.MakeNomArr('wh_mag_isupper')
        self.MakeAttrArr('xrt_column')
        self.MakeAttrArr('z')
        self.MakeNomArr('z_isupper')
        self.MakeNomArr('triggerid_str')
        self.MakeAttrArr('xrt_time_delta')
        self.MakeAttrArr('uvot_time_delta')
        
        # Self-created values
        if not self.class_updated:
            self.update_class()
        self.MakeAttrArr('FL_over_SQRT_T90')
        self.MakeNomArr('uvot_detection')
        self.MakeAttrArr('gal_EB_V')
        
        # Make the following Binary attributes
        keys_to_binary = ['v_mag_isupper','wh_mag_isupper','bat_is_rate_trig',
            'uvot_detection']
        for key in keys_to_binary:
            self.MakeBinArr(key,'yes')
        
    def DistHist(self,keylist):
        self.MakeAttrArr('T90',negerrkey='DT90',poserrkey='DT90')
        long_ind=numpy.nonzero(self.T90['array'] > 2.0) #indeces of long bursts
        remove_short = True
        for key in keylist:
            attr = getattr(self,key)
            if remove_short:
#                pylab.hist(attr['subarr'],bins=15)  # With Short Bursts
                testarr = attr['array'][long_ind] 
                pylab.hist(testarr[numpy.nonzero(numpy.nan_to_num(testarr))] ,bins=15)  # Without Short Bursts
            else:
                pylab.hist(attr['array'],bins=15)
            pylab.xlabel(key)

            figname = self.name
            figname += '_hist_'+key
            figname += '.png'
            figoutdir = storepath + 'figures/'+figname
            pylab.savefig(figoutdir)
            pylab.close()

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
        
    
    def removeOutliers(self,key,threshold=0.32):
        '''Given a key and error type, This will remove all values with 
        uncertainties larger than the threshold.  MUST Run self.MakeAllAttr() 
        first. '''
        
        if hasattr(self,key):
            
            keydict = getattr(self,key)
            if 'poserrarr' in keydict and 'negerrarr' in keydict:
            
                begin_num = len(RemoveNaN(keydict['array']))
                print "%s: Beginning with %i non-NaN values" % (key,begin_num)
            
                up_per_err = keydict['poserrarr']/keydict['array']
                lo_per_err = keydict['negerrarr']/keydict['array']
            
                #find indices of values where the percent error is larger than threshold
                up_ind = numpy.nonzero(up_per_err >= threshold)
                lo_ind = numpy.nonzero(lo_per_err >= threshold)
            
                # set indices = None, thus removing them
                keydict['array'][up_ind] = None
                keydict['array'][lo_ind] = None
                keydict['poserrarr'][up_ind] = None
                keydict['poserrarr'][lo_ind] = None
                keydict['negerrarr'][up_ind] = None
                keydict['negerrarr'][lo_ind] = None
                keydict['subarr'] = RemoveNaN(keydict['array'])
                
                # Now remove the outliers from any logged or normed arrays
                prefix_list = ['norm','log','norm_log']
                for prefix in prefix_list:
                    newkey = prefix + '_' + key
                    if hasattr(self,newkey):
                        print 'Also Removing for %s' % (newkey)
                        newkeydict = getattr(self,newkey)
                        # set indices = None, thus removing them
                        newkeydict['array'][up_ind] = None
                        newkeydict['array'][lo_ind] = None
                        try:
                            newkeydict['poserrarr'][up_ind] = None
                            newkeydict['poserrarr'][lo_ind] = None
                            newkeydict['negerrarr'][up_ind] = None
                            newkeydict['negerrarr'][lo_ind] = None
                        except(KeyError):
                            pass
                        newkeydict['subarr'] = RemoveNaN(keydict['array'])
                
            
                end_num = len(RemoveNaN(keydict['array']))
                print "Ending with %i non-NaN values" % (end_num)
                print "%i outliers removed." % (begin_num-end_num)
                print ''
            
                setattr(self,key,keydict)
            else:
                print 'Key %s does not have uncertainties; cannot remove outliers.' % (key)
        else:
            print 'Key %s does not exist' % (key)
        
        # In [21]: new=db.RT45['poserrarr']/db.RT45['array']
        # 
        # In [22]: numpy.nonzero(new>0.1)
        # Out[22]: 
        # (array([ 46,  54,  61,  71,  72,  95, 105, 127, 136, 142, 180, 260, 266,
        #        286, 293, 316, 320, 330, 334, 340, 343, 353, 359, 365, 381, 384,
        #        385, 400, 404, 441, 479, 512, 524, 530]),)
        # 
    
    
    def removeValues(self,key,argument):
        remove_list = []
        # Quick check to see if the argument is kind of formatted correctly
        allowed_arguments = ['>','<','=','!']
        arg_flag = 0
        for arg in allowed_arguments:
            if argument[0] == arg:
                arg_flag = 1
        if not arg_flag: 
            print 'Malformed argument; returning full array.'
            
        for ii in iter(self.dict):
            already_removed = False
            if key in self.dict[ii]:
                keyval = self.dict[ii][key]
            else:
                keyval = 'unknown'
            
            if not keyval == 'unknown':
                execstring = 'if keyval %s: remove_list.append(ii); already_removed=True '\
                    % (argument)
            else:
                execstring = 'pass'
            
            try:
                exec(execstring)
            except:
                print 'Cannot evaluate argument, returning'
                return 
        
        print 'initial length of dictionary: ' + str(len(self.dict))
        print ' removing ' + str(len(remove_list)) + ' GRBs'
        for key in remove_list:
            self.dict.pop(key)
        print 'new length of dictionary: ' + str(len(self.dict))
        #update the length
        self.length=len(self.dict)
        # Now remake the arrays
        self.MakeAllAttr()
        
        
    def removeShort(self,remove_no_redshift=True):
        '''Removes all short bursts with T90 < 2 and, if remove_no_redshift==
        True, those without a redshift'''
        remove_list = []
        for i in iter(self.dict):
            already_removed = False
            if 'T90' in self.dict[i]: 
                t90 = self.dict[i]['T90']
            elif 't90' in self.dict[i]:
                t90 = self.dict[i]['t90']
            else:
                t90 = 'unknown'
            if 'Z' not in self.dict[i]:
                z = None
            else:
                z = self.dict[i]['Z']
            if t90 < 2.0: 
                remove_list.append(i)
                already_removed=True
            if remove_no_redshift and not z and not already_removed:
                remove_list.append(i)
        print 'initial length of dictionary: ' + str(len(self.dict))
        print ' removing ' + str(len(remove_list)) + ' GRBs'
        for key in remove_list:
            self.dict.pop(key)
        print 'new length of dictionary: ' + str(len(self.dict))
        #update the length
        self.length=len(self.dict)
                    
                    
    def makeArffFromArray(self,time_list=['na','nfi_prompt', 'processed', 'bat_prompt', 'late_processed'],attrlist=['triggerid_str','Z','A','B','CHI2','CHI2_PC','CHI2_WT','CHI2_PC_LATE','DT_MAX_SNR','EP','EP0','FL','FLX_PC','FLX_PC_LATE','FLX_WT','GAM_PC','GAM_PC_LATE','GAM_WT','MAX_SNR','MODEL','NH_GAL','NH_PC','NH_PC_LATE','NH_WT','NU','PK_O_CTS','RT45','T50','T90','bat_bkg_inten','bat_image_signif','bat_img_peak','bat_inten','bat_is_rate_trig','bat_rate_signif','bat_trigger_dur','xrt_inten','xrt_signif','xrt_amplifier','xrt_tam0','xrt_tam1','xrt_tam2','xrt_tam3','fluence','peakflux','t90','v_mag_isupper','wh_mag_isupper','xrt_column'],inclerr=True):
        '''Create .arff file from array of attributes
        MUST Run self.MakeAllAttr() first.
        
        reduced_attr_list = ['A','B','EP0','FL','FLX_PC_LATE','GAM_PC','MAX_SNR','NH_PC','T90','bat_image_signif','bat_img_peak','bat_is_rate_trig','bat_trigger_dur','uvot_detection']
        
        '''
        
        # Open file
        arffpath = storepath+self.name+'.arff'
        arffpathpartial = arffpath + '_part'
        subpath = storepath+'redshiftdata'+'.txt'
        
        fmt = ''
        
        remove_outliers = True
        if remove_outliers:
            for attr in attrlist:
                self.removeOutliers(attr,threshold=0.32)
        
        helpdict = whatis('all')
        
        f=open(arffpathpartial,'w')

        # Create .arff header
        f.write('% 1. Title: Redshift Predictor for Swift GRBs\n')
        f.write('% \n')
        f.write('% 2. Sources:\n')
        f.write('%     (a) Creator: Adam N. Morgan\n')
        f.write('%     (b) Data From: GCN Notices (http://gcn.gsfc.nasa.gov/),\n')
        f.write('%                    http://astro.berkeley.edu/~nat/Swift/\n')
        f.write('%                    http://swift.gsfc.nasa.gov/docs/swift/archive/grb_table.html/\n')
        f.write('%     (c) Date: '+time.asctime()+'\n')
        f.write('% \n')
        f.write('% 3. This file was created automatically. \n')
        f.write('%    CHECK THE ATTRIBUTES before running Weka. \n')
        f.write('% \n')
        f.write('@RELATION swift_redshift\n')
        f.write('\n')
        
        fmt = ''
        nomkeystring = ''
        numkeystring = ''
        numattrlist = []
        nomattrlist = []
        
        # Create .arff attributes section 
        for keyitem in attrlist:
            
            # Check to see if the attribute exists
            if hasattr(self,keyitem):
                attrdict = getattr(self,keyitem)
                # Check if a numeric or nominal attribute
                if not 'type' in attrdict:
                    print 'type not in attribute dict.  Continuing...'
                    continue
                # Check if in the timedict.  Attributes made to binary values
                # Will have the same features, so just strip that word and 
                # look up the original in the helpdict
                if keyitem[-7:] == '_binary':
                    print keyitem
                    keyitem_lookup = keyitem.rstrip('binary').rstrip('_')
                    print keyitem_lookup
                else:
                    keyitem_lookup = keyitem
                helpdictitem = helpdict[keyitem_lookup]
                if not helpdictitem['speed'] in time_list:
                    print 'Speed for %s not in time_list; not including' % (keyitem)
                    continue
                if attrdict['type'] == 'numeric' or attrdict['type'] == 'binary':
                    numkeystring += ('@ATTRIBUTE %s NUMERIC\n') % keyitem
                    numattrlist.append(keyitem)
                    fmt += ',%f'
                    if inclerr and 'poserrarr' in attrdict and 'negerrarr' in attrdict:
                        posname = keyitem + '_poserr'
                        negname = keyitem + '_negerr'
                        numkeystring += ('@ATTRIBUTE %s NUMERIC\n') % posname
                        numkeystring += ('@ATTRIBUTE %s NUMERIC\n') % negname
                        fmt += ',%f,%f'
                elif attrdict['type'] == 'nominal':
                    # Do nominal Stuff
                    # WARNING: MIGHT NOT BE YES OR NO - MORE OPTIONS COULD BE PRESENT
                    f.write('% !CHECK ME:\n')
                    nomattrlist.append(keyitem)
                    fmt += ',%s'
                    nomkeystring += ('@ATTRIBUTE %s {yes, no}\n') % keyitem
                else:
                    print 'Attribute type is unknown (not nominal or numeric). Continuing..'
                # Write the keystring
            else:
                print 'Database does not have attribute %s' % (keyitem)

        f.write(numkeystring)
        f.write(nomkeystring)
            
        # Create .arff data section
        f.write('\n')
        f.write('@DATA\n')    
        
            
        if numattrlist:
            firstattr = numattrlist[0]
            numtotarr=numpy.array([getattr(self,firstattr)['array']])
            # Populate the total array
            for keyitem in numattrlist[1:]:
                attrdict = getattr(self,keyitem)
                numtotarr = numpy.concatenate((numtotarr,numpy.array([attrdict['array']])),axis=0)
                if inclerr and 'poserrarr' in attrdict and 'negerrarr' in attrdict:
                    numtotarr = numpy.concatenate((numtotarr,numpy.array([attrdict['poserrarr']])),axis=0)
                    numtotarr = numpy.concatenate((numtotarr,numpy.array([attrdict['negerrarr']])),axis=0)
            numarr2 = numtotarr
        else:
            numarr2 = numpy.array(['# NO NUMERIC ATTRIBUTES'])
            
        
        
        if nomattrlist:
            firstattr = nomattrlist[0]
            nomtotarr=numpy.array([getattr(self,firstattr)['array']])
            # Populate the total array
            for keyitem in nomattrlist[1:]:
                attrdict = getattr(self,keyitem)
                nomtotarr = numpy.concatenate((nomtotarr,numpy.array([attrdict['array']])),axis=0)
            nomarr2 = nomtotarr
        else:
            nomarr2 = numpy.array(['# NO NOMINAL ATTRIBUTES'])
        
        nomsubpath = subpath + 'nom'
        numsubpath = subpath + 'num'
        numnomsubpath = subpath + 'numnom'
        
        numpy.savetxt(numsubpath,numarr2,delimiter=',',fmt='%1.5e')
        numpy.savetxt(nomsubpath,nomarr2,delimiter=',',fmt='%s')
        
        cmd = 'cat %s %s > %s' %(numsubpath,nomsubpath,numnomsubpath)
        os.system(cmd)
        
        
        fixedarray = numpy.genfromtxt(numnomsubpath,delimiter=',',dtype='str')
        numpy.savetxt(subpath,fixedarray.T,delimiter=',',fmt='%s')
        
        
        f.close()
        
        # Find and replace all 'nan' strings with '?' recognized by weka
        stext = 'nan'
        rtext = '?'

        inputt = sys.stdin
        output = sys.stdout

        inputt = open(subpath)
        subpath2 = subpath + '_'
        output = open(subpath2, 'w')

        for s in inputt:
            output.write(s.replace(stext, rtext))

        inputt.close()
        output.close()
        
        # Combine the Header with the data
        cmd = 'cat %s %s > %s' %(arffpathpartial,subpath2,arffpath)
        os.system(cmd)
        


if __name__ == '__main__':
    collect()
    sys.exit(0)     