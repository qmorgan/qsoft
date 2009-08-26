#!/usr/bin/env python
# encoding: utf-8
"""
ParseNatCat.py
Author: Adam Morgan
Created: Aug 26, 2009
Last Updated: Aug 26, 2009
	Created
	
This program loads the fits tables from Nat's website http://astro.berkeley.edu/~nat/ 
as described in Butler et al. (2007; ApJ, 671, 656) and puts them in a dictionary.


"""
import cPickle as pickle
import os
import sys
try: 
    import pyfits
except:
    'Unable to load pyfits module'
    sys.exit['Cannot Import PyFITS Module']

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def read_nat_bat_cat(fitsname):
    hdulist = pyfits.open(fitsname)
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    

read_nat_bat_cat(storepath+'bat_catalog_07061275.fits')
read_nat_bat_cat(storepath+'bat_catalog_current.fits')