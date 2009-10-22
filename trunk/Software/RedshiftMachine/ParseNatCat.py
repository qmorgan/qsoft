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
If remove_zeroes == True, then remove the dictionary entries under the assumption 
that they were set to Zero in the catalog because the true value was unknown.

"""
import cPickle as pickle
import os
import sys
try: 
    import pyfits
except:
    'Unable to load pyfits module'
    sys.exit('Cannot Import PyFITS Module')

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

default_catlist = [storepath+'bat_catalog_07061275.fits',storepath+'bat_catalog_current.fits']

def read_nat_bat_cat(fitsname,remove_zeros=False):
    hdulist = pyfits.open(fitsname)
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    colnames = cols.names
    grblist = tbdata.field('GRB')
    nat_dict = {}
    for i in range(len(tbdata)):
        subdict = dict(zip(colnames,tbdata[i]))
        grbname = subdict.pop('GRB') # remove GRBname from dictionary
        # Had iteration problems trying to remove dictionary entries while
        # iterating over them, so just populate a new dictionary with values
        # that are NOT zero.
        if remove_zeros == True:
            subdict_zrem = {}
            for key,val in subdict.iteritems():
                if subdict[key] != 0.0:
                    subdict_zrem.update({key:val})
            subdict = subdict_zrem
        nat_dict.update({grbname:subdict})
    return nat_dict

def combine_natcats(natcatlist=default_catlist,remove_zeros=False):
    comb_nat_dict = {}
    for natcat in natcatlist:
        sub_nat_dict = read_nat_bat_cat(natcat,remove_zeros)
        comb_nat_dict.update(sub_nat_dict)
    # if remove_zeros = True:
    #     for i in iter(comb_nat_cat):
    #         for ii in iter(comb_nat_cat[i]):
    #             if comb
    return comb_nat_dict

