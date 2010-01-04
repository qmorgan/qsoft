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

default_bat_catlist = [storepath+'bat_catalog_07061275.fits',\
                    storepath+'bat_catalog_current.fits']

default_xrt_catlist = [storepath+'xrt_catalog_090831.fits']

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
        # that are NOT zero or 999.0 (which the n/a from the xrt cat are).
        if remove_zeros == True:
            subdict_zrem = {}
            for key,val in subdict.iteritems():
                if subdict[key] != 0.0 and subdict[key] != 999.0:
                    subdict_zrem.update({key:val})
            subdict = subdict_zrem
        nat_dict.update({grbname:subdict})
    return nat_dict

def combine_natcats(natcatlist=default_bat_catlist,remove_zeros=False):
    '''Use combine_netcats if one catalog goes up to a particular date and 
    another picks up where it left off - do NOT use it to combine, e.g., the
    XRT and the BAT catalogs.  Use to load EITHER the bat or xrt catalogs'''
    comb_nat_dict = {}
    for natcat in natcatlist:
        sub_nat_dict = read_nat_bat_cat(natcat,remove_zeros)
        comb_nat_dict.update(sub_nat_dict)
    # if remove_zeros = True:
    #     for i in iter(comb_nat_cat):
    #         for ii in iter(comb_nat_cat[i]):
    #             if comb
    return comb_nat_dict

def combine_natbatxrt(comb_bat_dict,comb_xrt_dict):
    '''Nat was nice enough to provide one catalog of xrt properties and
    another of the bat properties.  This function puts them both together
    after being parsed by the functions above.'''
    comb_nat_xrtbat_dict = {}
    # assume bat_dict > xrt_dict
    notxrtcount = 0
    notbatcount = 0
    # First grab the XRT values from the GRBs in the BAT cat, then vice versa
    for GRBgrb_str,catdict in comb_bat_dict.iteritems():
        try:
            catdict.update(comb_xrt_dict[GRBgrb_str])
        except:
            print '%s not in Nat XRT Dict' % (GRBgrb_str)
            notxrtcount += 1
        subdict = {GRBgrb_str:catdict}
        comb_nat_xrtbat_dict.update(subdict)
    for GRBgrb_str,catdict in comb_xrt_dict.iteritems():
        try:
            catdict.update(comb_bat_dict[GRBgrb_str])
        except:
            print '%s not in Nat BAT Dict' % (GRBgrb_str)
            notbatcount += 1
        subdict = {GRBgrb_str:catdict}
        comb_nat_xrtbat_dict.update(subdict)
    print "%i GRBs in Nat BAT catalog but not XRT" % (notxrtcount)
    print "%i GRBs in Nat XRT catalog but not BAT" % (notbatcount)
    return comb_nat_xrtbat_dict
    
def load_natcats(bat_catlist,xrt_catlist):
    '''Use to load BOTH the xrt and bat catalogs from Nat.'''
    xrtcat = combine_natcats(xrt_catlist,remove_zeros=True)
    print "length of xrt nat cat: ", len(xrtcat)
    batcat = combine_natcats(bat_catlist,remove_zeros=True)
    print "length of bat nat cat: ", len(batcat)
    batxrtcat = combine_natbatxrt(batcat,xrtcat)
    return batxrtcat