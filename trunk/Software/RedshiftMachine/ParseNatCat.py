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
import urllib2
import numpy 

try: 
    import pyfits
except:
    'Unable to load pyfits module'
    sys.exit('Cannot Import PyFITS Module')

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)

loadpath = os.environ.get("Q_DIR") + '/load/'
storepath = os.environ.get("Q_DIR") + '/store/'

default_bat_catlist = [loadpath+'bat_catalog_apj.671.656.fits',\
                    loadpath+'bat_catalog_apj.711.495.fits',\
                    loadpath+'bat_catalog_100706.fits']

default_xrt_catlist = [loadpath+'xrt_catalog_100706.fits']

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
    # Update to grab Nat's redshift probabilities
    batxrtcat = grab_nat_web_data(batxrtcat,filename='/bat/zprob.txt')
    return batxrtcat

def download_file(url,outfile):
    response = urllib2.urlopen(url)
    html = response.read()
    f = file(outfile,'w')
    f.write(html)
    f.close()

def grab_nat_web_data(natcatdict, filename='/bat/zprob.txt', clobber=False):
    base_url = 'http://astro.berkeley.edu/~nat/swift/'
    outpath = storepath + 'grbs/'
    for grbname in natcatdict.keys():
        fullurl = base_url + grbname + filename
        fullout = outpath + grbname + filename
        if not os.path.exists(os.path.dirname(fullout)):
            os.makedirs(os.path.dirname(fullout))
            
            # Try to download the file, if can't, skip to the next one
            try:
                download_file(fullurl,fullout)
            except:
                print 'Cannot download %s to %s' % (fullurl,fullout)
                continue 
        elif clobber == True:
            try:
                os.remove(fullout)
                download_file(fullurl,fullout)
            except:
                print 'Cannot download %s to %s' % (fullurl,fullout)
                continue
        else:
            pass
        # Grab probabilities if we're looking at that filename
        if filename == '/bat/zprob.txt':
            try:
                response = f = open(fullout,'r')
                html = f.readlines()
                z_pred_arr = []
                prob_arr = []
                ind = 0 # Grab the index to mark where redshifts are for summing probs
                for line in html:
                    if line[0] == '#':
                        continue
                    z_prob = line.rstrip('\n').split(' ')
                    z_pred_arr.append(z_prob[0])
                    prob_arr.append(float(z_prob[1]))
                    # What follows is very embarassing but I am lazy today
                    if line[0:3] == '1.0': ind1 = ind
                    if line[0:3] == '2.0': ind2 = ind
                    if line[0:3] == '3.0': ind3 = ind
                    if line[0:3] == '4.0': ind4 = ind
                    if line[0:3] == '5.0': ind5 = ind
                    ind += 1
                # Now sum up the probabilities to get the prob that lt some value
                norm_total = numpy.sum(numpy.array(prob_arr))
            
                prob_gt_1 = 1 - numpy.sum(numpy.array(prob_arr[0:ind1]))/norm_total
                prob_gt_2 = 1 - numpy.sum(numpy.array(prob_arr[0:ind2]))/norm_total
                prob_gt_3 = 1 - numpy.sum(numpy.array(prob_arr[0:ind3]))/norm_total
                prob_gt_4 = 1 - numpy.sum(numpy.array(prob_arr[0:ind4]))/norm_total
                prob_gt_5 = 1 - numpy.sum(numpy.array(prob_arr[0:ind5]))/norm_total
            
                natcatdict[grbname].update({'PROB_Z_GT_1':prob_gt_1})
                natcatdict[grbname].update({'PROB_Z_GT_2':prob_gt_2})
                natcatdict[grbname].update({'PROB_Z_GT_3':prob_gt_3})
                natcatdict[grbname].update({'PROB_Z_GT_4':prob_gt_4})
                natcatdict[grbname].update({'PROB_Z_GT_5':prob_gt_5})
            except:
                print 'Cannot load %s ' % (fullout)
                continue
    return natcatdict
