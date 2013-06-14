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
    # Update to grab Nat's redshift probabilities - [now do this through ParseNatWeb]
    #batxrtcat = grab_nat_web_data(batxrtcat,filename='/bat/zprob.txt')
    return batxrtcat

def download_file(url,outfile):
    response = urllib2.urlopen(url)
    html = response.read()
    f = file(outfile,'w')
    f.write(html)
    f.close()


def download_nat_web_file(grbname, webtrigid, filename='/bat/zprob.txt', clobber=False):
    
    base_url = 'http://butler.lab.asu.edu/Swift/'
    outpath = storepath + 'grbs/'
    
    fullurl = base_url + webtrigid + filename
    fullout = outpath + grbname + filename
    if grbname == ".DS_Store": return
    
    if not os.path.exists(os.path.dirname(fullout)):
        os.makedirs(os.path.dirname(fullout))
        
        # Try to download the file, if can't, skip to the next one
        try:
            download_file(fullurl,fullout)
            print 'Successful download of  %s to %s' % (fullurl,fullout)
        except:
            print 'Cannot download %s to %s' % (fullurl,fullout)
            return 
    elif clobber == True:
        try:
            if os.path.exists(fullout):
                os.remove(fullout)
            download_file(fullurl,fullout)
            print 'Successful overwrite of  %s to %s' % (fullurl,fullout)
        except:
            print 'Cannot override %s to %s' % (fullurl,fullout)
            return
    else:
        pass
    

def grab_nat_web_data(grbdict, filename='/bat/zprob.txt', clobber=False):
    base_url = 'http://butler.lab.asu.edu/Swift/'
    outpath = storepath + 'grbs/'
        
    # natzprobdict = {}
        
    # for grbname in os.listdir(outpath):
    for grbname in grbdict.keys():
        fullout = outpath + grbname + filename
        if grbname == ".DS_Store": continue
        webtrigid = grbdict[grbname]['webtrigid'].lstrip('id')
        
        # Download the file
        download_nat_web_file(grbname,webtrigid,filename=filename, clobber=clobber)
        # Grab probabilities if we're looking at that filename
        if filename == '/bat/zprob.txt':
            ##z prob
            #0.0 0.0
            #0.1 0.0173421269984
            #0.201000005007 0.0346842539969
            #0.300999999046 0.0685398958379
            try:
                response = f = open(fullout,'r')
                html = f.readlines()            
            except:
                  print 'Cannot load %s. File does not exist. ' % (fullout)
                  continue
            #error checking to see if the file is formatted correctly
            if len(html) < 103:
                print "%s malformed. Too few lines."
                continue
            elif len(html) < 10:
                print "%s malformed. Empty?"
                continue
            elif len(html) > 103:
                print "%s malformed. Too many lines."
                continue  
            z_pred_arr = []
            prob_arr = []
            ind = 0 # Grab the index to mark where redshifts are for summing probs
            for line in html:
                try:
                    int(line[0])
                except:
                    continue #checking if line is malformed
                if line[0] == '#':
                    continue
                z_prob = line.rstrip('\n').split(' ')
                z_pred_arr.append(float(z_prob[0]))
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
        
            # Prob less than is just 1-prob greater than, since these are normalized to 1
            prob_lt_1 = 1 - prob_gt_1
            prob_lt_2 = 1 - prob_gt_2
            prob_lt_3 = 1 - prob_gt_3
            prob_lt_4 = 1 - prob_gt_4
            prob_lt_5 = 1 - prob_gt_5
        
            # get the most probable z (maximum of the probability array)
            try:
                most_prob_z = z_pred_arr[prob_arr.index(max(prob_arr))]
            except(ValueError):
                most_prob_z = numpy.nan
            # natzprobdict[grbname]={}
        
            grbdict[grbname].update({'PROB_Z_GT_1':prob_gt_1})
            grbdict[grbname].update({'PROB_Z_GT_2':prob_gt_2})
            grbdict[grbname].update({'PROB_Z_GT_3':prob_gt_3})
            grbdict[grbname].update({'PROB_Z_GT_4':prob_gt_4})
            grbdict[grbname].update({'PROB_Z_GT_5':prob_gt_5})
        
            grbdict[grbname].update({'PROB_Z_LT_1':prob_lt_1})
            grbdict[grbname].update({'PROB_Z_LT_2':prob_lt_2})
            grbdict[grbname].update({'PROB_Z_LT_3':prob_lt_3})
            grbdict[grbname].update({'PROB_Z_LT_4':prob_lt_4})
            grbdict[grbname].update({'PROB_Z_LT_5':prob_lt_5})
        
            grbdict[grbname].update({'MOST_PROB_Z':most_prob_z})
    return grbdict
