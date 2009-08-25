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
from RedshiftMachine import ParseSwiftCat
from RedshiftMachine import LoadGCN

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def collect():
    swiftcatdict = ParseSwiftCat.parseswiftcat(storepath+'grb_table_1250801097.txt')
    collected_dict = {}
    for grb_str,catdict in swiftcatdict.iteritems():
        # If the swiftcat has a Redshift associated with it, grab the trigger id
        # For now, only collect if it has an associated redshift
        if 'z' in catdict:
            trigid_str = catdict['triggerid_str']
            print '\nNow loading ', trigid_str
            try:
                triggerid=int(trigid_str)
                loaded_gcn = LoadGCN.LoadGCN(triggerid)
                loaded_gcn.extract_values()
            except:
                print "Cannot load trigger %s for GRB %s" % (trigid_str,grb_str)
            catdict.update(loaded_gcn.pdict)
            
            subdict = {grb_str:catdict}
            collected_dict.update(subdict)
    
    print len(collected_dict), ' entries in the collected dictionary'        
    return collected_dict
        


if __name__ == '__main__':
    collect()
    sys.exit(0)     