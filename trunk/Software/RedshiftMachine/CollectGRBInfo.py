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

from AutoRedux import qImage
import sys
import os
import time
from RedshiftMachine import ParseSwiftCat
from RedshiftMachine import ParseNatCat
from RedshiftMachine import LoadGCN

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

make_finding_charts=True
def_natcats = [storepath+'bat_catalog_07061275.fits',storepath+'bat_catalog_current.fits']

def collect():
    print '\nNow loading Swift Online Catalog Entries'
    swiftcatdict = ParseSwiftCat.parseswiftcat(storepath+'grb_table_1250801097.txt')
    print "Now loading Nat's Catalog"
    natcatdict = ParseNatCat.combine_natcats(natcatlist=def_natcats,remove_zeros=True)
    collected_dict = {}
    failed_gcn_grbs = []
    failed_nat_grbs = []
    failed_finding_charts = []
    for grb_str,catdict in swiftcatdict.iteritems():
        # If the swiftcat has a Redshift associated with it, grab the trigger id
        # For now, only collect if it has an associated redshift
        if 'z' in catdict:
            trigid_str = catdict['triggerid_str']
            print '\nNow collecting GCN entries for trigger %s, GRB %s' % (trigid_str, grb_str)
            try:
                triggerid=int(trigid_str)
                loaded_gcn = LoadGCN.LoadGCN(triggerid)
                loaded_gcn.extract_values()
                if make_finding_charts:
                    try:
                        source_name = 'GRB' + grb_str
                        qImage.MakeFindingChart(ra=loaded_gcn.pdict['xrt_ra'],dec=loaded_gcn.pdict['xrt_dec'],uncertainty=loaded_gcn.pdict['xrt_pos_err'],src_name=source_name,pos_label='XRT',survey='dss2red')
                    except: 
                        failed_finding_charts.append('GRB'+grb_str+' ('+trigid_str+')')
            except:
                print "Cannot load GCN for trigger %s for GRB %s" % (trigid_str,grb_str)
                failed_gcn_grbs.append('GRB'+grb_str+' ('+trigid_str+')')
            catdict.update(loaded_gcn.pdict)
            
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
        
    print ''
    print len(collected_dict), ' entries in the collected dictionary'
    print 'GRBs failed to gather from GCN: ', failed_gcn_grbs   
    print "GRBs failed to gather from Nat's Catalogue: ", failed_nat_grbs  
    print "GRBs failed to obtain finding charts: ", failed_finding_charts           
    return collected_dict

def compare_z(mydict):
    for i in iter(mydict):
        try:
            print '\nGRB',i,'sw:',mydict[i]['z'],'nat:',mydict[i]['Z']
        except:
            pass

def createarff(outdict,keylist=['t90','fluence','peakflux','xrt_column','wh_mag_isupper','v_mag_isupper'],\
                    attributeclass='z_class',classlist=['high_z','medium_z','low_z']):
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
        datastring += outdict[entry][attributeclass]
        datastring += '\n'
        f.write(datastring)
    
    f.close()

if __name__ == '__main__':
    collect()
    sys.exit(0)     