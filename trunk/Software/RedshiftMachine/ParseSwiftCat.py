#!/usr/bin/env python
# encoding: utf-8
"""
ParseSwiftCat.py
Author: Adam Morgan
Created: July 20, 2009
Last Updated: July 20, 2009
    Created
    
This program takes a tab-delimated catalog of Swift parameters from 
http://swift.gsfc.nasa.gov/docs/swift/archive/grb_table.html/
and reads them into a Python Dictionary in preparation for doing some 
redshift estimates based on parameters with machine learning.  

Options on GRB Table lookup:
1:  Swift
2:  All Years
3:  Redshifts
4:  General: Redshift
    BAT Specific: T90 Duration, Fluence, 1-sec Peak Photon Flux
    XRT Specific: Location, Column Density (NH)
    UVOT Specific: V Magnitude, Other Filter Magnitudes
    
Then download table as tab-delimited text file
"""

import csv

def sex2dec(sex_pos):
    '''
    Convert sexagesimal position string tuple into decimal degree tuple
    
    '''
    if type(sex_pos).__name__ != 'tuple':
        print 'Was Epecting Tuple position'
        return
    elif len(sex_pos) != 2:
        print 'Was Expecting tuple of length two'
        return
    elif type(sex_pos[0]).__name__ != 'str' and type(sex_pos[1]).__name__ != 'str':
        print 'Sexagesimal entries need to be of type string'
        return
    elif len(sex_pos[0]) < 8 or len(sex_pos[1]) < 8:
        print "Positions not formatted correctly! ('12:34:56.7','-65:43:21.0')"
        return
    # TODO: Also Check declination    
    elif (sex_pos[0][2] != ':' or sex_pos[0][5] != ':'): 
        if (sex_pos[0][2] != ' ' or sex_pos[0][5] != ' '):
            print "Positions not formatted correctly: ('12:34:56.7','-65:43:21.0')"
            return
    # print 'yay it is formatted correctly'
    
    # Split up string into a list
    if sex_pos[0][2] == ':':
        ra_list=sex_pos[0].split(':')
        dec_list=sex_pos[1].split(':')
    if sex_pos[0][2] == ' ':
        ra_list=sex_pos[0].split(' ')
        dec_list=sex_pos[1].split(' ')
    
    # TODO: Add more error checking to make sure positions are in the right range.
    
    ra_ddeg = float(ra_list[0])*15 + float(ra_list[1])/60 + float(ra_list[2])/3600
    if dec_list[0][0] == '-':  #if it is a negative number
        dec_ddeg = float(dec_list[0]) - float(dec_list[1])/60 - float(dec_list[2])/3600
    else:
        dec_ddeg = float(dec_list[0]) + float(dec_list[1])/60 + float(dec_list[2])/3600
    
    ddeg_pos = (ra_ddeg,dec_ddeg)
    return ddeg_pos


# Read a tab delimited file
bork=csv.reader(open('grb_table_1248114367.txt'),delimiter='\t')
borklist=[]
for row in bork:
    borklist.append(row)

# This creates a list of all the GRBs 
# Check that there are the right number of entries in each list:
# for entries in bork:
#     print len(entries)
# all should come out as the same number

# Format:
# GRB, T90, Fluence, PeakFlux, XRT_RA, XRT_Dec, XRT_Column, UVOT V Mag,
# UVOT Other Mags, Redshift

grbdict = {}

# Now go through the list of objects read in and put them into a dictionary

for grbs in borklist:
    subdict={grbs[0]:{'t90_str':grbs[1],'fluence_str':grbs[2], 'peakflux_str':grbs[3], \
             'xrt_ra_str':grbs[4], 'xrt_dec_str':grbs[5], 'xrt_column_str':grbs[6], \
             'v_mag_str':grbs[7], 'uvot_list':grbs[8], 'z_str':grbs[9]}}
    grbdict.update(subdict)

# Update the dictonary to parse the crap and make it better
for entry in grbdict.keys():
    # Make XRT RA, Dec into decimal degree tuple.  This will create a tuple
    # keyword called 'xrt_pos' which is in decimal degrees.
    sex_pos_tup = (grbdict[entry]['xrt_ra_str'],grbdict[entry]['xrt_dec_str'])
    if sex_pos_tup[0] != 'TBD':
        xrt_pos = {'xrt_pos':sex2dec(sex_pos_tup)}
        grbdict[entry].update(xrt_pos)
    else:
        print 'COULD NOT PARSE XRT_POS for entry %s' % entry
    # TODO: Convert into a distance above the galactic plane

    # Convert T90 to float
    try:
        # Try to convert t90 to a float
        t90 = {'t90':float(grbdict[entry]['t90_str'])}
        grbdict[entry].update(t90)
    except: 
        try:
            # If that didn't work, assume that it starts with a ~ or > and try again
            if grbdict[entry]['t90_str'][0] == '~':
                print 'CONVERTING APPROXIMATE t90 TO ABSOLUTE for %s' % entry
            if grbdict[entry]['t90_str'][0] == '>':
                print 'CONVERTING LOWER LIMIT t90 TO ABSOLUTE for %s' % entry
            t90 = {'t90':float(grbdict[entry]['t90_str'][1:])}
            grbdict[entry].update(t90)
        except:
            print 'COULD NOT PARSE T90 for entry %s' % entry
    
    # Convert fluence to a float
    try:
        # Try to convert fluence to a float
        fluence = {'fluence':float(grbdict[entry]['fluence_str'])}
        grbdict[entry].update(fluence)
    except: 
        try:
            # If that didn't work, assume that it starts with a ~ or > and try again
            if grbdict[entry]['fluence_str'][0] == '~':
                print 'CONVERTING APPROXIMATE fluence TO ABSOLUTE for %s' % entry
            if grbdict[entry]['fluence_str'][0] == '>':
                print 'CONVERTING LOWER LIMIT fluence TO ABSOLUTE for %s' % entry
            fluence = {'fluence':float(grbdict[entry]['fluence_str'][1:])}
            grbdict[entry].update(fluence)
        except:
            print 'COULD NOT PARSE fluence for entry %s' % entry
            
    # Convert peakflux to a float
    try:
        # Try to convert peakflux to a float
        peakflux = {'peakflux':float(grbdict[entry]['peakflux_str'])}
        grbdict[entry].update(peakflux)
    except: 
        try:
            # If that didn't work, assume that it starts with a ~ or > and try again
            if grbdict[entry]['peakflux_str'][0] == '~':
                print 'CONVERTING APPROXIMATE peakflux TO ABSOLUTE for %s' % entry
            if grbdict[entry]['peakflux_str'][0] == '>':
                print 'CONVERTING LOWER LIMIT peakflux TO ABSOLUTE for %s' % entry
            peakflux = {'peakflux':float(grbdict[entry]['peakflux_str'][1:])}
            grbdict[entry].update(peakflux)
        except:
            print 'COULD NOT PARSE peakflux for entry %s' % entry
    
    # Convert xrt_column to a float
    try:
        # Try to convert xrt_column to a float
        xrt_column = {'xrt_column':float(grbdict[entry]['xrt_column_str'])}
        grbdict[entry].update(xrt_column)
    except: 
        try:
            # If that didn't work, assume that it starts with a ~ or < and try again
            if grbdict[entry]['xrt_column_str'][0] == '~':
                print 'CONVERTING APPROXIMATE xrt_column TO ABSOLUTE for %s' % entry
            if grbdict[entry]['xrt_column_str'][0] == '<':
                print 'CONVERTING UPPER LIMIT xrt_column TO ABSOLUTE for %s' % entry
            xrt_column = {'xrt_column':float(grbdict[entry]['xrt_column_str'][1:])}
            grbdict[entry].update(xrt_column)
        except:
            print 'COULD NOT PARSE xrt_column for entry %s' % entry
    
    # Convert UVOT_V_Mag to float & check if its upper limit
    if grbdict[entry]['v_mag_str'][0] == 'V':
        if grbdict[entry]['v_mag_str'][1] == '>':
            v_mag_isupper = {'v_mag_isupper':1}
            grbdict[entry].update(v_mag_isupper)
        elif grbdict[entry]['v_mag_str'][1] == '=':
            v_mag_isupper = {'v_mag_isupper':0}
            grbdict[entry].update(v_mag_isupper)
        elif grbdict[entry]['v_mag_str'][1] == '<':
            v_mag_isupper = {'v_mag_isupper':0}
            grbdict[entry].update(v_mag_isupper)
            print 'We have a LOWER limit for v_mag for grb %s' % entry
        else:
            print 'COULD NOT PARSE v_mag_str for entry %s.  Starts with V though.' % entry
      
    else:
        print 'COULD NOT PARSE v_mag_str for entry %s' % entry
        
    # Convert z to float and check if it is photometric
    # First we see if there's an absolute redshift available that's not 
    # upper limit or approximate. Preferentially use abs or emis over photo.
    # First split up via the | parser that the table gives for different z entries
    z_split = grbdict[entry]['z_str'].split('|')
    for z_ent in z_split:
        # if it's a photometric redshift and there's already a redshift, don't update
        if z_ent.find('photometric') != -1 and grbdict[entry].has_key('z'):
            continue
        else:
            # Split further into a list
            z_split_split = z_ent.split(' ')
            for z_ent_ent in z_split_split:
                # If there is a number in the sub string, assume it is the z
                try:
                    z={'z':float(z_ent_ent)}
                    grbdict[entry].update(z)
                    # If the redshift is photometric, mark it as such.
                    if z_ent.find('photometric') == -1:
                        z_isphot = {'z_isphot':0}
                    else:
                        z_isphot = {'z_isphot':1}
                    grbdict[entry].update(z_isphot)
                    z_isupper = {'z_isupper':0}
                    grbdict[entry].update(z_isupper)
                except:
                    pass
        # If that didn't find us a redshift, try again accepting ~,>,<
        if not grbdict[entry].has_key('z'):
            z_split_split = z_ent.split(' ')
            for z_ent_ent in z_split_split:
                iii = 0
                if z_ent_ent[0] == 'Z' or z_ent_ent[0] == 'z':
                    iii = 1
                if z_ent_ent[iii] == '~':
                    print 'CONVERTING APPROXIMATE redshift TO ABSOLUTE for %s' % entry
                    z_isupper = {'z_isupper':0}
                if z_ent_ent[iii] == '<':
                    print 'CONVERTING UPPER LIMIT redshift TO ABSOLUTE for %s' % entry
                    z_isupper = {'z_isupper':1}
                if z_ent_ent[iii] == '>':
                    print 'CONVERTING LOWER LIMIT redshift TO ABSOLUTE for %s' % entry
                    z_isupper = {'z_isupper':-1}
                # If the redshift is photometric, mark it as such
                if z_ent.find('photometric') == -1:
                    z_isphot = {'z_isphot':0}
                else:
                    z_isphot = {'z_isphot':1}
                try:
                    iii += 1 
                    z = {'z':float(z_ent_ent[iii:])}
                    grbdict[entry].update(z)
                    grbdict[entry].update(z_isphot)
                    grbdict[entry].update(z_isupper)
                except:
                    # cannot do anything..
                    pass 
    
    # Manually insert for the ultra-high-z 090423
    if entry == '090423':
        grbdict[entry]['z'] = 8.2
        grbdict[entry]['z_isphot'] = 0
        grbdict[entry]['z_isupper'] = 0
    
    # Now try again, filling in the approximate redshifts
    if grbdict[entry].has_key('z'):
        pass
    else:
        if not grbdict[entry].has_key('z'):
            print '**COULD NOT PARSE REDSHIFT for entry %s** This is bad' % entry
            
print '**Finished reading in Swift Catalog**'
