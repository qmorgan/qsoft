#!/usr/bin/env python
# encoding: utf-8
"""
ParseGCNNotice.py
Author: Adam Morgan
Created: July 29, 2009
Last Updated: July 29, 2009
	Created
	
This Program parses the online version of the GCN Notices 
posted online (http://gcn.gsfc.nasa.gov/gcn/swift_grbs.html
for a list) and returns a dictionary of relevant parameters.
See, e.g. http://gcn.gsfc.nasa.gov/gcn/other/358422.swift 
for an example list of GCN notices for a particular trigger.
"""
import sys
import os

class GCNNotice:
    '''Initializes a list of GCN Notices and populates a dictionary
    with the first set of keys being the title of the Notice and the 
    sub-keys being the entries in the Notice.  
    
    '''
    def __init__(self,triggerid,filetype="NoType"):
        self.filetype = filetype
        self.triggerid = triggerid
        self.dict = {}
        # If not already saved on disk, grab the gcn from web
        # Be sure to update if the web version has changed!
        self.grabgcnfromweb()
        # Once grabbed from web, create the dictionary
        self.createdict()
    
    def grabgcnfromweb(self):
        """
        Based on trigger number, goes to the web and grabs
        the GCN information and puts it in a list of strings.
        """
        import urllib2
        
        gcnaddress = 'http://gcn.gsfc.nasa.gov/gcn/other/%s.swift' % str(self.triggerid)
        gcnwebsite = urllib2.urlopen(gcnaddress)
        gcnstring = gcnwebsite.read()
        thedelimiter = '//////////////////////////////////////////////////////////////////////'
        # Split up the text file by GCN Notice - make it a list of strings
        gcn_notices = gcnstring.split(thedelimiter)
        num_of_gcns = len(gcn_notices)
        self.gcn_notices = gcn_notices
        self.num_of_gcns = num_of_gcns
        print "Finished loading GCN Notices from web for trigger %s" % self.triggerid
    
    def createdict(self):
        # If we don't already have the gcn list loaded, grab it from the web
        if hasattr(self,'gcn_notices') == False:
            self.grabgcnfromweb()
        commentstring = ''
        # for gcn in gcn_list:
        # 	# Search through each notice to determine type
        # 	gcnsplitlist = gcn.splitlines()
        gcn_type_list = []
        type_already_found = []
        gcndict = {}
        for gcn in self.gcn_notices:
            partialdict = {}
            # Make sure not a empty string
            if len(gcn) > 3:
                gcnsplit = gcn.splitlines()
                # Find what the notice type is  - 4th line
                typeline = gcnsplit[3]
                typesplit = typeline.split(':     ')
                if typesplit[0] != 'NOTICE_TYPE':
                    print 'THIRD LINE IS NOT NOTICE_TYPE'  
                    sys.exit()
                gcn_type_list.append(typesplit[1])
        # DETERMINE WHAT THE LATEST OF THAT TYPE IS
        # for gcn_type in gcn_type_list:
        for gcn_type in gcn_type_list:
            typecount = gcn_type_list.count(gcn_type)
            # Clearing out strings and sub dictionarys
            partialdict={}
            subdict={}
            commentstring=''
            if where(type_already_found,gcn_type) == []:
                # Use my defined 'where' function to return a list of indices
                gcn_wherelist = where(gcn_type_list,gcn_type)
                # Grab the latest index from the list; +1 because first is ''
                gcn_index = gcn_wherelist[-1] + 1
                if typecount > 1:
                    print "%s instances of %s found; choosing the latest" % (typecount, gcn_type)
                    type_already_found.append(gcn_type)
                else:
                    print "1 instance of %s found" % gcn_type
                gcnstring = self.gcn_notices[gcn_index]
                gcnlines = gcnstring.splitlines()
                for line in gcnlines:
                    # Strip to avoid splitting issues with ':  '
                    line = line.strip()
                    linelist = line.split(':  ')
                    if len(linelist) > 2:
                        print 'SOMETHINGS WRONG'
                    if len(linelist) == 2:
                        # Add to dictionary
                        if linelist[0] != 'COMMENTS':
                            subdict = {linelist[0]:linelist[1].strip()}
                            partialdict.update(subdict)
                        if linelist[0] == 'COMMENTS':
                            commentstring += linelist[1].strip() + ';'
                            subdict = {'COMMENTS':commentstring}
                            partialdict.update(subdict)
                subdict = {gcn_type:partialdict}
                self.dict.update(subdict)
        print "Finished populating dictionary for trigger %s" % self.triggerid
    

    def parse_positions(self, notice_type):
        '''Given the desired GCN Notice type, parses to find the positions 
        contained therein.  Returns a tuple of the RA, Dec, and positional
        error.  All units are in decimal degrees.
        '''
        if hasattr(self,'dict') == False:
            self.createdict(self.gcn_notices)
        if self.dict.has_key(notice_type):
            decstr = self.dict[notice_type]['GRB_DEC']
            rastr = self.dict[notice_type]['GRB_RA']
            errstr = self.dict[notice_type]['GRB_ERROR']
            dec = float(decstr.split('d ')[0])
            ra = float(rastr.split('d ')[0])
            err = float(errstr.split(' [')[0])
            if errstr.find('arcsec') != -1:
                err = err/3600
            elif errstr.find('arcmin') != -1:
                err = err/60
            else:
                sys.exit('Cannot understand positional error type')
            pos_tuple = (ra,dec,err)
            print "Parsed Positions from %s: %s" % (notice_type, str(pos_tuple))
            return pos_tuple
        else:
            pass
            #print "Dictionary does not have GCN Notice Type %s" % notice_type
    
    def get_positions(self, create_reg_file=False):
        if hasattr(self,'dict') == False:
            self.createdict(self.gcn_notices)
        # Maybe put these all in a big position dictionary instead?
        # if there is an update, overwrite 
        pos_list = ['Swift-BAT GRB Position','Swift-XRT Position', \
            'Swift-UVOT Position','Swift-XRT Position UPDATE', \
            'Swift-BAT GRB Position UPDATE', 'Swift-UVOT Position UPDATE']
        for item in pos_list:
            if self.dict.has_key(item):
                if item.find('BAT') != -1:
                    self.bat_pos = self.parse_positions(item)
                elif item.find('XRT') != -1:
                    self.xrt_pos = self.parse_positions(item)
                elif item.find('UVOT') != -1:
                    self.uvot_pos = self.parse_positions(item)
        
        if create_reg_file == True:
            # Creates a ds9 region file
            reg_name = str(self.triggerid) + '.reg'
            f=open(reg_name,'w')
            f.write('# Region file format: DS9 version 4.1\n')
            secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
            f.write(secondstr)
            f.write('fk5\n')
            if hasattr(self,'uvot_pos'):
                tmp_str = 'circle('+str(self.uvot_pos[0])+','+str(self.uvot_pos[1])\
                    +','+str(self.uvot_pos[2]*3600)+'") # color=blue text={UVOT}\n' 
                f.write(tmp_str)
            if hasattr(self,'xrt_pos'):
                tmp_str = 'circle('+str(self.xrt_pos[0])+','+str(self.xrt_pos[1])\
                    +','+str(self.xrt_pos[2]*3600)+'") # color=yellow text={XRT}\n' 
                f.write(tmp_str)
            if hasattr(self,'bat_pos'):
                tmp_str = 'circle('+str(self.bat_pos[0])+','+str(self.bat_pos[1])\
                    +','+str(self.bat_pos[2]*3600)+'") # color=red text={BAT}\n' 
                f.write(tmp_str)
            f.close
            print 'Created region file %s' % reg_name
    

def where(a,val,wherenot=False):
    """
    Analogous to the 'where' function in IDL
    See thread:
    http://mail.python.org/pipermail/python-list/2007-March/431322.html
    
    Returns a list of indices in the list 'a' where a[i] == val
    If wherenot=True, returns an index where a[i]!= val
    e.g. 
    >>> a=[0,1,2,5,42,2]
    >>> where(a,42)
    [4]
    >>> where(a,2)
    [2,5]
    >>> where(a,2,wherenot=True)
    [0,1,3,4]
    >>> where(a,999)
    []
    """
    if wherenot == False:
    	return [i for i in xrange(len(a)) if a[i]==val]
    else:
    	return [i for i in xrange(len(a)) if a[i]!=val]



