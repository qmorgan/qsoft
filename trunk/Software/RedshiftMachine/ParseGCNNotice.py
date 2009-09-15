#!/usr/bin/env python
# encoding: utf-8
"""
ParseGCNNotice.py
Author: Adam Morgan
Created: July 29, 2009
Last Updated: Aug 19, 2009
	
This Program parses the online version of the GCN Notices 
posted online (http://gcn.gsfc.nasa.gov/gcn/swift_grbs.html
for a list) and returns a dictionary of relevant parameters.
See, e.g. http://gcn.gsfc.nasa.gov/gcn/other/358422.swift 
for an example list of GCN notices for a particular trigger.
"""
import sys
import os

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

class GCNNotice:
    '''Initializes a list of GCN Notices and populates a dictionary
    with the first set of keys being the title of the Notice and the 
    sub-keys being the entries in the Notice.  
    
    '''
    def __init__(self,triggerid,filetype="NoType"):
        self.filetype = filetype
        self.triggerid = triggerid
        # Be sure to update if the web version has changed!
        # If not already saved on disk, grab the gcn from web
        try:
            self.grabgcnfromweb()
            # Once grabbed from web, create the dictionary
            self.createdict()
            self.successful_load = True
        except ValueError: 
            print "Cannot Load GCN Notice."
            self.successful_load = False
    
    def grabgcnfromweb(self):
        """
        Based on trigger number, goes to the web and grabs
        the GCN information and puts it in a list of strings.
        """
        import urllib2
        try:
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
        except:
            print "Cannot load GCN Notice from web."
    
    def createdict(self):
        '''Creates the dictionary From the web-based GCN notices; this function
        grabs the keys and string values from that GCN notice and puts them
        into a dictionary self.dict  
        '''
        # If we don't already have the gcn list loaded, grab it from the web
        if hasattr(self,'gcn_notices') == False:
            self.grabgcnfromweb()
        commentstring = ''
        # for gcn in gcn_list:
        # 	# Search through each notice to determine type
        # 	gcnsplitlist = gcn.splitlines()
        self.dict={}
        gcn_type_list = []
        type_already_found = []
        gcndict = {}
        add_to_where = 0
        for gcn in self.gcn_notices:
            partialdict = {}
            # Make sure not a empty string and check to make sure it is long enough
            # Q Edits 8/24/09
            gcnsplit = gcn.splitlines()
            if len(gcnsplit) > 3:
                # Find what the notice type is  - 4th line
                typeline = gcnsplit[3]
                typesplit = typeline.split(':     ')
                if typesplit[0] != 'NOTICE_TYPE':
                    print 'THIRD LINE IS NOT NOTICE_TYPE!' 
                    print gcnsplit[2:5]
                    add_to_where += 1
                    # sys.exit()
                else:
                    gcn_type_list.append(typesplit[1])
            else: 
                print "This line is not long enough."
                add_to_where += 1 
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
                # Grab the latest index from the list; + add_to_where because first is ''
                gcn_index = gcn_wherelist[-1] + add_to_where
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
        '''This function gets all available positions (BAT, UVOT, XRT) for a 
        particular burst trigger and creates object attributes for each.  If 
        desired, it will create a DS9 region file for this particular trigger
        and put it in the storage directory.
        '''
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
                elif item.find('XRT') != -1 and item.find('n UPDATE') == -1:
                    self.xrt_pos = self.parse_positions(item)
                elif item.find('XRT Position UPDATE') != -1:
                    self.xrt_pos_update = self.parse_positions(item)
                elif item.find('UVOT') != -1:
                    self.uvot_pos = self.parse_positions(item)
        
        if create_reg_file == True:
            # Creates a ds9 region file
            reg_name = storepath +'sw'+ str(self.triggerid) + '.reg'
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
            if hasattr(self,'xrt_pos_update'):
                tmp_str = 'circle('+str(self.xrt_pos_update[0])+','+str(self.xrt_pos_update[1])\
                    +','+str(self.xrt_pos_update[2]*3600)+'") # color=green text={XRT Update}\n' 
                f.write(tmp_str)
            if hasattr(self,'bat_pos'):
                tmp_str = 'circle('+str(self.bat_pos[0])+','+str(self.bat_pos[1])\
                    +','+str(self.bat_pos[2]*3600)+'") # color=red text={BAT}\n' 
                f.write(tmp_str)
            f.close
            print 'Created region file %s' % reg_name
            return reg_name
    
    def extract_values(self):
        '''The mother function to call all programmed parsing functions for 
        each notice type.  Each notice type will have its own "easyparselist" 
        which says how to extract the numerical values from the dictionary 
        created by createdict(), as well as how to extract any other attributes
        (e.g., from the comments) for that particular trigger.  
        
        The parseable_types attribute ocntains the list of programmed notice 
        types, and the command to call to parse them.
        
        Set up easyparselist for each notice type:
        [key_to_parse,[new_key_name,val_type,split_str,split_ind,lstrip_str,rstrip_str]]
        easyparselist should allow for the extraction of floats from most key items
        by allowing you to split, and then strip, to leave just the number behind.
        val_type is 'f','i', or 's' for float, integer, string
        '''
        if hasattr(self,'dict') == False:
            self.createdict(self.gcn_notices)
            subdict={}
        self.pdict={} # Parsed Dictionary  
        self.parseable_types = {"Swift-BAT GRB Position":"e_bat_pos",\
            "Swift-XRT Position":"e_xrt_pos","Swift-UVOT Position":"e_uvot_pos"}  
        for noticetype,noticedict in self.dict.iteritems():
            if self.parseable_types.has_key(noticetype):
                self.current_comment_string = noticedict['COMMENTS']   
                easyparselist = eval("self."+self.parseable_types[noticetype]+"()")
                self.ext_do_easy_parse(noticetype,noticedict,easyparselist)
                print "Parsed %s" % noticetype
            else:
                print "**Cannot yet parse %s" % noticetype
    
    def e_bat_pos(self):
        easyparselist=\
            [['GRB_DEC',['bat_dec','f','d ',0,'',''] ],\
             ['GRB_RA', ['bat_ra' ,'f','d ',0,'',''] ],\
             ['GRB_INTEN',['bat_inten','f','[cnts]    Image_Peak=',0,'',''],\
                          ['bat_img_peak','f','[cnts]    Image_Peak=',1,'','[image_cnts]']],\
             ['TRIGGER_DUR',['bat_trigger_dur','f','[sec]',0,'','']],\
             ['TRIGGER_INDEX',['bat_trig_ind','f','E_range:',0,'',''],\
                              ['bat_trig_ind_range','s','E_range:',1,'','']],\
             ['BKG_INTEN',['bat_bkg_inten','f','[cnts]',0,'','']],\
             ['BKG_DUR',['bat_bkg_dur','f','[sec]',0,'','']],\
             ['GRB_ERROR',['bat_pos_err','f',' [arcmin',0,'','']]
              ]
        # Now parse the not so trivial to extract values:   
        if self.current_comment_string.find('a rate trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'yes'}
            self.pdict.update(sub_dict)
        elif self.current_comment_string.find('an image trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'no'}   
            self.pdict.update(sub_dict)
        return easyparselist
    
    def e_xrt_pos(self):
        easyparselist=\
            [['GRB_DEC',['xrt_dec','f','d ',0,'',''] ],\
             ['GRB_RA', ['xrt_ra' ,'f','d ',0,'',''] ],\
             ['GRB_INTEN',['xrt_inten','f','[erg/cm2/sec]',0,'',''] ],\
             ['GRB_ERROR',['xrt_pos_err','f',' [arcsec',0,'','']],\
             ['GRB_SIGNIF',['xrt_signif','f',' [sigma]',0,'','']],\
             ['TAM[0-3]',['xrt_tam0','f',' ',0,'',''],\
                          ['xrt_tam1','f',' ',1,'',''],
                          ['xrt_tam2','f',' ',2,'',''],
                          ['xrt_tam3','f',' ',3,'','']
                          ],\
             ['AMPLIFIER',['xrt_amplifier','i','',0,'','']],\
             ['WAVEFORM',['xrt_waveform','i','',0,'','']],\
             ['SUN_DIST',['sun_dist','f','[deg]',0,'','']],\
             ['MOON_DIST',['moon_dist','f','[deg]',0,'','']],\
             ['MOON_ILLUM',['moon_illum','f',' [%]',0,'','']]\
              ]
        # Now parse the not so trivial to extract values:   
        if self.current_comment_string.find('a rate trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'yes'}
            self.pdict.update(sub_dict)
        elif self.current_comment_string.find('an image trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'no'}   
            self.pdict.update(sub_dict)
        return easyparselist
    
    def e_uvot_pos(self):
        easyparselist=\
            [['GRB_DEC',['uvot_dec','f','d ',0,'',''] ],\
             ['GRB_RA', ['uvot_ra' ,'f','d ',0,'',''] ],\
             ['GRB_ERROR',['uvot_pos_err','f',' [arcsec',0,'','']],\
             ['SUN_DIST',['sun_dist','f','[deg]',0,'','']],\
             ['MOON_DIST',['moon_dist','f','[deg]',0,'','']],\
             ['MOON_ILLUM',['moon_illum','f',' [%]',0,'','']]\
              ]
        # Now parse the not so trivial to extract values:   
        if self.current_comment_string.find('a rate trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'yes'}
            self.pdict.update(sub_dict)
        elif self.current_comment_string.find('an image trigger') != -1:
            sub_dict = {'bat_is_rate_trig':'no'}   
            self.pdict.update(sub_dict)
        return easyparselist
    
    def ext_do_easy_parse(self,noticetype,noticedict,easyparselist):
        '''ONLY CALL AS A FUNCTION OF self.extract_values()!!!!
        This does the "simple" parsing based on the easyparselist for each
        notice type.
        '''
        for parse_vals in easyparselist:
            key_to_parse=parse_vals[0]
            for sub_parse_vals in parse_vals[1:]:
                new_key_name=sub_parse_vals[0]
                val_type=sub_parse_vals[1]
                split_str=sub_parse_vals[2]
                split_ind=sub_parse_vals[3]
                lstrip_str=sub_parse_vals[4]
                rstrip_str=sub_parse_vals[5]
                try:
                    value_str = noticedict[key_to_parse]
                    if split_str != '':
                        almost_parsed_value = value_str.split(split_str)[split_ind]
                    else:
                        almost_parsed_value = value_str
                    parsed_value=almost_parsed_value.lstrip(lstrip_str).rstrip(rstrip_str)
                    if val_type=='f':
                        converted_value=float(parsed_value)
                    elif val_type=='i':
                        converted_value=int(parsed_value)
                    elif val_type=='s':
                        converted_value=parsed_value
                    else:
                        print 'NOT A VALID val_type'
                    sub_dict = {new_key_name:converted_value}
                    self.pdict.update(sub_dict)
                except:
                    print "Cannot parse key '%s' in notice type '%s' for trigger %s"\
                           % (key_to_parse,noticetype,self.triggerid)
                    print "Check your definition of easyparselist for this notice type."

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



def grabtriggeridfromrss(mail_reg=False,mail_toosci=False):
    '''Initial attempt to monitor for VOEvents, specifically for Swift GRBs,
    by checking an RSS feed for updates.  If a suitable RSS entry is found 
    (One indicating a Swift GRB Notice, for example), grab it and extract 
    information from it.  By just getting the triggerid, one can call the 
    functions in ParseGCNNotice to extract useful information about the event; 
    Eventually will want to parse the VOEvent xml file itself.
    
    It checks to see if a particular RSS entry has already been loaded by
    entering it in a sqlite database.
    
    This code can also email a region file to you containing the latest BAT,
    XRT, UVOT Positions.  
    
    # To keep checking this feed, put in infinite while loop with a set delay time
    # 
    # while(True):
    #     grabtriggeridfromrss(mail_reg=True)
    #     time.sleep(60)
    '''
    from time import strftime
    import sqlite3
    import LoadGCN
    import qImage
    import glob
    try:
        import feedparser
    except: 
        print "feedparser module not installed"
        print "visit http://www.feedparser.org/"
        sys.exit(1)
    
    # Database management code stolen from http://www.halotis.com/2009/07/01/rss-twitter-bot-in-python/
    DATABASE = storepath + 'gcn_rss_feed.sqlite'
    conn = sqlite3.connect(DATABASE)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    # Create the table if it doesn't exist
    c.execute('CREATE TABLE IF NOT EXISTS RSSContent (`url`, `title`, `dateAdded`, `content`, `xml_location`)')
    gcndict = {}
    
    rssinst = feedparser.parse("http://www.estar.org.uk/voevent/GCN/GCN.rdf")
    for entry in rssinst['entries']:
        if entry.title.find('SWIFT') != -1 and entry.title.count('-') == 1:
            
            # check for duplicates
            c.execute('select * from RSSContent where url=?', (entry.link,))
            if not c.fetchall():
                splitentry = entry.title.split('_')
                splitsplitentry = splitentry[-1].split('-')
                triggerid = splitsplitentry[0]
                # If the length of the triggerid is six, know its a grb and not a ToO
                if len(str(triggerid)) == 6:
    #                xml_file_str = entry.title.split('#')[-1]
    #                xml_file = 'http://www.estar.org.uk/voevent/GCN/nasa.gsfc.gcn/SWIFT/' \
    #                    + xml_file_str + '.xml'
                    xml_file = entry.link  # apparently the entry.link is the address I wanted
                    print xml_file
                    print triggerid
                    # If a new item exists, load the GCN and overwrite any pickle file
                    # That may already exist.
                
                    # There's got to be a better way to check whether a new position was found..
                    # For now, just use this hack
                    # gcn_had_bat_pos = gcndict.has_key('Swift-BAT GRB Position')
                    #         gcn_had_xrt_pos = gcndict.has_key('Swift-XRT Position')
                    #         gcn_had_uvot_pos = gcndict.has_key('Swift-UVOT Position')
                    #         
                    gcn = LoadGCN.LoadGCN(triggerid, clobber=True)
                
                    gcndict = gcn.dict
                    gcn.extract_values()
                    # gcn_has_bat_pos = gcndict.has_key('Swift-BAT GRB Position')
                    #                 gcn_has_xrt_pos = gcndict.has_key('Swift-XRT Position')
                    #                 gcn_has_uvot_pos = gcndict.has_key('Swift-UVOT Position')
                    #                 
                    if (True): # (gcn_had_bat_pos == False and gcn_has_bat_pos == True) or\
                    #                    (gcn_had_xrt_pos == False and gcn_has_xrt_pos == True) or\
                    #                    (gcn_had_uvot_pos == False and gcn_has_uvot_pos == True):
                        if mail_reg == True:
                            
                            regpath = storepath +'sw'+ str(triggerid) + '.reg'
                            
                            from AutoRedux import send_gmail
                
                            email_adam = 'amorgan@berkeley.edu'
                            email_adam_sub = 'Finding Chart for Swift trigger %i' % int(triggerid)
                            email_adam_body = 'Finding Chart for this trigger below.'
                            
                            email_to = email_adam
                            if mail_toosci == True: email_to += ' toosci@googlegroups.com'
                            email_subject = 'DS9 region files for Swift trigger %i' % int(triggerid)
                            email_body = 'Please find the latest region file for this burst below\n\n'
                        
                            reg_contents = 'Contains: '
                            if gcn.dict.has_key('Swift-BAT GRB Position'): reg_contents += 'BAT Position'
                            if gcn.dict.has_key('Swift-XRT Position'): 
                                reg_contents += ', XRT Position'
                                source_name = 'Swift_' + str(gcn.triggerid)
                                fc_list = qImage.MakeFindingChart(ra=gcn.pdict['xrt_ra'],dec=gcn.pdict['xrt_dec'],\
                                    uncertainty=gcn.pdict['xrt_pos_err'],src_name=source_name,pos_label='XRT',survey='dss2red')
                            if gcn.dict.has_key('Swift-XRT Position UPDATE'): reg_contents += ', XRT Position UPDATE'
                            if gcn.dict.has_key('Swift-UVOT Position'): reg_contents += ', UVOT Position'
                        
                            email_body += reg_contents
                                                        
                            # check to see if the new region file is actually different from the previous one
                            # denoted with a ~.  If it is, then assume the position has been updated, and email it.
                            # Make sure a position exists before sending (reg_contents != blank)
                            reg_file_path = gcn.get_positions(create_reg_file = True)
                            reg_check_path = reg_file_path+'~'
                            # If the region files already exist, add the word "updated" to subject line
                            if os.path.exists(reg_check_path):
                                email_subject = "UPDATED " + email_subject
                            # Make sure position already exists before sending email!!!
                            if (reg_contents != 'Contains: ') and hasattr(gcn,'bat_pos'):
                                if not os.path.exists(reg_check_path) or \
                                      (os.path.getsize(reg_check_path) != os.path.getsize(reg_file_path)):
                            
                                    send_gmail.domail(email_to,email_subject,email_body,[reg_file_path])
                                    if fc_list != []:
                                        send_gmail.domail(email_adam,email_adam_sub,email_adam_body,fc_list)
                                # Only make copy of region file if it is not blank
                                os.system('cp ' + reg_file_path + ' ' + reg_check_path)
                
                    shortened_link = xml_file
                    try:
                        t = (entry.link, entry.title, strftime("%Y-%m-%d %H:%M:%S", entry.updated_parsed), entry.summary, shortened_link)
                        c.execute('insert into RSSContent (`url`, `title`,`dateAdded`, `content`, `xml_location`) values (?,?,?,?,?)', t)
                    except:
                        print "Could not update database for trigger %i" % int(triggerid)
            conn.commit()
            