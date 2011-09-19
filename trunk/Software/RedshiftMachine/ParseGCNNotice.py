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
import time
from AutoRedux import send_gmail
from MiscBin.q import where
from MiscBin import qErr

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
    def __init__(self,triggerid,filetype="NoType",clobber=False):
        self.filetype = filetype
        self.triggerid = triggerid
        self.clobber = clobber
        self.parsed_types = []
        self.available_types = []
        # Be sure to update if the web version has changed!
        # If not already saved on disk, grab the gcn from web
        try:
            self.grabgcnfromweb()
            # Once grabbed from web, create the dictionary
            self.createdict()
            self.successful_load = True
        except ValueError: 
            errtitle="Cannot Load GCN Notice for trigger %s" % str(triggerid)
            qErr.qErr(errtitle=errtitle)
            self.successful_load = False
        except AttributeError:
            errtitle="ATTRIBUTE ERROR. Cannot Create dictionary for trigger %s?" % str(triggerid)
            qErr.qErr(errtitle=errtitle)    
            self.successful_load = False
    
    def grabgcnfromweb(self):
        """
        Based on trigger number, goes to the web and grabs
        the GCN information and puts it in a list of strings.
        """
        import urllib2
        gcnpath = storepath + 'gcn_notices/' + str(self.triggerid) + '.swift'
        gcnappendpath = storepath + 'gcn_notices/' + str(self.triggerid) + '.swift.append'
        if self.clobber or not os.path.exists(gcnpath):
            try:
                gcnaddress = 'http://gcn.gsfc.nasa.gov/gcn/other/%s.swift' % str(self.triggerid)
                gcnwebsite = urllib2.urlopen(gcnaddress)
                gcnstring = gcnwebsite.read()
                f = file(gcnpath,'w')
                f.write(gcnstring)
                f.close()
                thedelimiter = '//////////////////////////////////////////////////////////////////////'
                if os.path.exists(gcnappendpath):
                    g = file(gcnappendpath,'r')
                    gcn_append = g.read() 
                    gcnstring += '\n' + thedelimiter + '\n' + gcn_append
                    g.close()
                # Split up the text file by GCN Notice - make it a list of strings
                gcn_notices = gcnstring.split(thedelimiter)
                num_of_gcns = len(gcn_notices)
                self.gcn_notices = gcn_notices
                self.num_of_gcns = num_of_gcns
                print "Finished loading GCN Notices from web for trigger %s" % self.triggerid
            except:
                errtitle="Cannot load GCN Notice from web."
                errtext = "Cannot load GCN Notice %s" % self.triggerid
                qErr.qErr(errtitle=errtitle,errtext=errtext)
        else:
            try:
                f = file(gcnpath,'r')
                gcnstring = f.read()
                thedelimiter = '//////////////////////////////////////////////////////////////////////'
                # Split up the text file by GCN Notice - make it a list of strings
                if os.path.exists(gcnappendpath):
                    g = file(gcnappendpath,'r')
                    gcn_append = g.read() 
                    gcnstring += '\n' + thedelimiter + '\n' + gcn_append
                    g.close()
                gcn_notices = gcnstring.split(thedelimiter)
                num_of_gcns = len(gcn_notices)
                self.gcn_notices = gcn_notices
                self.num_of_gcns = num_of_gcns
                f.close()
                print "Finished loading GCN Notices from web for trigger %s" % self.triggerid
            except:
                errtitle = "Cannot read GCN Notice from file %s" % (gcnpath)
                qErr.qErr(errtitle=errtitle) 
    
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
        self.good_gcn_notices = []
        for gcn in self.gcn_notices:
            partialdict = {}
            
            # Pre July 2005 there was a change in the format of the GCN Notices
            # However, these values are all the same in the TDRSS. The true
            # Units are c/s according to the TDRSS help, even though it says otherwise
            if '[cnts/sec]' in gcn:
                gcn = gcn.replace('cnts/sec','image_cnts')
                gcn = gcn.replace(' Peak=',' Image_Peak=')
            
            # Make sure not a empty string and check to make sure it is long enough
            # Q Edits 8/24/09
            
            gcnsplit = gcn.splitlines()
            if '' in gcnsplit: gcnsplit.remove('')
            if ' ' in gcnsplit: gcnsplit.remove(' ')
            if len(gcnsplit) > 2:
                # Find what the notice type is  - 3rd line
                typeline = gcnsplit[2]
                typesplit = typeline.split(':     ')
                if typesplit[0] != 'NOTICE_TYPE':
                    print 'THIRD LINE IS NOT NOTICE_TYPE!' 
                    print gcnsplit[0:5]
                    qErr.qErr(errtitle='Third line is not notice type!',errtext=gcn)
                    raise Exception
                else:
                    # Append the GCN to the list of well-formatted notices
                    self.good_gcn_notices.append(gcn)
                    # Append the type of the GCN to the gcn type list
                    gcn_type_list.append(typesplit[1])
            else: 
                print "Split Line in GCN web Page is not long enough."
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
                gcn_index = gcn_wherelist[-1] #+ add_to_where 
                if typecount > 1:
                    print "%s instances of %s found; choosing the latest" % (typecount, gcn_type)
                    type_already_found.append(gcn_type)
                else:
                    print "1 instance of %s found" % gcn_type
                gcnstring = self.good_gcn_notices[gcn_index]
                gcnlines = gcnstring.splitlines()
                for line in gcnlines:
                    # Strip to avoid splitting issues with ':  '
                    line = line.strip()
                    linelist = line.split(':  ')
                    if len(linelist) > 2:
                        print 'SOMETHINGS WRONG - we have two instances of ":  " in this line on the gcn'
                        print line
                    if len(linelist) == 2:
                        # Add to dictionary
                        if linelist[0] != 'COMMENTS':
                            subdict = {linelist[0]:linelist[1].strip()}
                            partialdict.update(subdict)
                        if linelist[0] == 'COMMENTS':
                            commentstring += linelist[1].strip() + ';'
                            subdict = {'COMMENTS':commentstring}
                            partialdict.update(subdict)
                
                ########### Error checking############            
                print(partialdict['NOTICE_TYPE'],gcn_type)
                if not partialdict['NOTICE_TYPE'] == gcn_type:
                    qErr.qErr(errtitle='Notice Types do not match!')
                    raise Exception
                ######################################
                    
                subdict = {gcn_type:partialdict}
                self.dict.update(subdict)

                self.last_notice_loaded = gcn_type
        print "Finished populating dictionary for trigger %s" % self.triggerid
    
    def parse_positions(self, notice_type):
        '''Given the desired GCN Notice type, parses to find the positions 
        contained therein.  Returns a tuple of the RA, Dec, and positional
        error.  All units are in decimal degrees.
        '''
        if hasattr(self,'dict') == False:
            self.createdict()
        if self.dict.has_key(notice_type):
            decstr = self.dict[notice_type]['GRB_DEC']
            rastr = self.dict[notice_type]['GRB_RA']
            errstr = self.dict[notice_type]['GRB_ERROR']
            dec = float(decstr.split('d ')[0])
            ra = float(rastr.split('d ')[0])
            err = float(errstr.split(' [')[0])
            if errstr.find('arcsec') != -1:
                err = err#/3600
            elif errstr.find('arcmin') != -1:
                err = err*60#/60
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
        
        Somewhat vestigial!  Might want to incorporate this and all calls to 
        it into the pdict function.
        
        '''
        if hasattr(self,'dict') == False:
            self.createdict()
        # Maybe put these all in a big position dictionary instead?
        # if there is an update, overwrite 
        pos_list = ['Swift-BAT GRB Position','Swift-XRT Position', \
            'Swift-UVOT Position','Swift-XRT Position UPDATE', \
            'Swift-BAT GRB Position UPDATE', 'Swift-UVOT Position UPDATE']
        for item in pos_list:
            if self.dict.has_key(item):
                if item.find('BAT') != -1:
                    self.bat_pos = self.parse_positions(item)
                    if not hasattr(self,'best_pos') or (self.bat_pos[2] < self.best_pos[2]):
                        self.best_pos = self.bat_pos
                        self.best_pos_type = 'BAT'     
                elif item.find('XRT') != -1 and item.find('n UPDATE') == -1:
                    self.xrt_pos = self.parse_positions(item)
                    if not hasattr(self,'best_pos') or (self.xrt_pos[2] < self.best_pos[2]):
                        self.best_pos = self.xrt_pos
                        self.best_pos_type = 'XRT'
                elif item.find('XRT Position UPDATE') != -1:
                    self.xrt_pos_update = self.parse_positions(item)
                    if not hasattr(self,'best_pos') or (self.xrt_pos_update[2] < self.best_pos[2]):
                        self.best_pos = self.xrt_pos_update
                        self.best_pos_type = 'XRT Upd.'
                elif item.find('UVOT') != -1:
                    self.uvot_pos = self.parse_positions(item)
                    if not hasattr(self,'best_pos') or (self.uvot_pos[2] < self.best_pos[2]):
                        self.best_pos = self.uvot_pos
                        self.best_pos_type = 'UVOT'  
        # Put them in to a dictionary:
        sub_dict = {'best_ra':self.best_pos[0],'best_dec':self.best_pos[1],\
            'best_pos_error':self.best_pos[2],'best_pos_type':self.best_pos_type}
        self.pdict.update(sub_dict)
        
        
            
        if create_reg_file == True:
            # Creates a ds9 region file
            reg_name = storepath +'sw'+ str(self.triggerid) + '.reg'
            f=open(reg_name,'w')
            f.write('# Region file format: DS9 version 4.1\n')
            update_time = time.ctime(time.time())
            commentstr1 = '# Created Automatically from Swift GCN Notices on %s\n' % (update_time)
            commentstr2 = '# Contact Adam N. Morgan at qmorgan@gmail.com with problems/comments/suggestions\n'
            f.write(commentstr1)
            f.write(commentstr2) 
            secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
            f.write(secondstr)
            f.write('fk5\n')
            if hasattr(self,'uvot_pos'):
                tmp_str = 'circle('+str(self.uvot_pos[0])+','+str(self.uvot_pos[1])\
                    +','+str(self.uvot_pos[2])+'") # color=blue text={UVOT}\n' 
                f.write(tmp_str)
            if hasattr(self,'xrt_pos'):
                tmp_str = 'circle('+str(self.xrt_pos[0])+','+str(self.xrt_pos[1])\
                    +','+str(self.xrt_pos[2])+'") # color=yellow text={XRT}\n' 
                f.write(tmp_str)
            if hasattr(self,'xrt_pos_update'):
                tmp_str = 'circle('+str(self.xrt_pos_update[0])+','+str(self.xrt_pos_update[1])\
                    +','+str(self.xrt_pos_update[2])+'") # color=green text={XRT Update}\n' 
                f.write(tmp_str)
            if hasattr(self,'bat_pos'):
                tmp_str = 'circle('+str(self.bat_pos[0])+','+str(self.bat_pos[1])\
                    +','+str(self.bat_pos[2])+'") # color=red text={BAT}\n' 
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
            self.createdict()
            subdict={}
        self.pdict={} # Parsed Dictionary  
        self.parseable_types = {"Swift-BAT GRB Position":"e_bat_pos",\
            "Swift-XRT Position":"e_xrt_pos","Swift-UVOT Position":"e_uvot_pos"}  
        for noticetype,noticedict in self.dict.iteritems():
            if noticetype not in self.available_types:
                self.available_types.append(noticetype)
            if self.parseable_types.has_key(noticetype):
                self.current_comment_string = noticedict['COMMENTS']   
                easyparselist = eval("self."+self.parseable_types[noticetype]+"()")
                self.ext_do_easy_parse(noticetype,noticedict,easyparselist)
                print "Parsed %s" % noticetype
                if noticetype not in self.parsed_types:
                    self.parsed_types.append(noticetype)
            else:
                print "**Cannot yet parse %s" % noticetype
        print 'Parsed: %s' % str(self.parsed_types)
        sub_dict = {'notices_parsed':self.parsed_types}   
        self.pdict.update(sub_dict)
        sub_dict = {'notices_available':self.available_types}   
        self.pdict.update(sub_dict)
        
    def e_bat_pos(self):
        easyparselist=\
            [['GRB_DEC',['bat_dec','f','d ',0,'',''] ],\
             ['GRB_RA', ['bat_ra' ,'f','d ',0,'',''] ],\
             ['GRB_DATE', ['grb_date_tjd','i',';   ',0,'','TJD'],\
                          ['grb_date_doy','i',';   ',1,'','DOY'],\
                          ['grb_date_str','s',';   ',2,'','']],\
             ['GRB_TIME', ['grb_time_sod','f','SOD',0,'',''],\
                          ['grb_time_str','s','SOD',1,' {','} UT']],\
             ['GRB_INTEN',['bat_inten','f','[cnts]    Image_Peak=',0,'',''],\
                          ['bat_img_peak','f','[cnts]    Image_Peak=',1,'','[image_cnts]']],\
             ['TRIGGER_DUR',['bat_trigger_dur','f','[sec]',0,'','']],\
             ['TRIGGER_INDEX',['bat_trig_ind','f','E_range:',0,'',''],\
                              ['bat_trig_ind_range','s','E_range:',1,'','']],\
             ['RATE_SIGNIF',['bat_rate_signif','f','[sigma]',0,'','']],\
             ['IMAGE_SIGNIF',['bat_image_signif','f','[sigma]',0,'','']],\
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

