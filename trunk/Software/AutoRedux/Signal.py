#!/usr/bin/env python
# encoding: utf-8
"""
Signal.py
Author: Adam Morgan
Created: Oct 1, 2009
	
Collection of functions to watch and extract information from GCN Notices, 
VOEvent feeds, etc.

'There is only the truth of the signal.'
"""
import sys
import os
from AutoRedux import send_gmail
from AutoRedux import GRBHTML

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'




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
    import RedshiftMachine import LoadGCN
    from AutoRedux import qImage
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
                    
                    if mail_reg == True:
                        mail_region(mail_toosci)
                
                    shortened_link = xml_file
                    try:
                        t = (entry.link, entry.title, strftime("%Y-%m-%d %H:%M:%S", entry.updated_parsed), entry.summary, shortened_link)
                        c.execute('insert into RSSContent (`url`, `title`,`dateAdded`, `content`, `xml_location`) values (?,?,?,?,?)', t)
                    except:
                        print "Could not update database for trigger %i" % int(triggerid)
            conn.commit()


def mail_region(mail_toosci):
    
    regpath = storepath +'sw'+ str(triggerid) + '.reg'
    
    from AutoRedux import send_gmail

    email_adam = 'amorgan@berkeley.edu'
    email_adam_sub = 'Finding Chart for Swift trigger %i' % int(triggerid)
    email_adam_body = 'Finding Chart for this trigger below.'
    fc_list = []
    
    email_to = email_adam
    if mail_toosci == True: email_to += ' toosci@googlegroups.com'
    email_subject = 'DS9 region files for Swift trigger %i' % int(triggerid)
    email_body = 'Please find the latest region file for this burst below\n\n'

    # check to see if the new region file is actually different from the previous one
    # denoted with a ~.  If it is, then assume the position has been updated, and email it.
    # Make sure a position exists before sending (reg_contents != blank)
    reg_file_path = gcn.get_positions(create_reg_file = True)
    reg_check_path = reg_file_path+'~'

    reg_contents = 'Contains: '
    if gcn.dict.has_key('Swift-BAT GRB Position'): reg_contents += 'BAT Position'
    if gcn.dict.has_key('Swift-XRT Position'): 
        reg_contents += ', XRT Position'
        source_name = 'Swift_' + str(gcn.triggerid)
        fc_list = qImage.MakeFindingChart(ra=gcn.pdict['xrt_ra'],dec=gcn.pdict['xrt_dec'],\
            uncertainty=gcn.pdict['xrt_pos_err'],src_name=source_name,pos_label='XRT',survey='dss2red')
    if gcn.dict.has_key('Swift-XRT Position UPDATE'): 
        reg_contents += ', XRT Position UPDATE'
        source_name = 'Swift_' + str(gcn.triggerid)
        fc_list = qImage.MakeFindingChart(ra=gcn.xrt_pos_update[0],dec=gcn.xrt_pos_update[1],\
            uncertainty=gcn.xrt_pos_update[2]*3600.0,src_name=source_name,pos_label='XRT update',survey='dss2red')
    if gcn.dict.has_key('Swift-UVOT Position'): 
        reg_contents += ', UVOT Position'
        source_name = 'Swift_' + str(gcn.triggerid)
        fc_list = qImage.MakeFindingChart(ra=gcn.uvot_pos[0],dec=gcn.uvot_pos[1],\
            uncertainty=gcn.uvot_pos[2]*3600.0,src_name=source_name,pos_label='UVOT',survey='dss2red')

    email_body += reg_contents
                                
    # If the region files already exist, add the word "updated" to subject line
    if os.path.exists(reg_check_path):
        email_subject = "UPDATED " + email_subject
    # Make sure position already exists before sending email!!!
    if (reg_contents != 'Contains: ') and hasattr(gcn,'bat_pos'):
        if not os.path.exists(reg_check_path) or \
              (os.path.getsize(reg_check_path) != os.path.getsize(reg_file_path)):
    
            send_gmail.domail(email_to,email_subject,email_body,[reg_file_path])
            if fc_list != []:
                send_gmail.domail(email_to,email_adam_sub,email_adam_body,fc_list)
        # Only make copy of region file if it is not blank
        os.system('cp ' + reg_file_path + ' ' + reg_check_path)