#!/usr/bin/env python
# encoding: utf-8
"""
Signal.py
Author: Adam Morgan
Created: Oct 1, 2009
	
Collection of functions to watch and extract information from GCN Notices, 
VOEvent feeds, etc.

Initial attempt to monitor for VOEvents, specifically for Swift GRBs,
by checking an RSS feed for updates.  If a suitable RSS entry is found 
(One indicating a Swift GRB Notice, for example), grab it and extract 
information from it.  By just getting the triggerid, one can call the 
functions in ParseGCNNotice to extract useful information about the event; 
Eventually will want to parse the VOEvent xml file itself.

This code can also email a region file to you containing the latest BAT,
XRT, UVOT Positions

'You can't stop the signal.'
"""
import sys
import os
from AutoRedux import send_gmail
from AutoRedux import GRBHTML
from RedshiftMachine import LoadGCN
from MiscBin import qErr
import glob
import time

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def MonitorRSS(feed_url):
    '''
    This function checks to see if a particular RSS entry has already been loaded by
    entering it in a sqlite database.  
    
    To keep checking this feed, put in infinite while loop with a set delay time
    # 
    # while(True):
    #     sql_tuple_list = MonitorRSS("http://feedurl.xml")
    #     time.sleep(60)
    '''
    from time import strftime
    import sqlite3

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
    
    sql_entry_list=[]
    
    rssinst = feedparser.parse(feed_url)
    for entry in rssinst['entries']:
        if True:
            # check for duplicates
            c.execute('select * from RSSContent where url=?', (entry.link,))
            if not c.fetchall():
                if True:
                    xml_file = entry.link  # apparently the entry.link is the address I wanted
#                    print xml_file
                    shortened_link = xml_file
                    try:
                        sql_entry = (entry.link, entry.title, entry.summary, shortened_link)
                        print sql_entry
                        c.execute('insert into RSSContent (`url`, `title`, `content`, `xml_location`) values (?,?,?,?)', sql_entry)
                        sql_entry_list.append(sql_entry)
                    except:
                        print "Could not update RSS database for entry %s" % (entry.title)
            conn.commit()
            
    return sql_entry_list


def SwiftGRBFlow(incl_reg=True,incl_fc=True,\
                mail_reg=False, mail_to='amorgan@berkeley.edu',\
                make_html=True, html_path='/home/amorgan/www/swift/',\
                mail_html=True, feed_type = 'talons', tweet = True, force_mail=False,\
                feed_url="http://www.thinkingtelescopes.lanl.gov/rss/talons_swift.xml",
                out_url_path='http://swift.qmorgan.com/'):
    while(True):
        sql_tuple_list = MonitorRSS(feed_url)
        for sql_tuple in sql_tuple_list:
            print sql_tuple[1]
            entry_title = sql_tuple[1]
            xml_link = sql_tuple[0]
            if feed_type.lower() == 'estar' and entry_title.find('SWIFT') != -1 and entry_title.count('-') == 1:
                # Break up entry title into component parts to extract the triggerid
                splittitle = entry_title.split('_')
                triggerid = splittitle[-1].split('-')[0]
                # If the length of the triggerid is 6, know it is a GRB and not ToO
                # May want to revisit this conditionality to actually parse the VOEvent
            elif feed_type.lower() == 'talons' and entry_title.split()[1] == 'swift':
                # first six characters of 4th value in split link is the triggerid
                splittitle = entry_title.split(' ')
                triggerid = str(splittitle[3])
            else:
                triggerid = '0'
            if len(triggerid) == 6:
                _do_all_trigger_actions(triggerid, incl_reg=incl_reg, incl_fc=incl_fc, mail_reg=mail_reg,mail_to=mail_to, make_html=make_html, html_path=html_path, mail_html=mail_html, feed_type=feed_type, force_mail=force_mail, tweet=tweet, feed_url=feed_url, out_url_path=out_url_path)

        time.sleep(60)

def _do_all_trigger_actions(triggerid,  incl_reg=True,incl_fc=True,\
                        mail_reg=False, mail_to='amorgan@berkeley.edu',\
                        make_html=True, html_path='/home/amorgan/www/swift/',\
                        mail_html=True, feed_type = 'talons', tweet = True, force_mail=False,\
                        feed_url="http://www.thinkingtelescopes.lanl.gov/rss/talons_swift.xml",
                        out_url_path='http://swift.qmorgan.com/'):
    #out_url_path used to be 'http://astro.berkeley.edu/~amorgan/Swift/'
    
    triggerid = triggerid.lstrip('0')
    print 'Loading GCN for trigger %s' % (triggerid)
    gcn = LoadGCN.LoadGCN(triggerid, clobber=True)
    # Eventually want to depreciate the following function
    # and make a generic ds9 region file creating function
            #reg_file_path = gcn.get_positions(create_reg_file=True)
    if incl_reg:
        try:
            reg_path = _incl_reg(gcn)
            if not reg_path: 
                print '\nCOULDNT FIND REG PATH\n'
                qErr.qErr()
        except: qErr.qErr()
    if incl_fc:
        try:
            fc_path = _incl_fc(gcn,last_pos_check=True)
            if not fc_path: 
                print '\nCOULDNT FIND FC PATH\n'
                qErr.qErr()
        except: qErr.qErr()
    if mail_reg:
        try:
            mail_grb_region(gcn,mail_to,reg_path)
        except: qErr.qErr()
    if make_html:
        try:
            grbhtml = make_grb_html(gcn, html_path=html_path, reg_path=reg_path, fc_path=fc_path)
            if mail_html and grbhtml.successful_export:
                _mail_html(gcn,mail_to,clobber=force_mail,tweet=tweet,out_url_path=out_url_path,grbhtml=grbhtml)
        except: qErr.qErr()

def _mail_html(gcn,mail_to,clobber=False,tweet=True,out_url_path='http://swift.qmorgan.com/',grbhtml=None):
    if not hasattr(gcn,'mailed_web'):
        gcn.mailed_web = False
    if not gcn.mailed_web:
        email_subject = 'New Web Page for Swift trigger %i' % int(gcn.triggerid)
        email_body = '''Please visit %s%i/
        for information on this event, updated as more information arrives.''' % (out_url_path,int(gcn.triggerid)) 
        
        # Crap, the following doesn't work because of clobber=True when reloading new GCNs
        # EVENTUALL UPDATE SO WE DON'T RELOAD THE ENTIRE GCN EACH TIME, JUST ADD THE NEW NOTICES?
        # gcn.mailed_web = True
        # send_gmail.domail(mail_to,email_subject,email_body)
        # LoadGCN.SaveGCN(gcn)
        #
        # Temporary hack
        mailchkpath = storepath + '/.mlchk%i' % int(gcn.triggerid)
        if not os.path.exists(mailchkpath) or clobber == True:
            cmd = "echo y > %s" % mailchkpath 
            os.system(cmd)
            print 'Email with web link has not been sent yet; doing that now...'
            send_gmail.domail(mail_to,email_subject,email_body)
            if tweet:
                try:
                    # python-twitter requires some kind of oAuth authentication now which is a pain
                    # so just use tumblr, linked to the q.mailbot account.
                    tumblrmail = '661mafroil@tumblr.com'
                    # import twitter # requires http://code.google.com/p/python-twitter/
                    import tinyurl # requires http://pypi.python.org/pypi/TinyUrl/0.1.0 
                    bigurl = '%s%i/' % (out_url_path,int(gcn.triggerid))
                    littleurl = tinyurl.create_one(bigurl)
                    if grbhtml:
                        ra=str(gcn.best_pos[0]).rstrip('0')
                        dec=str(gcn.best_pos[1]).rstrip('0')
                        uncertainty=str(gcn.best_pos[2]).rstrip('0')
                        pos_label=gcn.best_pos_type,
                        twitsub = "New GRB! Swift Trigger %i" % (int(gcn.triggerid))
                        twittext = ''' *  RA = %s<br> *  Dec = %s<br> *  Uncertainty = %s %s<br> *  Visit %s for more info''' % (ra,dec,uncertainty,pos_label,littleurl)
                        
                    else:
                        twitsub = ''
                        twittext = 'New GRB! Swift Trigger %i. Visit %s for more info.' % (int(gcn.triggerid),littleurl)
                    # api = twitter.Api(username='qmorgan', password='twitme0bafgkm') 
                    # status = api.PostUpdate(twittext)
                    print 'Sending Tweet - %s' % (twittext)
                    send_gmail.domail(tumblrmail,twitsub,twittext,sig=False)
                except: qErr.qErr()
                
        else:
            print 'Email has already been sent for this trigger.'
        
def _incl_reg(gcn,clobber=False):
    searchpath = storepath + '*%s.reg' % str(gcn.triggerid)
    reg_list = glob.glob(searchpath)
    # If a region is found and the latest Notice is not a position, 
    # return the path of the already created region file
    if len(reg_list) == 1 and gcn.last_notice_loaded.find('Position') == -1:
        return reg_list[0]
    # if the latest gcn was a Position type, or if there is no region found,
    # create a new region and return the path
    elif (len(reg_list) == 1 and gcn.last_notice_loaded.find('Position') != -1)\
          or (len(reg_list) == 0):
        reg_path = gcn.get_positions(create_reg_file=True)
        return reg_path
    
def _incl_fc(gcn,src_name='',clobber=False, last_pos_check=False):
    from AutoRedux import qImage
    if not src_name: src_name = 'Swift_' + str(gcn.triggerid)
    searchpath = storepath + '%s_fc.png' % (src_name)
    fc_list = glob.glob(searchpath)
    
    if last_pos_check == True:
        if gcn.last_notice_loaded.find('Position') != -1:
            print 'Last notice was of Position type: %s' % gcn.last_notice_loaded
            clobber = True
            
    # If a fc is found and the latest Notice is not a position, 
    # return the path of the already created finding chart file
    if len(fc_list) != 0 and clobber == False:
        return fc_list[0]
    # if the latest gcn was a Position type, or if there is no f. chart found,
    # create a new f. chart and return the path. This function utilizes the
    # self.best_pos attribute for GCNNotice, if it exists.  
    elif hasattr(gcn,'best_pos'):
        fc_list = qImage.MakeFindingChart(ra=gcn.best_pos[0],dec=gcn.best_pos[1],\
              uncertainty=gcn.best_pos[2],src_name=src_name,pos_label=gcn.best_pos_type,\
              survey='dss2red',cont_str='',size="AUTO")
        fc_path = None
        for path in fc_list:
            if path.find('fc.png') != -1:
                fc_path = path
        return fc_path
    else:
        return None
    

def DownloadFile(base_url,file_name,out_path):
    from urllib2 import Request, urlopen, URLError, HTTPError
    
    #create the url and the request
    url = base_url + file_name
    req = Request(url)
    file_mode = 'b'
    
    # Open the url
    try:
        f = urlopen(req)
        print "downloading " + url
        
        # Open our local file for writing
        local_file = open(out_path, "w" + file_mode)
        #Write to our local file
        local_file.write(f.read())
        local_file.close()

    #handle errors
    except HTTPError, e:
        print "HTTP Error:",e.code , url
    except URLError, e:
        print "URL Error:",e.reason , url
        
def make_grb_html(gcn,html_path='/home/amorgan/www/swift',reg_path=None,fc_path=None):
    '''Create a GRB Page and Return the Path of the created GRB Webpage'''
    
    triggerid = gcn.triggerid
    
    # If gcn doesn't have correct position attributes, set missing attr to None
    pos_list = ['bat_pos','xrt_pos','uvot_pos']
    for pos_type in pos_list:
        if not hasattr(gcn,pos_type): setattr(gcn,pos_type,None)
    
    # Attempt to create grb_time string
    try:
        grb_time = gcn.pdict['grb_date_str'] + ' ' + \
                   gcn.pdict['grb_time_str'],
    except:
        grb_time = None
        
    grbhtml_inst = GRBHTML.MakeGRBPage(triggerid=triggerid,bat_pos=gcn.bat_pos,\
    xrt_pos=gcn.xrt_pos,uvot_pos=gcn.uvot_pos,reg_path=reg_path,\
    fc_path=fc_path, grb_time=grb_time, html_path=html_path)
    
    return grbhtml_inst
    
def mail_grb_region(gcn,mail_to,reg_file_path):
    from AutoRedux import qImage
    
    triggerid = gcn.triggerid 
    
    # Eventually want to depreciate the following - don't want get_positions()
    # to have to be run each time.  I created this function before I did the 
    # more general extract_values(), so it's somewhat vestigial... though maybe keep it
    if not os.path.exists(reg_file_path):
        reg_file_path = gcn.get_positions(create_reg_file=True)
    if not hasattr(gcn,'bat_pos'):
        gcn.get_positions()
    
    from AutoRedux import send_gmail
    
    email_fc = 'amorgan@berkeley.edu'
    email_fc_sub = 'Finding Chart for Swift trigger %i' % int(triggerid)
    email_fc_body = 'Finding Chart for this trigger below.'
    fc_list = []
    
    email_to = mail_to #email_fc
#    if mail_toosci == True: email_to += ' toosci@googlegroups.com'
    email_subject = 'DS9 region files for Swift trigger %i' % int(triggerid)
    email_body = 'Please find the latest region file for this burst below\n\n'

    # check to see if the new region file is actually different from the previous one
    # denoted with a ~.  If it is, then assume the position has been updated, and email it.
    # Make sure a position exists before sending (reg_contents != blank)
#    reg_file_path = gcn.get_positions(create_reg_file = True)
    reg_check_path = reg_file_path+'~'

    reg_contents = 'Contains: '
    if gcn.dict.has_key('Swift-BAT GRB Position'): reg_contents += 'BAT Position'
    if gcn.dict.has_key('Swift-XRT Position'): 
        reg_contents += ', XRT Position'
        source_name = 'Swift_' + str(gcn.triggerid)
        # fc_list = qImage.MakeFindingChart(ra=gcn.pdict['xrt_ra'],dec=gcn.pdict['xrt_dec'],\
        #     uncertainty=gcn.pdict['xrt_pos_err'],src_name=source_name,pos_label='XRT',survey='dss2red')
    if gcn.dict.has_key('Swift-XRT Position UPDATE'): 
        reg_contents += ', XRT Position UPDATE'
        source_name = 'Swift_' + str(gcn.triggerid)
        # fc_list = qImage.MakeFindingChart(ra=gcn.xrt_pos_update[0],dec=gcn.xrt_pos_update[1],\
        #     uncertainty=gcn.xrt_pos_update[2],src_name=source_name,pos_label='XRT upd.',survey='dss2red')
    if gcn.dict.has_key('Swift-UVOT Position'): 
        reg_contents += ', UVOT Position'
        source_name = 'Swift_' + str(gcn.triggerid)
        # fc_list = qImage.MakeFindingChart(ra=gcn.uvot_pos[0],dec=gcn.uvot_pos[1],\
        #     uncertainty=gcn.uvot_pos[2],src_name=source_name,pos_label='UVOT',survey='dss2red')

    email_body += reg_contents
    
    # check to see if there's a finding chart already created.
    fc_path = _incl_fc(gcn)

    # If the region files already exist, add the word "updated" to subject line
    if os.path.exists(reg_check_path):
        email_subject = "UPDATED " + email_subject
    # Make sure position already exists before sending email!!!
    if (reg_contents != 'Contains: ') and hasattr(gcn,'bat_pos'):
        if not os.path.exists(reg_check_path) or \
              (os.path.getsize(reg_check_path) != os.path.getsize(reg_file_path)):
    
            send_gmail.domail(email_to,email_subject,email_body,[reg_file_path])
            if fc_path:
                send_gmail.domail(email_to,email_fc_sub,email_fc_body,[fc_path])
        # Only make copy of region file if it is not blank
        os.system('cp ' + reg_file_path + ' ' + reg_check_path)
