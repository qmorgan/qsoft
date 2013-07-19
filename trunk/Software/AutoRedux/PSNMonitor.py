"""
PSNmonitor.py
Author: Adam Morgan
Created: July 17 2013
	
"""
import sys, os
import time
import datetime
from AutoRedux.PyRSS2Gen import RSS2
from AutoRedux.PyRSS2Gen import RSSItem
from AutoRedux.PyRSS2Gen import Guid
import feedparser
from MiscBin import qErr
from AutoRedux import send_gmail

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def _parse_psn_page(url):
    # download web page (~300kb)
    # split by tables
    # make dictionary with <a name="x"> as key
    # - deal with multiple anames - search if EITHER name is a key
    # if key exists, see if the entry has changed. send "Changed" email
    # if key doesnt exist, add it. send "New" email
    # - save pickle file with updated stuff
    pass
    
def _parse_psn_rss(rrsurl):
    #parse the 
    pass
    
def Monitor_PSN_RSS(feed_url="http://www.cbat.eps.harvard.edu/unconf/tocp.xml"):
    '''
    This function checks to see if a particular RSS entry has already been loaded by
    entering it in a sqlite database.  
    
    To keep checking this feed, put in infinite while loop with a set delay time
    # 
    # while(True):
    #     sql_tuple_list = Monitor_PSN_RSS("http://feedurl.xml")
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
    DATABASE = storepath + 'psn_rss_feed.sqlite'
    conn = sqlite3.connect(DATABASE)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    # Create the table if it doesn't exist
    c.execute('CREATE TABLE IF NOT EXISTS RSSContent (`updated`, `title`, `dateAdded`, `id`, `content`, `url`)')
    
    sql_entry_list=[]
    new_rss_entry_list=[]
    
    rssinst = feedparser.parse(feed_url)
    for entry in rssinst['entries']:
        if True:
            # check for duplicates
            c.execute('select * from RSSContent where updated=?', (entry.updated,)) #should be unique
            if not c.fetchall():
                if True:
                    xml_file = entry.link  # apparently the entry.link is the address I wanted
#                    print xml_file
                    shortened_link = xml_file
                    
                    if not 'link' in entry:
                        errtitle='link value not in RSS entry'
                        qErr.qErr(errtitle=errtitle)
                    if not 'title' in entry:
                        errtitle='title value not in RSS entry'
                        qErr.qErr(errtitle=errtitle)
                    if not 'summary' in entry:
                        errtitle='summary value not in RSS entry; using blank value'
                        print errtitle
                        summary = 'unknown'
                    else:
                        summary = entry.summary
                    if not 'id' in entry:
                        errtitle='id value not in RSS entry; using blank value'
                        print errtitle
                        entryid = 'unknown'
                    else:
                        entryid = entry.id
                        
                    try:
                        sql_entry = (entry.updated, entry.title, entryid, summary, entry.link)
                        print sql_entry
                        c.execute('insert into RSSContent (`updated`, `title`, `id`, `content`, `url`) values (?,?,?,?,?)', sql_entry)
                        sql_entry_list.append(sql_entry)
                        new_rss_entry_list.append(entry)
                    except:
                        qErr.qErr()
                        print "Could not update RSS database for entry %s" % (entry.updated)
            conn.commit()
        
    return new_rss_entry_list
    
def _do_new_entry_actions(new_entry):
    # being fed a parsed rss entry
    psn_id_full=new_entry.id.split('/followups/')[1].strip('"')
    # for some reason, the URL has a space when PSN label gets added
    # http://cbat.eps.harvard.edu/unconf/followups/PSN J15111485+4609115
    #u'PSN J15111485+4609115'
    psn_id = psn_id_full.split()[-1]
    #u'J15111485+4609115'
    
    # check if it's in the pickle file
    # if so, update it - add to summary list
    
    html_body = '<html><body><a href="http://cbat.eps.harvard.edu/unconf/followups/%s">%s</a><br><br>' % (psn_id,psn_id)
    html_body+= new_entry.summary
    html_body+= '<br><br><br></body></html>'
    
    print html_body
    # do email if new
    
    subject = "New Transient %s" % (psn_id_full)
    
    print "Sending email: '%s'" % (subject)
    send_gmail.domail('qmorgan@gmail.com',subject,html_body,html=True)
    
    # do separate email if updated

def PSNFlow():
    while(True):
        feed_url="http://www.cbat.eps.harvard.edu/unconf/tocp.xml"
        new_rss_list = Monitor_PSN_RSS(feed_url)
        if len(new_rss_list) > 10:
            errmsg = 'new_rss_list is huge; not doing trigger actions'
            qErr.qErr(errmsg)
        elif new_rss_list != []:
            for entry in new_rss_list:
                _do_new_entry_actions(entry)
        
        print time.ctime()
        time.sleep(720)
    
    
def MonitorFile(filename=None,maxsleep=100,check_interval=1):
    if not filename:
        filename = loadpath + "latest_grb_time.txt"
    
    f=file(filename,'r')
    old_grb_time_string = f.read().strip()
    f.close()

    longsleep=2
    while True:
        count = 0
        restoverride = 30
        while count < restoverride and count < longsleep:
            time.sleep(check_interval)
            f=file(filename,'r')
            new_grb_time_string=f.read().strip()
            if new_grb_time_string != old_grb_time_string:
                # if the latest_grb_time.txt file has changed, this should mean something has been updated.
                old_grb_time_string = new_grb_time_string
                print "Sense that Files have changed! Downloading faster."
                longsleep = 2
            restoverridedatetime=datetime.datetime.strptime(new_grb_time_string, '%m/%d/%y %H:%M:%S')
            resttimedelta = datetime.datetime.utcnow() - restoverridedatetime
            seconds = resttimedelta.seconds + 86400*resttimedelta.days
            print "It has been %i seconds since last grb. Not downloading for %i seconds" % (seconds, longsleep-count) 
            f.close()
            count += 1
        print  "Downloading file " + str(seconds) # Equivalent to downloading the file    
        longsleep *=2
        # if count < restoverride and restoverride < longsleep:
        #     # count being less than restoverride is due to a change of the timeout file
        #     # treat as equivalent to a GRB going off, so redownload faster.
        #     print "Sense that GRB has gone off! Downloading faster."
        #     longsleep = 2
        if longsleep > maxsleep:
            longsleep = maxsleep

    
