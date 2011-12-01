"""
qRSS.py
Author: Adam Morgan
Created: June 5, 2011
	
"""
import sys, os
import time
import datetime
from AutoRedux.PyRSS2Gen import RSS2
from AutoRedux.PyRSS2Gen import RSSItem
from AutoRedux.PyRSS2Gen import Guid
import feedparser
from MiscBin import qErr

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def UpdateRSS(xmlpath,rsstitle="Q's Feed",rsslink='http://qmorgan.com/',
    rssdescription='Unknown',entrytitle='Title',entrylink='http://qmorgan.com/newpage.html',
    entrydescription='Null',clear_rss=False,make_readable=True):
    
    old_rss_file = feedparser.parse(xmlpath)
    items = []
    
    if not clear_rss:
        #populate the old items
        for entry in old_rss_file.entries:
            items.append(RSSItem(
            title=entry['title'],
            link=entry['link'],
            description=entry['description'],
            guid=Guid(entry['link']),
            pubDate=entry['updated']))
    
    #add the new item
    items.append(RSSItem(
         title = entrytitle,
         link = entrylink,
         description = entrydescription,
         guid = Guid(entrylink),
         pubDate = datetime.datetime.now()))
        
    rss = RSS2(
        title = rsstitle,
        link = rsslink,
        description = rssdescription,
        lastBuildDate = datetime.datetime.now(),
        items = items)
    
    xmlopen = open(xmlpath,'w')
    rss.write_xml(xmlopen)
    xmlopen.close()
    
    if make_readable:
        rsstext = file(xmlpath,'r')
        lines = rsstext.readlines()
        rsstext.close()
        rsswrite = file(xmlpath,'w')
        for line in lines:
            newline = line.replace('<item>','\n<item>')
            rsswrite.write(newline)
        rsswrite.close()
    
def HeartBeatServer(xmlpath='/home/amorgan/www/swift/heartbeat.xml',sleeptime=1200):
    '''Run on the computer you want to see is still alive'''
    
    while(True):
        title=str(datetime.datetime.now())
        description='Server is still alive as of ' + title
        UpdateRSS(xmlpath,rsstitle="Q Heartbeat",rsslink='http://qmorgan.com/',
            rssdescription='Feed to ensure the Q Server is up and running',
            entrytitle=title,entrylink='http://qmorgan.com',
            entrydescription=description,clear_rss=True,make_readable=True)
        print 'Heartbeat sent at ' + title
        time.sleep(sleeptime)

def HeartBeatCheck(rssurl='http://swift.qmorgan.com/heartbeat.xml',deadtime=1800):
    
    #Open logging file
    f=open(storepath+'heartbeat.log','a')
    
    
    deadtime = datetime.timedelta(0,deadtime) # Convert deadtime into timedelta object
    rssinst = feedparser.parse(rssurl)
    if len(rssinst['entries']) > 1:
        raise ValueError('Your RSS feed should only have one entry.')
        
    try:
        assert 'bozo_exception' not in rssinst, 'Server might be dead! Cannot load xml page.'
    except:
        errtitle='Cannot load RSS URL %s. Server down?' % (rssurl)
        f.write(errtitle+'\n')
        qErr.qErr(errtitle=errtitle)
        
    updatedtime = datetime.datetime.strptime(rssinst.entries[0]['updated'],'%a, %d %b %Y %H:%M:%S %Z')
    nowtime = datetime.datetime.now()

    delta = nowtime - updatedtime

    print delta
    f.write('At ' + str(nowtime) + ' it has been ' + str(delta)+ ' since heartbeat\n')
    
    #adding four minutes buffer to account for clock differences on computers
    fourminutes = datetime.timedelta(0,240)
    comparetime = nowtime + fourminutes
    
    try:
        asserttext = 'WARNING: updated time seems to be > 4 minutes in the future. Clocks out of sync?'
        assert updatedtime < comparetime, asserttext
    except:
        f.write(asserttext+'\n')
        qErr.qErr(errtitle='Server/Client Clocks out of Sync!')

    try:
        asserttext = 'Server might be dead; has been ' + str(delta) + ' since last heartbeat'
        assert delta < deadtime, asserttext
    except:
        f.write(asserttext+'\n')
        qErr.qErr(errtitle='Server might be dead!')
    f.close()    
    
def HeartBeatMonitor(rssurl='http://swift.qmorgan.com/heartbeat.xml',checktime=600,deadtime=1800):
    '''Run on a different computer to check in on the server
    checktime is how often you want to check the rss feed to make sure it is updated
    deadtime is how hold the latest rss update would have to be in order to trigger an email
    '''
    deadtime = datetime.timedelta(0,deadtime) # Convert deadtime into timedelta object
    
    # TODO: Add ability to check whether it can even access the rss feed 
    
    while(True):
        HeartBeatCheck(rssurl=rssurl,deadtime=deadtime)
        time.sleep(checktime)
        
    
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

    
