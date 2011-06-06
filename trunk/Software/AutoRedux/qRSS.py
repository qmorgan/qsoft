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

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

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
    