import sys
sys.path.append('/Library/Frameworks/Python.framework/Versions/6.0.0/lib/python2.6/site-packages/feedparser-4.1-py2.6.egg')
sys.path.append('/Users/amorgan/qrepo/trunk/Software/')
from AutoRedux import qRSS

qRSS.HeartBeatCheck(rssurl='http://swift.qmorgan.com/heartbeat.xml',deadtime=1800)
