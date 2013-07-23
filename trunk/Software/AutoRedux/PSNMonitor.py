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
from MiscBin.q import sex2dec
from AutoRedux import send_gmail
import datetime
from BeautifulSoup import BeautifulSoup
import urllib2
from MiscBin import qPickle

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

def _get_last_entry():
    '''For testing purposes, get last entry of of the psn feed'''
    last_entry_outpath = storepath + 'psn_last_entry.pkl'
    last_entry = qPickle.load(last_entry_outpath)
    if last_entry == None:
        Monitor_PSN_RSS()
        last_entry = qPickle.load(last_entry_outpath)
    return last_entry

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
    last_entry = rssinst['entries'][0] # saving this for testing purposes
    last_entry_outpath = storepath + 'psn_last_entry.pkl'
    qPickle.save(last_entry,last_entry_outpath,clobber=True)
    duplicate_count = 0
    for entry in rssinst['entries']:
        if duplicate_count < 3:
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
            else:
                duplicate_count += 1
            conn.commit()
        else:
            # break the loop if more than 3 duplicates; really only need to 
            # see one duplicate to break the loop, but adding this just in case
            # (since newer feed entries are at the top, no need to loop through
            # every single one. if there are no new ones, you should know immediately)
            break 
    return new_rss_entry_list
    
def _do_new_entry_actions(new_entry,email='psnmonitor@googlegroups.com'):
    # being fed a parsed rss entry
    psn_id_full=new_entry.id.split('/followups/')[1].strip('"')
    # for some reason, the URL has a space when PSN label gets added
    # http://cbat.eps.harvard.edu/unconf/followups/PSN J15111485+4609115
    #u'PSN J15111485+4609115'
    psn_id = psn_id_full.split()[-1]
    #u'J15111485+4609115'
    
    # check if it's in the pickle file
    # if so, update it - add to summary list
    
    psn_url="http://cbat.eps.harvard.edu/unconf/followups/%s" % (psn_id)
    psn_string = _download_and_obtain_psn_string(psn_url)
    if psn_string != None:
        psn_dict = _parse_psn_format(psn_string)
    else:
        psn_dict = None
    
    if psn_dict:
        pretty_output='''
<br><br>
<table border="0">
<tr><td>Object:</td><td>%s</td></tr>
<tr><td>Designation:</td><td>%s</td></tr>
<tr><td>Discovery date:</td><td>%s</td></tr>
<tr><td>Mag at date:</td><td>%s</td></tr>
<tr><td>Filter:</td><td>%s</td></tr>
<tr><td>RA:</td><td>%s (= %f)</td></tr>
<tr><td>Dec:</td><td>%s (= %f)</td></tr>
<tr><td>Presumed host:</td><td>%s</td></tr>
<tr><td>Offset from host:</td><td>%s, %s (arcsec)</td></tr>
<tr><td>Discoverer:</td><td>%s</td></tr>
<tr><td>Obs. arc:</td><td>%s</td></tr>
</table>
        ''' %  (psn_dict['obj_type'],psn_dict['designation'],
        psn_dict['date_string'].replace(' ','-').replace('2013','UT2013'),
        psn_dict['mag'],psn_dict['filter'],
        psn_dict['ra'],psn_dict['ra_deg'],psn_dict['dec'],psn_dict['dec_deg'],
        psn_dict['locale'],psn_dict['ra_offset'],psn_dict['dec_offset'],
        psn_dict['discoverer'],psn_dict['arc'])
    else:
        pretty_output = 'Cannot parse PSN Message.'
    
    print pretty_output
    
    html_body = '''<html><body>
    <a href="%s">%s</a><br><br>''' % (psn_url,psn_id)
    if psn_dict:
        html_body += psn_dict['dss_html']
        html_body += psn_dict['sdss_html']
        html_body += pretty_output
    html_body+= new_entry.summary
    html_body+= '<br><br><br></body></html>'
    
    print html_body
    # do email if new
    
    subject = "New Transient %s" % (psn_id_full)
    
    print "Sending email: '%s'" % (subject)
    send_gmail.domail(email,subject,html_body,html=True)
    
    # do separate email if updated

def _parse_psn_format(psn_string):
    '''
    Columns 1-21: Designation (three-letter designation in columns 1-3 describes 
    the type of variable: PSN = (possible) supernova; PNV = (possible) nova; 
    TCP = some other type of variable (or unknown).

    Columns 25-39: Date in Universal Time (given as Year Month Date).

    Column 40: Note column (* = discovery observation; # or lower-case letters 
    indicate follow-up observation, with # = single line, a = line 1 of multiple 
    lines, b = line 2 of multiple lines, etc.).

    Columns 43-65: Postion (right ascension and declination, for equinox 2000.0), 
    given to 0s.01 of R.A. and to 0".1 of Decl.

    Columns 68-71: magnitude of object at time specified (column 73 gives the 
    bandpass, thus: U = unfiltered CCD; v = visual; V = CCD V-band; R = CCD 
    R-band; etc.)

    Columns 75-79: offset of potential supernova in R.A. from presumed host 
    galaxy, in arc seconds (maximum 9999), with column 79 the direction 
    (E = east, W = west).

    Columns 80-84: offset of potential supernova in Decl. from presumed host 
    galaxy, in arc seconds (maximum 9999), with column 79 the direction 
    (N = north, S = south).

    Columns 87-95: "Locale", meaning the presumed host galaxy if a non-Milky-Way 
    variable, or the 3-letter IAU constellation abbreviation if a Milky-Way variable. 
    For presumed host galaxies, only use single upper-case letters in column 87, 
    followed by several digits as appropriate, and only give galaxies for these 
    catalogues, and in this order of usage: 
    M = Messier; N = NGC; I = IC; U = UGC; G = MCG; P = PGC; E = ESO.

    Column 97: 1-digit character specifying either the experience of the discoverer 
    (0 = no previous CBAT-confirmed discoveries, 1 = one previous CBAT-confirmed 
    discovery, ..., 9 = nine or more previous CBAT-confirmed discoveries) or the 
    discovery group (given as Roman letters, to be assigned on a case-by-case 
    basis as groups request that they be given a letter code to note that their 
    group made the discovery). A key to such letter codes is given here:
        B = Tom Boles (Coddenham, England)
        C = CHASE program (Cerro Tololo, Chile)
        D = Catalina Real-time Transient Survey
        H = Kamil Hornoch (Ondrejov Observatory, Czech Rep.)
        I = Koichi Itagaki (Yamagata, Japan)
        J = Brazilian Supernovae Search (Cristovao Jacques et al.)
        L = Lick Observatory Supernova Search
        M = Berto Monard (Pretoria, South Africa)
        N = Koichi Nishiyama and Fujio Kabashima (Japan)
        P = Tim Puckett's Supernova Search Program
        R = Guoyuo Sun and Jiangao Ruan (China)
        S = La Sagra Sky Survey (Spain)

    Column 99: 1-digit character to specify the number of separate nights (arc) 
    that the discoverer has positive images of the discovered object (may include 
    images by other observers that the discoverer knows of), with a dash (-) meaning 
    a single image on a single night, a zero (0) meaning multiple images on a 
    single night, a "1" indicating a one-day arc (i.e., two nights), ..., and a 
    "9" indicating an arc of nine or more days.
    
             1         2         3         4         5         6         7         8         9
    123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 

    Object Designation      Date (UT)           R.A. (2000.0) Decl.    Mag. p   Offset    Locale    D A
    PSN J15111485+4609115   2013 07 13.91  *  15 11 14.85 +46 09 11.5  17.8 U   15E   8N  U9761     9 1
    '''
    prefix = psn_string[0:3]
    if prefix == 'PSN': obj_type = '(possible) supernova' 
    elif prefix == 'PNV': obj_type =  '(possible) nova' 
    else: obj_type = 'unknown' 
    
    designation = psn_string[4:21]
    date_string = psn_string[24:39]
    date_ymd = date_string[0:10]
    fraction_of_day = int(date_string.split('.')[-1].strip())/100.
    date_parsed = datetime.datetime.strptime(date_ymd,'%Y %m %d') + datetime.timedelta(fraction_of_day)
    ra = psn_string[42:54]
    dec = psn_string[54:65]
    ra_deg, dec_deg = sex2dec((ra,dec))
    
    mag = psn_string[67:71]
    filt = psn_string[72]
    if filt == 'U': filt = 'unfiltered'
    
    ra_offset_string = psn_string[74:79].strip()
    try:
        ra_offset_value = int(psn_string[74:78])
    except:
        ra_offset_value = 'Unknown'
    ra_offset_direction = psn_string[78]
    dec_offset_string = psn_string[79:84].strip()
    try:
        dec_offset_value = int(psn_string[79:83])
    except:
        dec_offset_value = 'Unknown'
    dec_offset_direction = psn_string[83]
    
    locale = psn_string[86:95]
    if locale.strip() == '':
        locale = 'UNKNOWN'
    
    discoverer_key = psn_string[96]
    disc_dict={
    "B":"Tom Boles (Coddenham, England)",
    "C":"CHASE program (Cerro Tololo, Chile)",
    "D":"Catalina Real-time Transient Survey",
    "H":"Kamil Hornoch (Ondrejov Observatory, Czech Rep.)",
    "I":"Koichi Itagaki (Yamagata, Japan)",
    "J":"Brazilian Supernovae Search (Cristovao Jacques et al.)",
    "L":"Lick Observatory Supernova Search",
    "M":"Berto Monard (Pretoria, South Africa)",
    "N":"Koichi Nishiyama and Fujio Kabashima (Japan)",
    "P":"Tim Puckett's Supernova Search Program",
    "R":"Guoyuo Sun and Jiangao Ruan (China)",
    "S":"La Sagra Sky Survey (Spain)"
    }
    if discoverer_key in disc_dict.keys():
        discoverer_string = disc_dict[discoverer_key]
    elif discoverer_key == '9':
        discoverer_string = 'Unknown observer with nine or more previous CBAT-confirmed discoveries'
    elif discoverer_key in '012345678':
        discoverer_string = 'Unknown observer with %s previous CBAT-confirmed discoveries' % (discoverer_key)
    else:
        discoverer_string = 'Unknown discoverer; cannot parse discoverer_key'
    
    arc_key = psn_string[98]
    if arc_key == '-':
        arc_string = 'A single image on a single night'
    elif arc_key == '0':
        arc_string = 'Multiple images on a single night'
    elif arc_key == '9':
        arc_string = 'An arc of nine or more days'
    elif arc_key in '12345678':
        arc_string = 'A %s day arc' % (arc_key)
    else:
        arc_string = 'Unknown arc; cannot parse arc_key'
    
    dss_url = "http://fc.qmorgan.com/fcserver.py?ra=%f&dec=%f&uncertainty=2&err_shape=combo&incl_scale=yes&size=4&src_name=%s&pos_label=Pos&cont_str=&survey=dss2red" % (ra_deg,dec_deg,designation)
    dss_html = "<a href='%s'>DSS Finding Chart</a><br>" % (dss_url)
    sdss_url = "http://fc.qmorgan.com/fcserver.py?ra=%f&dec=%f&uncertainty=2&err_shape=combo&incl_scale=yes&size=4&src_name=%s&pos_label=Pos&cont_str=&survey=sdss" % (ra_deg,dec_deg,designation)
    sdss_html = "<a href='%s'>SDSS Finding Chart</a> (May not be available)<br>" % (sdss_url)
    
    psn_dict = {
    'dss_html':dss_html,
    'sdss_html':sdss_html,
    'ra':ra,
    'dec':dec,
    'ra_deg':ra_deg,
    'dec_deg':dec_deg,
    'prefix':prefix,
    'obj_type':obj_type,
    'designation':designation,
    'date_string':date_string.strip(),
    'date_parsed':date_parsed,
    'mag':mag,
    'filter':filt,
    'locale':locale,
    'discoverer':discoverer_string,
    'arc':arc_string,
    'psn_string':psn_string,
    'dec_offset':dec_offset_string,
    'ra_offset':ra_offset_string
    }    
    
    return psn_dict

def _download_and_obtain_psn_string(followup_url):
    try:
        sock = urllib2.urlopen(followup_url)
    except:
        errmsg = '%s does not exist. Returning None.' % (followup_url)
        qErr.qErr(errtext=errmsg)
        return None
    html = sock.read()
    sock.close()
    soup=BeautifulSoup(html)
    psn = soup.find('pre')
    if psn == None:
        errmsg = 'Cannot find PSN string in %s. Returning None.' % (followup_url)
        qErr.qErr(errtext=errmsg)
        return None
    elif len(psn.contents) > 1:
        errmsg = '%s has more than one PSN String. Look into this.' % (followup_url)
        qErr.qErr(errtext=errmsg)
    psn_string = psn.contents[0]
    psn_string = psn_string.strip()
    if len(psn_string) != 99:
        errmsg = 'PSN string is of the wrong length. Expected 99, got %i' % len(psn_string)
        errmsg += ' %s\n%s\n' % (followup_url,psn_string)
        qErr.qErr(errtext=errmsg)
        return None
    return str(psn_string)
    
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
        time.sleep(1800)
    
    
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

    
