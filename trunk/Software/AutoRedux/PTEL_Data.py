import os, sys
import pyfits
import glob
import time
from RedshiftMachine import ParseSwiftCat

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def RawToDatabase(raw_path,objtype='GRB'):
    '''Given a path to raw pairitel data and and object ID (could be '*'),
    load some basic information from the pairitel data and, if a swift trigger,
    attempt to load extra information about the trigger and put it in a database
    (for now, just a csv file)
    
    How to get all the raw data from lyra?  Maybe just search for all p0-0 files
    with GRB string; this still may be thousands of files..
    
    '''
    swift_cat_path = storepath+'grb_table_1251400549.txt'
    if not os.path.exists(swift_cat_path): print "WARNING: %s does not exist." % (swift_cat_path)
    # Feed it a raw data folder, grab a list of all the raw p0-0.fits files
    swiftcatdict = {}
    swiftcatdict = ParseSwiftCat.parseswiftcat(swift_cat_path)
    
    pteldict = {}
    
    if not os.path.isdir(raw_path): sys.exit('Not a Directory. Exiting.')
    globstr = 'r20*' + objtype + '*p0-0.fits'
    raw_list = glob.glob(globstr)
    for filename in raw_list:
        semester = ''
        burst_time_str = ''
        grb = ''
        time_delta_hours_str = ''
        comments = ''
        
        
        filepath = raw_path + '/' + filename
        # Open the fits file
        hdulist = pyfits.open(filepath)
        # Get the primary header.  For raw ptel data, should only be 1 header.
        prihdr = hdulist[0].header
        # Get the observation ID, e.g. swift-123456, integral-3501
        target_id = prihdr['TRGTNAME']
        splittarget_id = target_id.split('-')
        # If it was able to be split into two, the format is mission-triggerid
        if len(splittarget_id) == 2:
            mission = splittarget_id[0]
            triggerid = splittarget_id[1]
        else:
            mission = 'Unknown'
            triggerid = 'Unknown'
        # GRB.10000.1
        object_id = prihdr['OBJECT']
        
        
        # Get the time of the first (p0-0) observation for a particular ID
        # '2006-09-29 08:45:31.824422'
        ptel_time = prihdr['STRT_CPU']
        # Split off the microseconds at the end since I don't know how to deal with them atm
        ptel_time_split = ptel_time.split('.')
        # The followin is the format in the ptel header with the microseconds stripped off
        fmt = '%Y-%m-%d %H:%M:%S'
        ptel_time_tuple=time.strptime(ptel_time_split[0],fmt)
        # Convert to seconds since the epoch (sse)
        ptel_time_sse = time.mktime(ptel_time_tuple)
        
        
        if target_id not in pteldict:
            targdict = {target_id:{'mission':mission,'triggerid':triggerid,'obs':{object_id:{'first_obs_time_sse':ptel_time_sse,'first_obs_time':ptel_time_split[0]}}}}
            pteldict.update(targdict)
        else:
            pteldict[target_id]['obs'].update({object_id:{'first_obs_time_sse':ptel_time_sse,'first_obs_time':ptel_time_split[0]}})
        
        if 'ptel_time_sse' not in pteldict[target_id]:
            pteldict[target_id].update({'ptel_time_sse':pteldict[target_id]['obs'][object_id]['first_obs_time_sse']})
            pteldict[target_id].update({'ptel_time':pteldict[target_id]['obs'][object_id]['first_obs_time']})
        else: # if new object observation time is less than the old recorded first time, then subtract
            if pteldict[target_id]['ptel_time_sse'] > pteldict[target_id]['obs'][object_id]['first_obs_time_sse']:
                pteldict[target_id]['ptel_time_sse'] = pteldict[target_id]['obs'][object_id]['first_obs_time_sse']
                pteldict[target_id]['ptel_time'] = pteldict[target_id]['obs'][object_id]['first_obs_time']
        
        # If mission == swift, grab info from published Swift Catalog
        if mission == 'swift' and swiftcatdict != {}:
            found_id = False
            # THE FOLLOWING IS A VERY INEFFICIENT LOOP.  But it should work in the interim.
            for grb_str,catdict in swiftcatdict.iteritems():
                # if the triggerids match, then grab the info
                if triggerid == catdict['triggerid_str']:
                    found_id = True
                    grb = grb_str
                    # YYMMDD
                    grb_ymd = grb_str[0:6]
                    # HH:MM:SS.??
                    burst_time = catdict['burst_time_str']
                    burst_time_split = burst_time.split('.')
                    # YYMMDDHH:MM:SS
                    burst_time_toparse = grb_ymd + burst_time_split[0]
                    bt_fmt = '%y%m%d%H:%M:%S'
                    burst_time_tuple = time.strptime(burst_time_toparse,bt_fmt)
                    # Convert to seconds since last epoch
                    burst_time_sse = time.mktime(burst_time_tuple)
                    burst_time_str = time.strftime(fmt,burst_time_tuple)
                    
                    # Get difference from PTEL time from Burst Time in seconds
                    time_delta = ptel_time_sse - burst_time_sse
                    time_delta_hours = time_delta/3600.0
                    time_delta_hours_str = str(time_delta_hours)
                    
                    pteldict[target_id].update({'time_delta':time_delta_hours,'grb_time_sse':burst_time_sse,'grb_time':burst_time_str})
                    
            if found_id == False:
                print 'COULD NOT FIND ID %s in SWIFT CATALOG' % (triggerid)
                    
        else:
            print 'Cannot yet grab extra info for %s.' % (target_id)
        
        
        string_output = '%s,%s,%s,%s,%s,%s,%s,%s' % (semester,grb,object_id,target_id,ptel_time_split[0],burst_time_str,time_delta_hours_str,comments)
        print string_output
        hdulist.close()
        
        
        #TODO: If two object_ids have the same target_id, combine them.
    return pteldict


def testraw2db():
    RawToDatabase('/Users/amorgan/Data/PAIRITEL/tmp/10637/raw/','GRB')