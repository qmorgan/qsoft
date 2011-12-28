import os, sys
import glob
import time
from Phot import extinction
from MiscBin import q
from RedshiftMachine import ParseSwiftCat

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

swift_cat_path = storepath+'grb_table_current.txt'

def update_dict(pteldict,burst):
    '''If a burst is specified on pteldict, this function updates pteldict by adding a formated_time key, which gives date in the form of "http://skycam.mmto.arizona.edu/skycam/YYYYMMDD/night_movie.avi, and an extinction key.'''
    # Pierre- add the snippet of code from dictdate below.  Do it on a per-
    # burst basis, and then we can use update_all_dicts to loop over everything
    if pteldict.has_key(burst):
        time1 = pteldict[burst]['ptel_time']
        time2 = time1[0:4]+time1[5:7]+time1[8:10]
        formated = 'http://skycam.mmto.arizona.edu/skycam/'+time2+'/night_movie.avi'
        pteldict[burst]['formated_time'] = formated 
        firstobs = pteldict[burst]['obs'].keys()[0]
        ra = pteldict[burst]['obs'][firstobs]['scope_ra']
        dec = pteldict[burst]['obs'][firstobs]['scope_dec']
        decpos = (ra, dec)
        secpos = q.dec2sex(decpos)
        ext = extinction.extinction(secpos[0], secpos[1])
        # Just take the 0th entry of the returned extinction list; the rest 
        # of the info is redundant positional information.
        pteldict[burst]['extinction'] = ext[0]
    else: print 'Cannot find specified burst in pteldict!'


def update_all_dicts(pteldict):
    '''Update the dictionary for all keys in the pteldict.'''
    for key in pteldict:
        try:
            update_dict(pteldict,key)
        except:
            print 'Cannot update dictionary for %s. Uh oh!' % (key)
        

def RawToDatabase(raw_path,objtype='GRB',pteldict={},swiftcatdict={}):
    '''Given a path to raw pairitel data and and object ID (could be '*'),
    load some basic information from the pairitel data and, if a swift trigger,
    attempt to load extra information about the trigger and put it in a database
    (for now, just a csv file)
    
    How to get all the raw data from lyra?  Maybe just search for all p0-0 files
    with GRB string; this still may be thousands of files..
    
    '''
    import pyfits
    if not os.path.exists(swift_cat_path): print "WARNING: %s does not exist." % (swift_cat_path)
    # Feed it a raw data folder, grab a list of all the raw p0-0.fits files
    if swiftcatdict=={}:
        swiftcatdict = ParseSwiftCat.parseswiftcat(swift_cat_path)
    
#    pteldict = {}
    
    if not os.path.isdir(raw_path): sys.exit('Not a Directory. Exiting.')
    globstr = raw_path + 'r20*' + objtype + '*p0-0.fits'
    raw_list = glob.glob(globstr)
    print raw_list
    for filepath in raw_list:
        semester = ''
        burst_time_str = ''
        grb = ''
        time_delta_hours_str = ''
        comments = ''
         
#        filepath = raw_path + '/' + filename
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
        
        try:
            ra = prihdr['RA']
            dec = prihdr['DEC']
        except:
            ra = 0.00
            dec = 0.00
        
        if target_id not in pteldict:
            targdict = {target_id:{'mission':mission,'triggerid':triggerid,'obs':{object_id:{'scope_ra':ra,'scope_dec':dec,'filename':filepath,'first_obs_time_sse':ptel_time_sse,'first_obs_time':ptel_time_split[0]}}}}
            pteldict.update(targdict)
        else:
            # Check to see if there is already an object id with that name
            # if there is, make a new object id name
            if object_id not in pteldict[target_id]['obs']:
                pteldict[target_id]['obs'].update({object_id:{'scope_ra':ra,'scope_dec':dec,'filename':filepath,'first_obs_time_sse':ptel_time_sse,'first_obs_time':ptel_time_split[0]}})
            else:
                new_obj_id = str(object_id) + '_' + str(ptel_time_sse)
                pteldict[target_id]['obs'].update({new_obj_id:{'scope_ra':ra,'scope_dec':dec,'filename':filepath,'first_obs_time_sse':ptel_time_sse,'first_obs_time':ptel_time_split[0]}})
                
        # if 'ptel_time_sse' not in pteldict[target_id]:
        #     pteldict[target_id].update({'ptel_time_sse':pteldict[target_id]['obs'][object_id]['first_obs_time_sse']})
        #     pteldict[target_id].update({'ptel_time':pteldict[target_id]['obs'][object_id]['first_obs_time']})
        # else: # if new object observation time is less than the old recorded first time, then subtract
        #     if pteldict[target_id]['ptel_time_sse'] > pteldict[target_id]['obs'][object_id]['first_obs_time_sse']:
        #         pteldict[target_id]['ptel_time_sse'] = pteldict[target_id]['obs'][object_id]['first_obs_time_sse']
        #         pteldict[target_id]['ptel_time'] = pteldict[target_id]['obs'][object_id]['first_obs_time']
        
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
                    
                    # pteldict[target_id].update({'time_delta':time_delta_hours,'grb_time_sse':burst_time_sse,'grb_time':burst_time_str})
                    pteldict[target_id].update({'grb_time_sse':burst_time_sse,'grb_time':burst_time_str,'grb_name':grb})
            if found_id == False:
                print 'COULD NOT FIND ID %s in SWIFT CATALOG' % (triggerid)
                    
        else:
            print 'Cannot yet grab extra info for %s.' % (target_id)
        
        string_output = '%s,%s,%s,%s,%s,%s,%s,%s' % (semester,grb,object_id,target_id,ptel_time_split[0],burst_time_str,time_delta_hours_str,comments)
        print string_output
        hdulist.close()
        
        # loop through observations in pteldict to find the earliest one for time purposes
        for target,targdict in pteldict.iteritems():
            min_ptel_time = 2.0E10  # Reset value
            min_ptel_time_string = ''
            time_delta = 0.0
            for observation,obsdict in targdict['obs'].iteritems():
                if obsdict['first_obs_time_sse'] < min_ptel_time:
                    min_ptel_time = obsdict['first_obs_time_sse']
                    min_ptel_time_string = obsdict['first_obs_time']
                    # if we have the GRB time, calculate the time_delta in hours
                    if 'grb_time_sse' in targdict:
                        time_delta = (min_ptel_time - targdict['grb_time_sse'])/3600.0
                        targdict.update({'time_delta':time_delta})
                    targdict.update({'ptel_time_sse':min_ptel_time,'ptel_time':min_ptel_time_string})
        #TODO: If two object_ids have the same target_id, combine them.
    return pteldict

def testraw2db():
    RawToDatabase('/Users/amorgan/Data/PAIRITEL/tmp/10637/raw/','GRB')

def CrawlThruLyraData(basepath='/PAIRITEL/'):
    swiftdict = ParseSwiftCat.parseswiftcat(swift_cat_path)
    rawpaths=[]
    ptel_dict={}
    error_paths=[]
    
    globstr = basepath + 'sem20???/Dir20??-???-??/'
    rawpaths = glob.glob(globstr)
    for path in rawpaths:
        try:
            ptel_dict = RawToDatabase(path,objtype='GRB',pteldict=ptel_dict,swiftcatdict=swiftdict)
        except:
            error_paths.append(path)
    
    globstr = '/Volumes/BR2/Bloom/PAIRITEL-DATA/sem20???/Dir20??-???-??/Raw/'
    rawpaths = glob.glob(globstr)
    for path in rawpaths:
        try:
            ptel_dict = RawToDatabase(path,objtype='GRB',pteldict=ptel_dict,swiftcatdict=swiftdict)
        except:
            error_paths.append(path)
    
    globstr = '/Volumes/BR2/Bloom/PAIRITEL-DATA/sem20???/Dir20??-???-??/raw/'
    rawpaths = glob.glob(globstr)
    for path in rawpaths:
        try:
            ptel_dict = RawToDatabase(path,objtype='GRB',pteldict=ptel_dict,swiftcatdict=swiftdict)
        except:
            error_paths.append(path)
    
    
    print 'error paths:', error_paths
    return ptel_dict

def SwiftTargUnderTime(pteldict,range=(0.0,24.0),savecsv=False):
    count = 0
    countlist=[]
    grblist=[]
    timedeltalist=[]
    grbtimelist=[]
    pteltimelist=[]
    newfilelist=[]
    badtrigger=[]
    for target in pteldict.keys():
        if pteldict[target]['mission'] == 'swift':
            try:
                if range[0] <= pteldict[target]['time_delta'] < range[1]:
                    count += 1
                    countlist.append(target)
                    grblist.append(pteldict[target]['grb_name'])
                    timedeltalist.append(pteldict[target]['time_delta'])
                    grbtimelist.append(pteldict[target]['grb_time'])
                    pteltimelist.append(pteldict[target]['ptel_time'])
            except:
                badtrigger.append(target)
    print '%i Swift Targets observed between %f and %f hours' % (count,range[0],range[1])
    print 'Triggers: ', countlist
    print 'GRB List: ', grblist
    print 'Bad Triggers (not a GRB): ', badtrigger
    zipped = zip(countlist,grblist,grbtimelist,pteltimelist,timedeltalist)
    if savecsv:
        f = open('targundertime.csv','w')
        for line in zipped:
            newline = str(line).lstrip('(').rstrip(')') + '\n'
            f.write(newline)
        f.close()
        print 'file saved to targundertime.csv'
    return zipped
    
def MakeReduxScript(pteldict,trigid,pl3base='/home/ptelreducer/PyPAIRITELfullredux_distribution/',outdirbase='/home/ptelreducer/storage/amorgan/redux/'):
    '''Designed for use on Betsy!
    
    To make a key script for all,
    
    for key in db.keys():
        PTEL_Data.MakeReduxScript(db,key)
    
    copy all the resultant scripts over to betsy /home/ptelreducer/storage/amorgan/redux/scripts/
    
    move the desired script over to /home/ptelreducer/PyPAIRITELfullredux_distribution/
    
    ./script_name.txt
    
    when finished, move script to /home/ptelreducer/storage/amorgan/redux/scripts/finished_scripts
    
    '''
    script_str = 'mkdir %s/%s\n' % (outdirbase,trigid)
    for objid, objidval in pteldict[trigid]['obs'].iteritems():
        myobjid = objid.split('_')[0]
        rawpath = os.path.dirname(objidval['filename']) # directory containing raw files
        ssestr = str(objidval['first_obs_time_sse'])
        fulloutpath = outdirbase + '/' + trigid + '/' + ssestr + '/'
        script_str += 'mkdir %s\n' % (fulloutpath)
        script_str += 'python2.5 %sPyPAIRITELfullredux_local.py -r ptelreducer@lyra.berkeley.edu:%s -o %s -d %s -f >> reduxout_%s.txt\n' % (pl3base,rawpath,myobjid,fulloutpath,trigid)
    filename = storepath + 'reduxscript_'+trigid+'.txt'
    f = open(filename,'w')
    f.write(script_str)
    cmd = 'chmod +x %s' % (filename)
    os.system(cmd)
    f.close()

def CopyRaw(pteldict,trigid=None,grbname=None,outlocation='/Volumes/MyPassport/Data/FullDB/',clobber=False):
    '''
    Specifying either a target ID or a grb name, copy the raw data over to the
    desired location.
    '''
    obsoutlocation = outlocation + '/%s' % (trigid)
    
    if os.path.exists(obsoutlocation) and clobber == False:
        print "Already copied these files and clobber = False; not copying."
        return
    if not trigid and not grbname:
        print 'No trigid or grbname specified; doing nothing'
        return
    if trigid and grbname:
        print 'Please only specify either trigid or grbname. Doing nothing.'
        return
    if trigid:
        print 'transferring.. %s' % (trigid)
        for objid, objidval in pteldict[trigid]['obs'].iteritems():
            print objid
            rawpath = os.path.dirname(objidval['filename']) # directory containing raw files
            # Determine GRB ID - above we might have added a _sse, but now we
            # will split and make sure we are left with just the actual objid 
            myobjid = objid.split('_')[0]
            globstr = rawpath + '/r20*' + myobjid + '*.fits'
            if not os.path.exists(obsoutlocation):
                try:
                    os.mkdir(obsoutlocation)
                except:
                    print 'Cannot make the directory.  Crap.'
            cmd = 'scp -P 10222 amorgan@lyra.berkeley.edu:%s %s/.' % (globstr, obsoutlocation)
            os.system(cmd)
            # globlist = glob.glob(globstr)
            # print globlist
            # for filename in globlist:
            #     try: 
            #         cmd = 'scp amorgan@lyra.berkeley.edu:%s %s' % (filename, outlocation)
            #         os.system(cmd)
            #     except:
            #         print 'Cannot copy file %s' % filename
            
            
    if grbname:
        pass


def RenameRaw(inpaths,outpath='~/Data/PAIRITEL/tmpraw/',id_list=[],newid='GRB.999.1',newdirs=True):
    '''Copy all files in a directory (such as created by copyraw) into 
    a new folder elsewhere, renaming all files to a common GRB ID.  Useful
    when needing to tack observations together for a common reduction
    
    Put a list of all raw observation folders to be combined into inpaths
    
    If newdirs = True, it will further split on days and create folders for each one.
    '''
    if isinstance(inpaths,str):
        inpaths = [inpaths]

    if not isinstance(inpaths,list):
        raise ValueError('inpaths needs to be of type "list"')
    
    if not id_list:
        print 'No ids provided...'
        
    rawfilelist = []
    for inpath in inpaths:
        if not os.path.exists(inpath):
            verr = "Path %s does not exist. Exiting." % (inpath)
            raise ValueError(verr)
        else:
            for grbid in id_list:
                rawfilelist += glob.glob(inpath+'/*'+str(grbid)+'*')

    newfilelist = []
    # extract 
    datelist=[]
    for line in rawfilelist:
        # Grab the dates
        date = os.path.basename(line)[1:12]
        if date not in datelist:
            datelist.append(date)
            
        fbase = os.path.basename(line)
        fdir = os.path.dirname(line)
        
        fs = fbase.split('-')
        # Check if a raw file
        if len(fs) == 7 and fs[0][0] == 'r':
            fbasenew = '%s-%s-%s-%s-%s-%s-%s' % (fs[0],fs[1],fs[2],fs[3],newid,fs[5],fs[6])
            newfilelist.append(fbasenew)
            
            #check the path exists
            if not os.path.exists(outpath):
                try:
                    os.mkdir(outpath)
                except:
                    print 'Cannot make the directory %s.  Crap.' % (outpath)
            
            outfull = outpath + '/' + fbasenew
            
            cmd = 'cp %s %s' % (line, outfull)
            os.system(cmd)
        
        else:
            print 'Line %s did not pass critera; not moving' % (line)    
    
    if newdirs:
        for date in datelist:
            rawfilelist = glob.glob(os.path.dirname(outfull)+'/?'+date+'*')
            newdatepath = outpath + '/' + date
            if not os.path.exists(newdatepath):
                try:
                    os.mkdir(newdatepath)
                except:
                    print 'Cannot make the directory %s. Crap.' % (newdatepath)
                    
            for rawfile in rawfilelist:
                cmd = 'mv %s %s' % (rawfile, newdatepath)
                os.system(cmd)
        

def PlotHist(pteldict):
    '''Plot a histogram of the response times of all available data.'''
    import pylab
    timedeltas=[]
    # Load all the timedeltas into a list
    for target in pteldict.keys():
        if pteldict[target]['mission'] == 'swift':
            try:
                timedeltas.append(pteldict[target]['time_delta'])
            except:
                pass
    pylab.subplot(311)
    pylab.title('PAIRITEL Response Times to Swift Triggers')
    pylab.ylabel('# of Events')
    pylab.hist(timedeltas,bins=[0,1,2,3,4,5,6,7,8,9,10],facecolor='grey',edgecolor='none')
    pylab.hist(timedeltas,bins=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],facecolor='black',edgecolor='none')
    pylab.xlim(0,10)
    pylab.subplot(312)
    pylab.ylabel('# of Events')
    pylab.hist(timedeltas,bins=[0,1,2,3,4,5,6,7,8,9,10],facecolor='grey',edgecolor='none')
    pylab.hist(timedeltas,bins=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],facecolor='black',edgecolor='none')
    pylab.hist(timedeltas,bins=[0,.01,.02,.03,.04,.05,.06,.07,.08,.09,.10],facecolor='red',edgecolor='none')
    pylab.xlim(0,1)
    pylab.ylim(0,25)
    pylab.subplot(313)
    pylab.hist(timedeltas,bins=[0,1,2,3,4,5,6,7,8,9,10],facecolor='grey',edgecolor='none')
    pylab.hist(timedeltas,bins=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],facecolor='black',edgecolor='none')
    pylab.hist(timedeltas,bins=[0,.01,.02,.03,.04,.05,.06,.07,.08,.09,.10],facecolor='red',edgecolor='none')
    pylab.ylim(0,8)
    pylab.xlim(0,0.1)
    pylab.ylabel('# of Events')
    pylab.xlabel('Hours Post-Trigger to First Observation')
    pylab.savefig('PTELHist.eps')
    pylab.show()


def SaveSortedDictTxt(pteldict):
    import pprint
    f = open('lyracrawl.txt','w')
    keys = pteldict.keys()
    keys.sort()
    sorted_vals = map(pteldict.get,keys)
    for item in keys:
        myindex=keys.index(item)
        mystr = '\n\n** %s **\n' % (item)
        f.write(mystr)
        prettydict = pprint.pformat(sorted_vals[myindex])
        f.write(prettydict)
    f.close
