import glob
import os
import time 
import pyfits 
from MiscBin.q import sex2dec  
from MiscBin.q import sphere_dist

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'



pathbase = os.environ.get('Q_DIR')
anetlocation = pathbase + '/trunk/Software/PTELCoadd/anet.py'

def _GetMostRecent(allday=False):
    '''Will download all the quickobs downloaded from the last 0-10 minutes,
    depending on what time it is run (searches for HH-M?-??).  Once downloaded,
    will return the path to the most recent of those.
    '''
    mostrecentitem = ['None',1]

    # Get the current gmt YYYY-Mmm-DD-HH-M? string
    thishourstr = time.strftime("%Y-%b-%d-%Hh",time.gmtime())
    minutesplitstr = time.strftime("%M",time.gmtime())[0]
    timestr = thishourstr + minutesplitstr
    if allday:
        # grab all the files taken today
        timestr = time.strftime("%Y-%b-%d",time.gmtime())
    # Download the jrr* files from pteld within this hour
    command = 'scp obs@pteld.sao.arizona.edu:/tmp/jrr%s*.fits %s.' % (timestr,storepath)
    os.system(command)
    
    gsearch = '%s/jrr%s*.fits' % (storepath,timestr)
    glist = glob.glob(gsearch)

    for path in glist:
        print path, os.path.getmtime(path)
        if os.path.getmtime(path) > mostrecentitem[1]:
            mostrecentitem[0] = path
            mostrecentitem[1] = os.path.getmtime(path)

    return mostrecentitem[0]

def _FitWcs(filename):
    '''First Run astrometry.net on file and update header.                   
    For some reason, running within python doesn't update the header info 
    if the astrometry is solved.  However, the command line will.
    '''

    os.system("python "+anetlocation+" "+filename)
    # anet.anet([mostrecentitem[0]])

 
def Compare(filename):
    '''For a downloaded PAIRITEL quick reduced file, compare the location
    that the telescope thinks it is at to the actual location determined 
    by astrometry.net.
    '''
    # Load file header info
    try:
        header = pyfits.getheader(filename)
    except: 
        print 'Cannot obtain header info.  Exiting'
        return
    
    # Grab telescope pointing RA and Dec from header
    # in the quickreduced file - RA and DEC are overwritten by some 
    # untrustworthy wcs fit attempt.  Instead, take the RAS and DECS and 
    # Convert by using q.sex2dec.  Verified from raw file that this is the 
    # same that the RA, DEC used to be
    try:
        point_ra_s = header['RAS']
        point_dec_s = header['DECS']
    except:
        print 'File does not contain RAS and DECS (pointing) keywords.'
        return
    
    point_pos = sex2dec((point_ra_s,point_dec_s))
    point_ra = point_pos[0]
    point_dec = point_pos[1]
    
    try: 
        anet_id = header['AN_JOBID']
    except:
        print "Astrometry.net failed to solve the field.  You're on your own."
        return
    
    # Obtain RA and Dec of center pixel? from astrometry.net info
    # Use xy2sky from wcstools.  -d option outputs in decimal degrees.
    xy2sky_command = "xy2sky -n 6 -d %s 128 128" % (filename)
    try:
        (a,b,c) = os.popen3(xy2sky_command)
    except:
        print "Cannot run xy2sky.  Make sure wcstools installed."
        return
    a.close()
    c.close()
    xy2sky_output=b.readlines()

    ref_radec = " ".join([x.split("J2000")[0] for x in xy2sky_output])
    ref_xy    = " ".join([x.split("J2000")[1] for x in xy2sky_output])
    
    fit_ra = float(ref_radec.split("  ")[0])
    fit_dec = float(ref_radec.split("  ")[1])

    print "Telescope pointing: ", point_ra, point_dec
    print "Actual location: ", fit_ra, fit_dec
    dist_asec = sphere_dist(fit_ra,fit_dec,point_ra,point_dec)
    dist_amin = dist_asec/60.0
    print "Distance: %f arcsecs" % (dist_asec)
    print "This is %f arcmin away; half the FOV is 4.25'." % (dist_amin)
    
def PCheck(remove=True):
    ''' First click the thumbnail on status.pairitel to induce a quick reduction
    runs obs@pteld:/home/jbloom/public_html/show_postage.php
    Then run PCheck() to download the latest jrr file and run anet.py to 
    check position
    '''
    filename = _GetMostRecent()
    _FitWcs(filename)
    Compare(filename)
    
    if remove:    
        thishourstr = time.strftime("%Y-%b-%d-%H",time.gmtime())
        cmd = 'rm %s/jrr%s*.fits' % (storepath,thishourstr)
        os.system(cmd)
