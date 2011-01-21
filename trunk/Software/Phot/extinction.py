#!/bin/env python

"""
extinction
  Grabs the extinction in a given direction off
  the NED calculator:

    http://nedwww.ipac.caltech.edu/forms/calculator.html
    See also 
    http://www.astro.princeton.edu/~schlegel/dust/dustpub/CodeC/README.C
"""

__version__ = "0.2"
__author__  = "jbloom@astro.berkeley.edu"

import httplib
import urllib
import time
import sys

urls = [ \
    ['nedwww.ipac.caltech.edu',\
     '/forms/calculator.html','/cgi-bin/nph-calc',
     '<pre>',2,0],\
    ['ned.ipac.caltech.edu',\
     '/forms/calculator.html','/cgi-bin/nph-calc',
     '<pre>',2,0],\
    ]

headers = {"Content-type": \
           "application/x-www-form-urlencoded", \
           "Accept": "text/plain"}

def extinction(lon="12:12:12",lat="-12:12:12",\
               in_equinox="J2000.0",\
               system_in="Equatorial",out_equinox="J2000.0",\
               system_out="Galactic",obs_epoch="2005.0"):

    ind = 0
    for url in urls:
        # try to connect
        conn = httplib.HTTPConnection(url[0])
        conn.request("GET", url[1])
        r1 = conn.getresponse()
        if (r1.status == 200) and (r1.reason == "OK"):
            # this site is active
            break
        
        ind += 1
    if ind > len(urls) + 1:
        print "Sorry. No URLs with extinction are responding."
        if __name__ == "__main__":
            sys.exit()
        else:
            return -1

    conn.close()
    ## ok, we're going to use the url with index ind
    params = urllib.urlencode(\
        {'in_csys': system_in,'in_equinox': in_equinox, \
         'obs_epoch': obs_epoch, 'lon': lon, 'lat': lat,
         'pa': 0.0, 'out_csys': system_out, \
         'out_equinox': out_equinox})
    
    
    conn = httplib.HTTPConnection(url[0])
    conn.request("POST", url[2], params, headers)
    response = conn.getresponse()
    data = response.read()
    extfield = \
             (((data.split(url[3]))[url[4]]).splitlines())[url[5]]
    ebv = float(((extfield.split("=")[1]).split("mag."))[0])
    dlines = data.splitlines()

    nline = 0
    outline = '? ? ?'
    for l in dlines:
        if ((l.find(system_out) != -1) and (l.find("Output") != -1)):
            outline = dlines[nline + 1] 
        nline += 1

    outl = outline.split()

    nline = 0
    inline = '? ? ?'
    for l in dlines:
        if ((l.find("RA or Longitude") != -1)):
            inline = dlines[nline + 2] 
        nline += 1

    inl = inline.split()
    
    return [ebv,outl[0],outl[1],inl[0],inl[1]]

def usage():

    print """
usage: extinction <RA> <DEC> [options]

  RA, DEC can be in decimal degrees/hours (e.g. 16.378),
                    colon separated (e.g., 12:34:23.5)
                    letter separated (e.g., 12h34m23s)
                    or spaced with quotes (e.g., '12 12 34.1')
  where options are:
     -ie J|B          In equinox   (default val = J (2000))
                                                  B (1950))
     -oe J|B          Out equinox  (default val = J)
     -si ga|ec|eq|sg  In system (default eq = Equatorial
                                         ga = Galactic
                                         sg = Supergalactic
                                         ec = Ecliptic)
     -so ga|ec|eq|sg  Out system (default eq = Galactic)
     -v               Verbose output, shows E(B-V)
                         and approximate extinction in
                         other bands
                         
  """
    sys.exit()
    
if __name__ == "__main__":

    # invoked from the command line
    lon="12:12:12"
    lat="-12:12:12"
    in_equinox="J2000.0"
    system_in="Equatorial"
    out_equinox="J2000.0"
    system_out="Galactic"
    obs_epoch="2005.0"
    verbose = 0
    
    if len(sys.argv) == 1:
        print usage()
        sys.exit(0)
        
    if ((sys.argv[1] == "-h") or (sys.argv[1] == "--h")):
        print usage()
        sys.exit(0)

    lon = (sys.argv[1]).strip()
    lat = (sys.argv[2]).strip()

    eqd  = {"J": "J2000.0", "B": "B1950.0"}
    sysd = {"ga": "Galactic", "sg": "SuperGalactic",
            "eq": "Equatorial", "ec": "Ecliptic"}
    

    if len(sys.argv) > 3:
        try:
            t3 = float(sys.argv[3])
            print "WARNING! Malformed input"
        except:
            pass
        
        opts = sys.argv[3:]
        i = 0
        while (i<len(opts)):
            arg=opts[i]
            if (arg.find("-v") > -1):
                verbose = 1
            elif (arg.find("-ie") > -1):
                tmp = arg[i+1]
                if tmp in eqd.keys():
                    in_equinox = eqd[tmp]
            elif (arg.find("-oe") > -1):
                tmp = arg[i+1]
                if tmp in eqd.keys():
                    out_equinox = eqd[tmp]
            elif (arg.find("-so") > -1):
                tmp = arg[i+1]
                if tmp in sysd.keys():
                    system_out = sysd[tmp]
            elif (arg.find("-si") > -1):
                tmp = arg[i+1]
                if tmp in sysd.keys():
                    system_in = sysd[tmp]
            else:
                pass

            i += 1
    ret = extinction(lon=lon,lat=lat,\
                     in_equinox=in_equinox,
                     system_in=system_in,
                     out_equinox=out_equinox,
                     system_out=system_out,
                     obs_epoch=obs_epoch)
    
    if type(ret) == type(int):
        print "Malformed query. Try again"
        sys.exit()

    if ret == None:
        sys.exit()

    if len(ret) != 5:
        print "Bad return %s" % repr(ret)
        sys.exit()
        
    if verbose:
        filtd = {\
            'U': [3372,5.434],\
            'B': [4404,4.315],\
            'g': [5244,3.476],\
            'V': [5428,3.315],\
            'R': [6509,2.673],\
            'r': [6707,2.590],\
            'i': [7985,1.991],\
            'I': [8090,1.940],\
            'z': [9055,1.540],\
            'J': [12660,0.902],\
            'H': [16732,0.576],\
            'K': [22152,0.367]}
            

        print "At location: %s %s (%s; %s) = %s %s (%s; %s):" % \
              (ret[3],ret[4],in_equinox,system_in,\
               ret[1],ret[2],out_equinox,system_out)
        print "E(B-V) = %4.3f mag" % (ret[0])
        print "Extinction in other filters:"
        for f in filtd.keys():
            print "\t%s (%i Ang):\t %6.3f" % \
                  (f,filtd[f][0],filtd[f][1]*ret[0])
        print "See also http://wwwmacho.mcmaster.ca/JAVA/Acurve.html"

    else:
        print ret[0]

    sys.exit()