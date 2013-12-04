#!/usr/bin/env python

import pyfits
import sys
import os
'''
Emulates MissFITS header replacement behaviour using pyfits.  
Assumes the following configuration:
HEADER SUFFIX	.head
OUTFILE TYPE	SAME
SAVE TYPE	    BACKUP
'''
def usage():

    print """
    usage: missfits.py filename.fits

    Requires a filename.head file output by Scamp
    Will move filename.fits to filename.back, and 
    change the headers in filename.fits
  
  """
    sys.exit()

def MissFITS(fname):
    if not os.path.exists(fname):
        print 'File %s does not exist. Exiting.' % (fname)
        sys.exit()
    cmd = 'cp %s %s.back' % (fname, fname)
    
    headname = fname.rstrip('fits') + 'head'
    hdulist = pyfits.open(fname, mode='update')
    prihdr = hdulist[0].header
    
    openedfile = open(headname)
    lines = openedfile.readlines()
    for line in lines:
        keyvalcomm = line.split('/')
        if len(keyvalcomm) < 2:
            print "not parsing " + line
            continue
        comment_list = keyvalcomm[1:]
        comment = comment_list[0]
        # Rest of comment is truncated?
        keyval = keyvalcomm[0]
        key = keyval.split('=')[0].replace(' ','')
        val = keyval.split('=')[1].replace(' ','')
        if val == '1':
            val = 1
        elif val == '0':
            val = 0
        elif val == 'T':
            val = True
        elif val == 'F':
            val = False
        else:
            try:
                val = float(val)
            except:
                try:
                    val = val.strip("'")
                except:
                    print 'Warning - cannot parse'
        prihdr.update(key,val,comment)
    hdulist.close()
    
if __name__ == "__main__":

    # invoked from the command line
    if len(sys.argv) == 1 or len(sys.argv) > 2:
        print usage()
        sys.exit(0)
        
    if ((sys.argv[1] == "-h") or (sys.argv[1] == "--h")):
        print usage()
        sys.exit(0)
        
    fname = sys.argv[1]
    MissFITS(fname)
    sys.exit()