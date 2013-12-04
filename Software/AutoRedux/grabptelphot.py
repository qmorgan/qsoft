#!/usr/bin/env python
# encoding: utf-8
"""
grabptelphot.py

Author: Adam Morgan

Given sky position and data location, grabs the photometry from ptel_astrophot.py
Currently only works with mosaics output directly from the pipeline; will 
fix later to deal with any directory/file
"""

import ptel_astrophot
import time
import os
import sys

class ImgFile:
    '''Finds the Location of the images, moves the files to a writable location'''

    def __init__(self,sem,year,month,day,project,object,obs=1):
        datestr = time.strftime("%Y-%b-%d",(year,month,day,0,0,0,0,1,-1))
	self.dirpath = "/Bloom/PAIRITEL-DATA/%s/Dir%s/outs_2/Mosaics/" % (sem,datestr)
	datestr2 = time.strftime("%Y%b%d",(year,month,day,0,0,0,0,1,-1))
	intfilestr = "%s.%i.%i-%s" % (project,object,obs,datestr2)

	self.imagename_j = "mosj%s.fits" % intfilestr
	self.imagename_h = "mosh%s.fits" % intfilestr
	self.imagename_k = "mosk%s.fits" % intfilestr
	self.weightname_j = "mosj%s_w.fits" % intfilestr
	self.weightname_h = "mosh%s_w.fits" % intfilestr
	self.weightname_k = "mosk%s_w.fits" % intfilestr

	self.imagepath_j = "%s%s" % (self.dirpath, self.imagename_j)
	self.weightpath_j = "%s%s" % (self.dirpath, self.weightname_j)
	self.imagepath_h = "%s%s" % (self.dirpath, self.imagename_h)
	self.weightpath_h = "%s%s" % (self.dirpath, self.weightname_h)
	self.imagepath_k = "%s%s" % (self.dirpath, self.imagename_k)
	self.weightpath_k = "%s%s" % (self.dirpath, self.weightname_k)

	self.tmpdir = "/Volumes/oBloom/amorgan/tmpphot/"

	if not os.path.exists(self.tmpdir):
	    print 'The temp directory "%s" does not exist. Exiting..' % self.tmpdir
	    sys.exit(1)

	if not os.path.exists(self.imagepath_j):
	    print 'The file "%s" does not exist.  Exiting..' % self.imagepath_j
	    sys.exit(1)
	if not os.path.exists(self.imagepath_h):
	    print 'The file "%s" does not exist.  Exiting..' % self.imagepath_h
	    sys.exit(1)
	if not os.path.exists(self.imagepath_k):
	    print 'The file "%s" does not exist.  Exiting..' % self.imagepath_k
	    sys.exit(1)
	
	print '(Initialized ImgFile %s)' % self.imagepath_j
	print '(Initialized ImgFile %s)' % self.imagepath_h
	print '(Initialized ImgFile %s)' % self.imagepath_k

    def MoveFiles(self):
        '''Moves the Images from the default directory to a temp directory we can do photometry in'''
	os.system("cp %s %s." % (self.imagepath_j,self.tmpdir))   #copy j images
	os.system("cp %s %s." % (self.imagepath_h,self.tmpdir))   #copy h images
	os.system("cp %s %s." % (self.imagepath_k,self.tmpdir))   #copy k images

	os.system("cp %s %s." % (self.weightpath_j,self.tmpdir))  #copy j weights
	os.system("cp %s %s." % (self.weightpath_h,self.tmpdir))  #copy h weights
	os.system("cp %s %s." % (self.weightpath_k,self.tmpdir))  #copy k weights

	self.imagepath_j = "%s%s" % (self.tmpdir, self.imagename_j)
	self.weightpath_j = "%s%s" % (self.tmpdir, self.weightname_j)
	self.imagepath_h = "%s%s" % (self.tmpdir, self.imagename_h)
	self.weightpath_h = "%s%s" % (self.tmpdir, self.weightname_h)
	self.imagepath_k = "%s%s" % (self.tmpdir, self.imagename_k)
	self.weightpath_k = "%s%s" % (self.tmpdir, self.weightname_k)

	print "6 Files moved to %s" % self.tmpdir

    def RemoveTmpFiles(self):
        '''Removes all files in the temporary directory after the photometry is done'''
	os.system ("rm -f %smos*" % self.tmpdir)
	print "Removed Temporary Files from %s" % self.tmpdir


def grabphot(filepath,src_ra,src_dec):

    dblock = ptel_astrophot.PTEL_data_block(filepath, req_filts=['j','h','k'])
    dblock.get_phot_at_pos(pos=(src_ra,src_dec),force_photom=False,req_lim=3) #returns dictionary pos_results

    return  dblock.pos_results

def autophot(ra,dec,semester,year,month,day,projectid,objectid,obsid=1):

    ptelimg = ImgFile(semester,year,month,day,projectid,objectid,obsid)
    ptelimg.MoveFiles()
    photdict = grabphot(ptelimg.imagepath_j, ra, dec)
    ptelimg.RemoveTmpFiles()

    print "J Mag = %f +/- %f" % (photdict['j']['mag'], photdict['j']['merr'])
    print "H Mag = %f +/- %f" % (photdict['h']['mag'], photdict['h']['merr'])
    print "K Mag = %f +/- %f" % (photdict['k']['mag'], photdict['k']['merr'])

    return photdict

def grb090313test():

    ra = 198.4007
    dec = 8.097223

    semester = "sem2009a"
    year = 2009
    month = 3
    day = 13
    projectid = "GRB"
    objectid = 388

    obsid = 2 #will default to 1 if not specified

    autophot(ra,dec,semester,year,month,day,projectid,objectid,obsid)
