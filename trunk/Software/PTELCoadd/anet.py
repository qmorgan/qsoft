#!/bin/env python

"""
anet -- command line interface to astrometry.net. Finds WCS solutions to your FITS images
   USAGE: python anet.py {fitsfile.fits | *.fits}
   
"""
import os, urllib, sys, datetime, copy
import urllib2, cookielib
import threading
import mimetools, mimetypes
import os, stat
from cStringIO import StringIO
import pyfits 
import time


class Callable:
    def __init__(self, anycallable):
        self.__call__ = anycallable

# Controls how sequences are uncoded. If true, elements may be given multiple values by
#  assigning a sequence.
doseq = 1

def anet(imlist):
	a = AstrometrySolver()
	a.get_wcs(imlist=imlist)
	print a
	
global b
class MultipartPostHandler(urllib2.BaseHandler):
    handler_order = urllib2.HTTPHandler.handler_order - 10 # needs to run first

    def http_request(self, request):
        data = request.get_data()
        if data is not None and type(data) != str:
            v_files = []
            v_vars = []
            try:
                 for(key, value) in data.items():
                     if type(value) == file:
                         v_files.append((key, value))
                     else:
                         v_vars.append((key, value))
            except TypeError:
                systype, value, traceback = sys.exc_info()
                raise TypeError, "not a valid non-string sequence or mapping object", traceback

            if len(v_files) == 0:
                data = urllib.urlencode(v_vars, doseq)
            else:
                boundary, data = self.multipart_encode(v_vars, v_files)

                contenttype = 'multipart/form-data; boundary=%s' % boundary
                if(request.has_header('Content-Type')
                   and request.get_header('Content-Type').find('multipart/form-data') != 0):
                    print "Replacing %s with %s" % (request.get_header('content-type'), 'multipart/form-data')
                request.add_unredirected_header('Content-Type', contenttype)

            request.add_data(data)
        return request

    def multipart_encode(vars, files, boundary = None, buf = None):
        if boundary is None:
            boundary = mimetools.choose_boundary()
        if buf is None:
            buf = StringIO()
        for(key, value) in vars:
            buf.write('--%s\r\n' % boundary)
            buf.write('Content-Disposition: form-data; name="%s"' % key)
            buf.write('\r\n\r\n' + str(value) + '\r\n')
        for(key, fd) in files:
            file_size = os.fstat(fd.fileno())[stat.ST_SIZE]
            filename = fd.name.split('/')[-1]
            contenttype = mimetypes.guess_type(filename)[0] or 'application/octet-stream'
            buf.write('--%s\r\n' % boundary)
            buf.write('Content-Disposition: form-data; name="%s"; filename="%s"\r\n' % (key, filename))
            buf.write('Content-Type: %s\r\n' % contenttype)
            # buffer += 'Content-Length: %s\r\n' % file_size
            fd.seek(0)
            buf.write('\r\n' + fd.read() + '\r\n')
        buf.write('--' + boundary + '--\r\n\r\n')
        buf = buf.getvalue()
        return boundary, buf
    multipart_encode = Callable(multipart_encode)

    https_request = http_request

opener = urllib2.build_opener(MultipartPostHandler)

results = []

class ImageContainer:
	astrometry_dot_net_url = "http://live.astrometry.net/"
	
	def __init__(self,call_string,name,verbose=True):
		self.reqid       = None
		self.name  = name
		self.call_string = call_string
		self.status      = "unknown"
		self.stat        = "unknown"
		self.newhead = None
		self.verbose = verbose

		self._make_request()
		self._get_req_id()
		self._get_job_status()
		self._get_new_wcs()
		self._replace_new_wcs()
		
		results.append((self.name,self.status,self.stat,self.reqid))
		
	def _make_request(self):
		
		self.time_started= datetime.datetime.now()
		self.status      = "submitted"
		if self.verbose:
			print "** Submitting WCS request for image = %s" % self.name
		self.req  = opener.open(self.astrometry_dot_net_url + "index.php", self.call_string)
		if self.verbose:
			print "   (Finished uploading image = %s)" % self.name
		self.status      = "returned"
	
	def _get_req_id(self):
		if self.status != "returned":
			self.status = "failed"
		tmp = self.req.read().splitlines()
		gotit = False
		for i in range(len(tmp)):
			if tmp[i].find("<title>") != -1:
				gotit = True
				break
		if gotit:
			tmp = tmp[i+1].split("Job ")
			self.reqid = tmp[1].split()[0]
			self.status = "got req id"
		return
	
	def _get_job_status(self,timeout=200.0):
		if self.status != "got req id":
			print "bad job status"
			return
			
		got_status = False
		start = datetime.datetime.now()
		call = self.astrometry_dot_net_url + "status.php?" + urllib.urlencode({"job": self.reqid})
		timeout = datetime.timedelta(seconds=timeout)
		if self.verbose:
			print "   If you'd like to check the status of %s, go to: \n    %s" % (self.name,call)
		while not got_status and datetime.datetime.now() - start < timeout:
			f = urllib.urlopen(call)
			tmp = f.readlines()
			for i in range(len(tmp)):
				if tmp[i].find("<tr><td>Status:</td><td>") != -1:
					self.stat = tmp[i+1].split("</td>")[0]
					if self.stat in ["Failed", "Solved"]:
						got_status = True
						break
			time.sleep(1)
		# print "   Status of file %s (req id = %s).... %s" % (self.name,self.reqid,self.stat)

	def _get_new_wcs(self):
		#print "here1 (%s)" % self.stat
		
		if self.stat != "Solved":
			return
		call = self.astrometry_dot_net_url + "status.php?" + urllib.urlencode({"job": self.reqid, "get": "wcs.fits"})
		self.newhead = "wcs-" + self.reqid + ".fits"
		urllib.urlretrieve(call,self.newhead)
	
	def _replace_new_wcs(self,delnew=True):

		if self.newhead is None or self.stat != "Solved":
			return
		
		if self.name.find(".fits") == -1:
			## not a fits image
			self.status = "wcs=%s" % self.newhead
			return
		wascompressed = False
		if self.name.endswith(".gz"):
			os.system("gunzip " + self.name)
			wascompressed = True
			self.name = self.name.split(".gz")[0]
		tmp = pyfits.open(self.name,"update")
		tmp1 = pyfits.open(self.newhead,"readonly")
		tmp2 = tmp1[0].header
		del tmp2["SIMPLE"]
		del tmp2["BITPIX"]
		del tmp2["NAXIS"]
		
		## copy the header over
		tmp1.close()
		for c in tmp2.ascardlist():
			try:
				tmp[0].header.update(c.key,c.value,c.comment)
			except:
				print 'ERROR updating header for following:'
				print 'Key: ', c.key
				print 'Val: ', c.value
				print 'Com: ', c.comment		
				
		try:
		    tmp.verify("silentfix")
		except:
		    print "CANNOT VERIFY SOMETHING!"
		tmp.close(output_verify='warn')
		if delnew:
			os.remove(self.newhead)
		if wascompressed:
			os.system("gzip " + self.name)
			self.name += ".gz"

		if self.verbose:	
			print "Finished WCS request for image %s (%s)" % (self.name,self.stat)
		
		
class AstrometrySolver:
	
	job_dict = {"uname": "Adam Morgan","email": "amorgan@astro.berkeley.edu", "fsunit" :"arcsecperpix",\
		"fstype-ul": 1, "fsu": 2.2, "fsl": 0.85, "xysrc": "img", "parity": 2, "index": "10arcmin", "tweak": 1,\
		"tweak_order": 2, "imgfile": "","submit": "Submit"}

	def __init__(self,verbose=True):
		self.verbose= verbose
		self.threads = []
		
	def _make_request(self,imgfile=None,pixel_size_range = [0.1,2.5], tweak_astrometry=True):
		
		if imgfile is None or not os.path.isfile(imgfile):
			if self.verbose:
				print "! imgfile is bad"
			return
		tmp =copy.copy(self.job_dict)
		tmp.update({"imgfile": open(imgfile,"rb"), "fsl": pixel_size_range[0], "fsu": pixel_size_range[1], "tweak": int(tweak_astrometry)})
		if imgfile.find(".fits") == -1:
			tmp.update({"index": "auto"})
		self.threads.append(threading.Timer(0.0,ImageContainer,args=[tmp,imgfile],kwargs={'verbose': self.verbose}))
		self.threads[-1].start()

		#print opener.open(self.astrometry_dot_net_url, tmp).read()
		
		#params = urllib.urlencode(tmp)
		#print self.astrometry_dot_net_url + "?" + params
		#f = url
	
	def get_wcs(self,imlist=None,howmany_at_a_time=5):
		
		print "Verbose is set to %s" % repr(self.verbose)
		if imlist is None:
			return
		
		if type(imlist) == type("a"):
			self._make_request(imlist)
			self.threads[-1].join()
			
		if type(imlist) == type([]):
			nsets = len(imlist)/howmany_at_a_time + 1
			for i in range(nsets):
				# print (i, (i*howmany_at_a_time),((1 + i)*howmany_at_a_time))
				for im in imlist[(i*howmany_at_a_time):((i+1)*howmany_at_a_time)]:
					if im.find(".fits") == -1:
						## probably cannot trust the image scale to be small
						self._make_request(im,pixel_size_range=[0.1,500],tweak_astrometry=False)
					else:
						self._make_request(im)
				## wait until the last guy finishes before firing off more
				self.threads[-1].join()
		
		self.threads[-1].join()
	
	def __str__(self):
		a = "RESULTS OF THE SUBMITTED JOBS\n"
		a += "%-45s\t%-10s\t%-10s\t%-15s\n" % ("name","status","stat","reqid")
		a += "*"*100 + "\n"
		for r in results:
			a += "%-45s\t'%-10s'\t%-10s\t%-15s\n" % (os.path.basename(r[0]),r[1],r[2],r[3])
		return a
		#results.append((self.name,self.status,self.stat,self.reqid))


	
if __name__ == "__main__":
	
	if len(sys.argv) <= 1:
		print __doc__
	else:
		a = AstrometrySolver()
		a.get_wcs(imlist=sys.argv[1:])
		print a
	
		
