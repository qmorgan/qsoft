#!/usr/bin/env python
# encoding: utf-8
"""
ParseGCNNotice.py
Author: Adam Morgan
Created: July 29, 2009
Last Updated: July 29, 2009
	Created
	
This Program parses the online version of the GCN Notices 
posted online (http://gcn.gsfc.nasa.gov/gcn/swift_grbs.html
for a list) and returns a dictionary of relevant parameters.
See, e.g. http://gcn.gsfc.nasa.gov/gcn/other/358422.swift 
for an example list of GCN notices for a particular trigger.
"""
import sys

class GCNNotice(object):
	"""docstring for GCNNotice"""
	def __init__(self, gcnstring):
		super(GCNNotice, self).__init__()
		self.dict={}
		commentstring = ''
		# for gcn in gcn_list:
		# 	# Search through each notice to determine type
		# 	gcnsplitlist = gcn.splitlines()
		gcnlines = gcnstring.splitlines()
		for line in gcnlines:
			linelist = line.split(':  ')
			if len(linelist) > 2:
				print 'SOMETHINGS WRONG'
			if len(linelist) == 2:
				# Add to dictionary
				if linelist[0] != 'COMMENTS':
					subdict = {linelist[0]:linelist[1].strip()}
					self.dict.update(subdict)
				if linelist[0] == 'COMMENTS':
					commentstring += linelist[1].strip() + ';'
		subdict = {'COMMENTS':commentstring}
		self.dict.update(subdict)
		print "Initialized GCN Notice:", self.dict['NOTICE_TYPE']
			

def where(a,val,wherenot=False):
	"""
	Analogous to the 'where' function in IDL
	See thread:
	http://mail.python.org/pipermail/python-list/2007-March/431322.html
	
	Returns a list of indices in the list 'a' where a[i] == val
	If wherenot=True, returns an index where a[i]!= val
	e.g. 
	>>> a=[0,1,2,5,42,2]
	>>> where(a,42)
	[4]
	>>> where(a,2)
	[2,5]
	>>> where(a,2,wherenot=True)
	[0,1,3,4]
	>>> where(a,999)
	[]
	"""
	if wherenot == False:
		return [i for i in xrange(len(a)) if a[i]==val]
	else:
		return [i for i in xrange(len(a)) if a[i]!=val]

def grabwebgcn(triggerid=358422):
	"""
	Based on trigger number, goes to the web and grabs
	the GCN information and puts it in a list of strings.
	"""
	import urllib2
	
	gcnaddress = 'http://gcn.gsfc.nasa.gov/gcn/other/%s.swift' % str(triggerid)
	gcnwebsite = urllib2.urlopen(gcnaddress)
	gcnstring = gcnwebsite.read()
	thedelimiter = '//////////////////////////////////////////////////////////////////////'
	# Split up the text file by GCN Notice - make it a list of strings
	gcn_notices = gcnstring.split(thedelimiter)
	num_of_gcns = len(gcn_notices)
	# Only keep the GCN notices we care about
	# First search for which GCN notices have a certain string, 
	# THEN parse the latest one using the "where" function defined above
	# This avoids having to parse each case where it's true and avoids
	# parsing each item found in the for loop
	# BATPositionilist=[]
	# BATPositionstring='Swift-BAT GRB Position'
	# FOMWillObserveilist=[]
	# FOMWillObservestring='Swift-FOM Will_Observe'
	# BATLightcurveilist=[]
	# BATLightcurvestring='Swift-BAT GRB Lightcurve'
	# 
	# for gcn in gcn_notices:
	# 	# Go through the list and Pick out only the most recent
	# 	# GRB Parameters from GCN Notices of the same type.
	# 	# Get a list of all GCNs that have the search string:
	# 	BATPositionilist.append(gcn.find(BATPositionstring))
	# 	FOMWillObserveilist.append(gcn.find(FOMWillObservestring))
	# 	BATLightcurveilist.append(gcn.find(BATLightcurvestring))
	# # Grab only the latest notice for a particular type - this assumes 
	# # the last index is the most recent
	
	# SEARCH THROUGH EACH NOTICE
	# FIND WHAT THE NOTICE TYPE IS - 4th line
	gcn_type_list = []
	for gcn in gcn_notices:
		if len(gcn) > 3:
			gcnsplit = gcn.splitlines()
			# Find what the notice type is
			typeline = gcnsplit[3]
			typesplit = typeline.split(':     ')
			if typesplit[0] != 'NOTICE_TYPE':
				print 'THIRD LINE IS NOT NOTICE_TYPE'  
				sys.exit()
			gcn_type_list.append(typesplit[1])
	print gcn_type_list
	# DETERMINE WHAT THE LATEST OF THAT TYPE IS
	
	BATPosition = GCNNotice(gcn_notices[1])
	print BATPosition.dict['GRB_INTEN']
	print BATPosition.dict['COMMENTS']
	# BATPositionindex = where(BATPositionilist,-1,wherenot=True)[-1]
	
		

	
	
	
	

