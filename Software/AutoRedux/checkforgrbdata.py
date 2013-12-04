#! /usr/bin/env python
# encoding: utf-8
"""
checkforgrbdata.py
Checking for the existence of GRB data in PAIRITEL directories
QMorgan June 10 2009 UNTESTED
"""
 
import time
import os
import sys
import glob
import shutil
import cPickle as pickle 

req_num_of_files = 5

### DETERMINE OR CREATE TEMPORARY AUTOREDUX DIRECTORY
autoreduxdir = os.path.expanduser("~/ptel_auto_redux/")
if not os.path.exists(autoreduxdir):  #create storage directory if it doesn't exist already
    try:
        os.system("mkdir " + autoreduxdir)
        os.system("mkdir " + autoreduxdir + "tmp")
        print "Creating Temporary directory %s and moving files there." % (autoreduxdir)
    except:
        print "CANNOT MAKE TEMP DIRECTORY.  Will try to continue anyway, but will likely fail."
        pass
else: 
    print "Copying temporary files to %s" % (autoreduxdir)

pklflagpath = autoreduxdir+'ptelgrbcheck.pkl'
reduxdir = autoreduxdir+'/tmp/'

bodystring = ''
errmsg = ''

### DETERMINE THE SEARCH PATH ###
### The place to search for files is based on the current date
searchstring = "pjrrGRB*.fits"
searchpath = determinesearchpath()

### If the search path exists, coninue.  If not, end program.
if os.path.exists(searchpath):

    ### LOAD or CREATE PICKLE FLAG FILE 
    # The pickle file here contains flags about what has been reduced, etc
    if os.path.exists(pklflagpath): 
        storefile=open(pklflagpath)
        pklflagdict = pickle.load(storefile)
        storefile.close()
    else:
        ### Create an empty pickle file if it doesn't exist
        storefile = open(pklflagpath,'w')
        thedict = {'flag':0,'dateflag':0}
        pklflagdict=thedict
        pickle.dump(pklflagdict,storefile)
        storefile.close
        bodystring += "No Pickle file existed, so one was created\n\n"

        ### DETERMINE NUMBER OF FILES IN DIRECTORY
        globfiles = glob.glob(searchpath+searchstring)
        numfiles = len(globfiles)
    
        ### START COMPOSITION OF EMAIL
        bodystring +="There are now %s files in %s" % (str(numfiles), searchpath)
        bodystring +='\n'

        ### DETERMINE IF CRITERIA TO TAKE ACTION ARE SATISFIED
        ### Need X number of files, and for this to be a different date than before 
        if (numfiles > req_num_of_files) and (pklflagdict['dateflag'] != dateint): 

        	### MOVE FILES TO TEMPORARY PATH FOR REDUCTION
        	### Since files can't be reduced in their initial directory,
        	### Move to a location where the user has write permission 
        	if os.path.exists(reduxdir): 
        	    for i in range(0,numfiles):
                    shutil.copy(globfiles[i],reduxdir)
                    bodystring +=  "\nCopied files to %s" % (reduxdir)
        	else: bodystring += "\nCould not move files, the path %s does not exist." % (reduxdir)

    	### SEND EMAIL TO APPROPRIATE PEOPLE
    	try:  # Only do if the domail module exists
    	    import domail
    	    domail.domail('qmorgan@gmail.com','Files Updated',bodystring)  
    	except:
    	    print bodystring
    	    print "Could not send email; no domail.py module found". 

    	### UPDATE DICTIONARY FLAGS
    	pklflagdict['dateflag'] = dateint #Only perform actions once per date
	
        # Is another flag useful for something?
        # Maybe have it be a init_is_moved or something, but need new pkl files for each 
        if (numfiles > 4):   
    	    pklflagdict['flag'] = 1
        else: pklflagdict['flag'] = 0

        ### SAVE UPDATED INFO IN PICKLE FILE
        storefile = open(pklflagpath,'w')
        pickle.dump(pklflagdict,storefile)
        storefile.close

    ### END

### This is just a definition to reset the flags in the flag dictionary for testing purposes 
def resetflags():
    if os.path.exists(pklflagpath):
        storefile = open(pklflagpath,'w')
	    pklflagdict = {'flag':0,'dateflag':0}
	    pickle.dump(pklflagdict,storefile)
	    storefile.close

def determinesearchpath(sem="sem2009a",year=None,month=None,day=None):

    if (year == None) or (month == None) or (day == None):
	    tlist=time.gmtime()
	    year = tlist[0] ; month = tlist[1] ; day = tlist[2] 
	print "Using current UT date as directory search path"
        sem = "sem2009a"
        datestr = time.strftime("%Y-%b-%d",(year,month,day,0,0,0,0,1,-1))
        dirpath = "/Bloom/PAIRITEL-DATA/%s/Dir%s/" % (sem,datestr)
    return dirpath
