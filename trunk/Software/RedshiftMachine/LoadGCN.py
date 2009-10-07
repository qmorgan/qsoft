#!/usr/bin/env python
# encoding: utf-8
"""
loadGCN.py
Author: Adam Morgan
Created: Aug 5, 2009
Last Updated: Aug 5, 2009
	Created
	
This was created as a wrapper around ParseGCNNotice.py since apparently you
cannot create and load a pickle file for an instance of a class correctly
from within the program that defines that class.  For more info (though I 
didn't read too thoroughly), see the following website:
http://stefaanlippens.net/pickleproblem

This program loads an instance of GCNNotice based on the trigger ID.  If it
has never been loaded before, it saves the resultant instance of that class
in a pickle file so that the web won't need to be called 
"""
import cPickle as pickle
import os
import sys
from ParseGCNNotice import GCNNotice
from MiscBin import qPickle

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

# May want to set clobber = True if, e.g., new notices have been added since
# we last checked.
def LoadGCN(triggerid, clobber=False):
    ### LOAD or CREATE PICKLE STORAGE FILE 
    # Attempt to load pickle file
    pklpath = storepath+'sw'+str(triggerid)+'GCN.pkl'
    loadedgcn = qPickle.load(pklpath)
    # If couldn't load, or clobber == True, create a new instance of the class
    if clobber or not loadedgcn:
        loadedgcn = GCNNotice(triggerid)
        try:
            loadedgcn.extract_values()
            loadedgcn.get_positions()
            if loadedgcn.successful_load:
                qPickle.save(loadedgcn,pklpath,clobber=True)
            else:
                print 'Could not succesfully load GCN.'
                return
        except:
            print "Could not Extract Values for GCN."
    return loadedgcn
        
    # pklpath = storepath+'sw'+str(triggerid)+'GCN.pkl'
    # if os.path.exists(pklpath) and clobber==False:
    #     loadedgcn = qPickle.load(pklpath)
    #     print "GCNs for trigger %s already downloaded;" % triggerid
    #     print "Loaded pickle file for this trigger."
    # else:
    #     ### Create a pickle file if it doesn't exist, or if clobber is set
    #     loadedgcn = GCNNotice(triggerid)
    #     try:
    #         loadedgcn.extract_values()
    #         loadedgcn.get_positions()
    #     except:
    #         print "Could not Extract Values for GCN."
    #     
    #     if loadedgcn.successful_load==True:
    #         qPickle.save(loadedgcn,pklpath,clobber=True)
    #     else:
    #         print 'Could not succesfully load the GCN!!'
    # return loadedgcn