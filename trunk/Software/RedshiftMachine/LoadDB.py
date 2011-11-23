#!/usr/bin/env python
# encoding: utf-8
"""
LoadDB.py
Author: Adam N Morgan

Load and save qDatabases
"""

import sys
import os
import copy
from MiscBin import qPickle
from MiscBin import qErr

#from MiscBin.q import Standardize

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'


def LoadDB(name, clobber=False, redownload_gcn=False,incl_reg=True,incl_fc=False):
    ### LOAD or CREATE PICKLE STORAGE FILE 
    # Attempt to load pickle file
    pklpath = storepath+'DB_'+str(name)+'.pkl'
    loadeddb = qPickle.load(pklpath)
    # If couldn't load, or clobber == True, create a new instance of the class
    if clobber or not loadeddb:
        loadeddb = None
        # Create new instance of db Notice
        # loadeddb = GRBdb(name,redownload_gcn=redownload_gcn,incl_reg=incl_reg,incl_fc=incl_fc)
        # try:
        #     if loadeddb.successful_load:
        #         # Save new Pickle file
        #         qPickle.save(loadeddb,pklpath,clobber=True)
        #     else:
        #         print 'Could not succesfully load db.'
        #         return
        # except:
        #     print "Could not Extract Values for db."
        #     qErr.qErr()
    return loadeddb

def SaveDB(loadeddb):
    # If the attributes of the loaded GCN have changed since loading, use 
    # this to save the new version of the GCN.
    pklpath = storepath+'DB_'+str(getattr(loadeddb,'name'))+'.pkl'
    qPickle.save(loadeddb,pklpath,clobber=True)