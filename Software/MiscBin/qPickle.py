#!/usr/bin/env python
# encoding: utf-8
"""
qPickle.py
Author: Adam N Morgan
Created: Oct 7, 2009

This is a generalized wrapper around the pickle module, allowing for easier
loading and saving of pickle files.

"""
import cPickle as pickle
import os

def load(pklpath):
    '''Given an input object and an output path, load a pickle file.'''
    if os.path.exists(pklpath):
        storefile=open(pklpath)
        loadedpkl = pickle.load(storefile)
        storefile.close()
        print "Loaded pickle file for this object."
        return loadedpkl
    else:
        print "Pickle file %s does not exist" % pklpath
        return None

def save(input,outpath,clobber=False):
    path_existed = os.path.exists(outpath)
    if path_existed and not clobber:
        print '%s already exists and clobber == False; not overwriting.' % (outpath)
        return
    else:
        storefile = open(outpath,'w')
        pickle.dump(input,storefile)
        storefile.close
        if not path_existed:
            print "No Pickle file existed, so one was created:"
            print outpath
        else:
            print "Overwrote Pickle file:"
            print outpath
            