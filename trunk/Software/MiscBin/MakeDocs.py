#!/usr/bin/env python
# encoding: utf-8
"""
MakeDocs.py
Author: Adam Morgan

This Program uses pydoc to create documentation for 
"""
import sys
import os
import glob

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
docspath = os.environ.get("Q_DIR") + '/docs/'
softpath = os.environ.get("Q_DIR") + '/trunk/Software/'

def MakeDocs():
    modlist = glob.glob(softpath + '*/')
    print modlist

    cmd1 = 'cd %s' % (softpath)
    os.system(cmd1)


    for mod in modlist:
        modname=os.path.dirname(mod).split('/')[-1]
        print modname
    
        cmd2 = 'pydoc -w %s' % (modname)
        os.system(cmd2)
    
        proglist = glob.glob(mod+'/*.py')
        for prog in proglist:
            progname=os.path.basename(prog)
            if progname[0:2] != '__':
                outname = modname+'.'+progname.rstrip('py').rstrip('.')
                cmd2 = 'pydoc -w %s' % (outname)
                os.system(cmd2)
        

    cmd3 = 'mv *.html %s' % (docspath)
    os.system(cmd3)



