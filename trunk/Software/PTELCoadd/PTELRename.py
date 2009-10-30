'''
PTELRename.py
Author: Adam N. Morgan


'''

import os, sys
import shutil
import glob
pypath = "/Library/Frameworks/Python.framework/Versions/Current/bin/python"

def pl3mos2pl2mos(path='./',outpath=None, weightpath='./',weightoutpath=None, move=False):
    '''
    Rename the output files of PAIRITEL Pipeline3 to those of PAIRITEL Pipeline2

    e.g. 

    j_long_GRB.11081.1_coadd.fits -> mosjGRB.11081.1_2009Jul09.fits 
    j_long_GRB.11081.1_coadd.weight.fits -> mosjGRB.11081.1_2009Jul09_w.fits
    
    if move = False, then it will just make a copy of the images.

    If outpaths are not defined, then just place them in the same directory

    '''
    
    # Eventually want to grab the datestr from the header
    datestr = '2999Jan01'
    
    if path:
        globstr = path + '?_long_*.fits'
        globlist = glob.glob(globstr)
    else:
        globlist = []
    if weightpath:
        wglobstr = weightpath + '?_long*.weight.fits'
        wgloblist = glob.glob(wglobstr)
    else:
        wgloblist = []
    
    if path and not outpath:
        outpath = path
    if weightpath and not weightoutpath:
        weightoutpath = weightpath 
    
    for pl3path in globlist:
        if pl3path.find('weight') == -1:
            pl2path = os.path.basename(pl3path)
            # j_long_GRB.11081.1_coadd.fits
            pl2path = pl2path.replace('j_long_','mosj')
            pl2path = pl2path.replace('h_long_','mosh')
            pl2path = pl2path.replace('k_long_','mosk')
            # mosjGRB.11081.1_coadd.fits
            dashdatestr = '-' + datestr
            pl2path = pl2path.replace('_coadd',dashdatestr)
            # mosjGRB.11081.1-2999Jan01.fits
            pl2path = outpath + pl2path
            if move: 
                cmd = 'mv %s %s' % (pl3path,pl2path)
                printstr = 'Moving '
            elif not move:
                cmd = 'cp %s %s' % (pl3path,pl2path)
                printstr = 'Copying '
            else:
                sys.exit('Invalid entry for "move" keyword')
            print "%s %s to %s" % (printstr,pl3path,pl2path)
            os.system(cmd)
        
    for pl3path in wgloblist:
        pl2path = os.path.basename(pl3path)
        # j_long_GRB.11081.1_coadd.weight.fits
        pl2path = pl2path.replace('j_long_','mosj')
        pl2path = pl2path.replace('h_long_','mosh')
        pl2path = pl2path.replace('k_long_','mosk')
        # mosjGRB.11081.1_coadd.fits
        dashdatestr = '-' + datestr + '_w'
        pl2path = pl2path.replace('_coadd.weight',dashdatestr)
        # mosjGRB.11081.1-2999Jan01.fits
        pl2path = weightoutpath + pl2path
        if move: 
            cmd = 'mv %s %s' % (pl3path,pl2path)
            printstr = 'Moving '
        elif not move:
            cmd = 'cp %s %s' % (pl3path,pl2path)
            printstr = 'Copying '
        else:
            sys.exit('Invalid entry for "move" keyword')
        print "%s %s to %s" % (printstr,pl3path,pl2path)
        os.system(cmd)
        