'''
CoaddWrap.py
Author: Adam N. Morgan

A wrap around Chris's mosaic_maker.py, to create mosaics of every N 
triplestacks created by PAIRITEL Pipeline3.

To run, put this file, along with mosaic_maker.py, anet.py, and 
pairitel_redux.swarp into the folder containing the 
PROJ.OBJ.OBS-reduction_output folder, e.g. GRB.388.1-reduction_output.

This string, "PROJ.OBJ.OBS" ("GRB.388.1") is the "obsid" as defined by this 
program.

The parameter numcoadd specifies how many total observations to coadd in a 
given epoch (say you only want to coadd the first 40 of a 90 observation 
epoch, set numcoadd=40).  By default, numcoadd=None coadds all of them.  

Start python, then do the following:
>>> import CoaddWrap
>>> CoaddWrap.prep("GRB.388.1")
>>> # 3 files (?_long_triplestacks_full.txt) should have been created
>>> CoaddWrap.coadd("GRB.388.1",max_sum=4,dowcs=False)
>>> # This will loop through mosaic_maker, coadding every max_sum images
>>> # They will be renamed ?_long_GRB.388.1_coadd_N-M.fits where N is the 
>>> # First coadded image and M is the last; and M-N = max_sum - 1 (except for
>>> # the last image, which will just be the coaddition of all remaining)
>>> CoaddWrap.cleanup("GRB.388.1")
>>> # Will remove all resultant images and move them to a folder with optional
>>> # naming string: obsid + opt_str + '_mosaics'

'''

import os, sys
import shutil
import glob
pypath = "/Library/Frameworks/Python.framework/Versions/Current/bin/python"

def prep(obsid):
    '''Given a string or list of obsids, combine them into a text file
    containing all of the observations to coadd.
    '''
    # if obsid is a string, convert it into a list
    if isinstance(obsid,str):
        obsid = [obsid]
    # if obsid is not a list now, raise exception
    if not isinstance(obsid,list):
        raise TypeError('obsid is of invalid type; needs to be list or str')
    for oid in obsid:
        globstr = oid + '-reduction_output'
        if not glob.glob(globstr):
            sys.exit('search directory does not exist')
        prepstr = pypath + " mosaic_maker.py -o " + oid + " -p"
        os.system(prepstr)
        globstr = '?_long_triplestacks.txt'
        text_list = glob.glob(globstr)
        if not text_list:
            sys.exit('search string does not exist')
        for item in text_list:
            new_item = item.replace('stacks.txt','stacks_full.txt')
            syscmd = 'cat %s >> %s' % (item,new_item)
            os.system(syscmd)
#            shutil.move(item,new_item)

def cleanup(obsid,opt_str=''):
    dirname = obsid + opt_str + '_mosaics'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    wdirname = dirname + '/weights'
    if not os.path.exists(wdirname):
        os.mkdir(wdirname)
    wdirout = wdirname + '/.'
    dirout = dirname + '/.'
    wglobstr = '?_long_' + obsid + '*weight.fits'
    wcommand = 'mv %s %s' % (wglobstr, wdirout)
    os.system(wcommand)
    globstr = '?_long_' + obsid + '*.fits'
    command = 'mv %s %s' % (globstr, dirout)
    os.system(command)
    command = 'rm ?_long_triplestacks*.txt'
    os.system(command)
    print "Files moved to %s" % dirout
        
    
def coadd(obsid,max_sum=4,dowcs=False,coadd_range=None):
    
    j_filename_new = "j_long_triplestacks.txt"
    j_filename_old = "j_long_triplestacks_full.txt"
    h_filename_new = "h_long_triplestacks.txt"
    h_filename_old = "h_long_triplestacks_full.txt"
    k_filename_new = "k_long_triplestacks.txt"
    k_filename_old = "k_long_triplestacks_full.txt"

    j_file = file(j_filename_old,"r")
    j_list_full = j_file.readlines()
    h_file = file(h_filename_old,"r")
    h_list_full = h_file.readlines()
    k_file = file(k_filename_old,"r")
    k_list_full = k_file.readlines()

    if len(k_list_full) != len(j_list_full) or len(k_list_full) != len(h_list_full):
        print "WARNING - list lengths not identical"
    
    # If the range of observations to coadd is not specified, just coadd 
    # all of them by default.
    if not coadd_range: 
        numiter = len(j_list_full)
        j_list = j_list_full
        k_list = k_list_full
        h_list = h_list_full
    else:
        if not isinstance(coadd_range,tuple) or not len(coadd_range) == 2:
            sys.exit('coadd_range needs to be a tuple of length 2 ')
        i_start = coadd_range[0] - 1
        i_stop = coadd_range[1]
        numiter = i_stop - i_start
        if not isinstance(i_start,int) or not isinstance(i_stop,int):
            sys.exist('coadd_range values need to be of type integer')
        j_list = j_list_full[i_start:i_stop]
        h_list = h_list_full[i_start:i_stop]
        k_list = k_list_full[i_start:i_stop]
    
    ii = 0
    kk = 0
    j_file_new = file(j_filename_new,"w")
    h_file_new = file(h_filename_new,"w")
    k_file_new = file(k_filename_new,"w")

    for item in j_list:
        print item
        if kk == 0:
            start_seg = item.split("-p")[-1].split('.fits')[0]
        ind = j_list.index(item)
        j_file_new.write(item)
        h_file_new.write(h_list[ind])
        k_file_new.write(k_list[ind])
        ii += 1 
        kk += 1
        if kk == max_sum or ii == numiter:
            end_seg = item.split("-p")[-1].split('.fits')[0]
            j_file_new.close()
            h_file_new.close()
            k_file_new.close()
            
            print 'Now Coadding triplestacks ' + start_seg + '-' + end_seg + ".."
           
            coaddstr = pypath + " mosaic_maker.py -o " + obsid           
            if dowcs:
                coaddstr += ' -w'
            os.system(coaddstr)
            
            j_file_new = file(j_filename_new,"w")
            h_file_new = file(h_filename_new,"w")
            k_file_new = file(k_filename_new,"w")
            kk = 0
            
            globstr = '?_long_'+obsid+'_coadd.*fits'
            seg_str = '_'+start_seg+'-'+end_seg    
            replace_str = '_coadd'+seg_str
            text_list = glob.glob(globstr)
            if not text_list:
                sys.exit('search string does not exist')
            for item in text_list:
                new_item = item.replace('_coadd',replace_str)
                shutil.move(item,new_item)
    j_file_new.close()
    h_file_new.close()
    k_file_new.close()
    
    
