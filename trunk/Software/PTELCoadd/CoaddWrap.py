'''
CoaddWrap.py
Author: Adam N. Morgan

A wrap around Chris's mosaic_maker.py, to create mosaics of every N 
triplestacks created by PAIRITEL Pipeline3.

To run, put this file, along with mosaic_maker.py, anet.py, and 
pairitel_redux.swarp into the folder containing the 
PROJ.OBJ.OBS-reduction_output folder, e.g. GRB.388.1-reduction_output.

This string, "PROJ.OBJ.OBS" ("GRB.388.1") is the "obsid" as defined by this 
program.  Actually you can define a list of obsids, though note that the 
computer might not have enough memory to swarp them all together.  It was OK
with 1.5 hours of data, but choked on 3 hours.  I may change the program to 
do some intermediate swarping in the future to deal with this problem.

The parameter coadd_range specifies the range of images to coadd together. 
Setting a coadd_range=(1,40) coadds the first to the 40th image together.
By default, coadd_range=None coadds all of them.  

Start python, then do the following:
>>> import CoaddWrap
>>> CoaddWrap.prep(["GRB.388.1","GRB.388.2"])
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
pypath = "~/Programs/epd-6.1-1-rh5-x86/bin/python"


def smartStack(obsidlist):
    '''I need to comment this more.  This is a rough first go at doing 
    smart coaddition  
    '''
    prep(obsidlist)
    firstid = obsidlist[0] # the first id is all you need for coadd() and cleanup()
    
    j_filename_old = "j_long_triplestacks_full.txt"
    j_file = file(j_filename_old,"r")
    j_list_full = j_file.readlines()
    
    total_length = len(j_list_full)
        
    doubling_time = 4
    initial_sum_length = 0
    initial_obs_number = 1
    
    doubling_count=0
    length_count=0
    
    sum_length = initial_sum_length
    obs_num_i = initial_obs_number
    obs_num_f = 0
    
    while obs_num_f < total_length:
                    
        doubling_count += 1
        obs_num_f = obs_num_i + sum_length
        myrange = (obs_num_i, obs_num_f)
        
        
        # if we're running out of observations, tack the rest on to the end
        if obs_num_f + sum_length > total_length:
            obs_num_f = total_length
            myrange = (obs_num_i, obs_num_f)
        
        obs_num_i = obs_num_f + 1
        
        if doubling_count == doubling_time: 
            doubling_count = 0
            if sum_length == 0:
                sum_length = 1
            else:
                sum_length *= 2
        
        
        print myrange 
        coadd(firstid,coadd_range=myrange)

    cleanup(firstid)

def prep(obsid, exclude=False):
    '''Given a string or list of obsids, combine them into a text file
    containing all of the observations to coadd. Exclude keyword will 
    exclude triplestacks which times are specified 
    (eg:['06h10m32s', '06h11m08s', '06h11m44s']).  
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
            # First sort the text file since the new version of the pipeline
            # Doesn't seem to have these lines sorted by default.
            f=open(item,'r')
            linelist=f.readlines()
            linelist.sort()
            new_item = item.replace('stacks.txt','stacks_full.txt')
            f.close()
            f=open(new_item,'w')
            for line in linelist:
                f.write(line)
            f.close()
            # syscmd = 'cat %s >> %s' % (item,new_item)
            # os.system(syscmd)
            syscmd = 'rm %s' % (item)
            os.system(syscmd)
#            shutil.move(item,new_item)
    # erasing lines that are excluded
    if not exclude:
        pass
    else:
        for TStime in exclude:
            globlist = glob.glob('?_long_triplestacks_full.txt')
            for globname in globlist:
                syscmd = 'sed -e \'/%s/d\' %s >> j_temp.txt' % (TStime,globname)
                os.system(syscmd)
                os.remove(globname)
                os.rename('j_temp.txt', globname)


def cleanup(obsid,opt_str=''):
    dirname = obsid + opt_str + '_mosaics'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    wdirname = dirname + '/weights'
    # if not os.path.exists(wdirname):
    #     os.mkdir(wdirname)
    # wdirout = wdirname + '/.'
    dirout = dirname + '/.'
    wglobstr = '?_long_' + obsid + '*weight.fits'
    # wcommand = 'mv %s %s' % (wglobstr, wdirout)
    wcommand = 'mv %s %s' % (wglobstr, dirout)
    os.system(wcommand)
    globstr = '?_long_' + obsid + '*.fits'
    command = 'mv %s %s' % (globstr, dirout)
    os.system(command)
    command = 'rm ?_long_triplestacks*.txt'
    os.system(command)
    moscommand = 'mv *_mosaics.txt ./' + dirname 
    os.system(moscommand)
    print "Files moved to %s" % dirout
        
    
def coadd(obsid,max_sum=None,dowcs=False,coadd_range=None):
    
    fj = open('j_mosaics.txt' , 'w')
    fh = open('h_mosaics.txt' , 'w')
    fk = open('k_mosaics.txt' , 'w')
    
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
    if not max_sum:
        max_sum = len(j_list_full)
    else:
        if not isinstance(max_sum,int) or max_sum < 1:
            raise TypeError('max_sum needs to be a positive integer')
    
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
                print new_item
                if 'weight' in new_item:
                    pass
                else:
                    if 'h_long' in new_item:
                        strname = new_item + '\n'
                        fh.write(strname)
                    elif 'j_long' in new_item:
                        strname = new_item + '\n'
                        fj.write(strname)
                    elif 'k_long' in new_item:
                        strname = new_item + '\n'
                        fk.write(strname)
                shutil.move(item,new_item)
    j_file_new.close()
    h_file_new.close()
    k_file_new.close()
    fj.close()
    fh.close()
    fk.close()
