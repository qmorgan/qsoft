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
import pyfits
pypath = "python"
mosaic_maker_path = '$Q_DIR/trunk/Software/PTELCoadd/mosaic_maker.py'
swarp_bin = "swarp"
sethead_bin = "sethead"

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

# Simple function allowing parallelization of mosaicing with SWarp.
def run_swarp(command):
    os.system(command)
    return
#-------------------------------------------------------------------------------
### ADAM - PUT THE COADD STUFF IN HERE
### RECORD INFORMATION ABOUT CHANGING THE SEXTRACTOR PARAM FILES
### NOTE - TEST HOW CALIB STAR TRENDS CHANGE WITH APERTURE SIZE

def Coadd(filelist, outname):
    try:
        from multiprocessing import Pool
        from multiprocessing import cpu_count as cpuCount
        doparallel = 1
    except:
        print "Parallel processing library not installed. Not paralellizing."
        doparallel = 0
        
    if doparallel == 1:
        numprocessors = cpuCount()
    else:
        numprocessors = 1
        
    from time import time
    # Begin program execution timing.
    start_time = time()
    
    # OBTAIN WORKING DIRECTORIES 
    # reduction_output_directory = str(obs_string) + "-reduction_output"
    # 
    # ## Obtain Obs String
    # globlist = glob.glob(obs_string+'*-reduction_output')
    j_long_list = file("j_long_mosaics.txt", "w")

    j_long_list_weights = file("j_long_mosaic_weights.txt", "w")


    j_stop_list = []
    j_start_list = []

    for item in filelist:
        # Requires filename to be something.fits, and the weight file to be something.weight.fits
        j_path = item
        if j_path.find('.fits') == -1:
            print j_path
            raise ValueError('File must end in .fits')
        if j_path.find('.weight.fits') != -1:
            print j_path
            raise ValueError('Cannot coadd weight files')
        j_w_path = item.split('.fits')[0] + '.weight.fits'
        j_long_list.write(j_path+'\n')
        j_long_list_weights.write(j_w_path+'\n')
    
        # Obtain the start and stop times of the image
        j_hdulist = pyfits.open(j_path)
        j_header = j_hdulist[0].header
        j_stop_list.append(str(j_header["STOP_CPU"]))
        j_start_list.append(str(j_header["STRT_CPU"]))
        j_hdulist.close()
    
    ## Sort the lists of start and stop times
    j_start_list.sort()
    j_stop_list.sort()

    ## Take the first start time and the last stop time
    j_earliest_start = j_start_list[0]
    j_latest_stop = j_stop_list[-1]

    # Close the relevant files
    j_long_list.close()
    j_long_list_weights.close()
    
    if outname.find('.fits') == -1:
        print outname
        raise ValueError('Output file must end in .fits')
    outweightname = outname.split('.fits')[0] + '.weight.fits'
    
    # Run the mosaicing. 
    # Make list of swarp_commands.
    swarp_commands = [
        swarp_bin + " @j_long_mosaics.txt " + 
        "-c " + loadpath + "pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @j_long_mosaic_weights.txt " + 
        "-IMAGEOUT_NAME " + outname +
        " -WEIGHTOUT_NAME "+ outweightname
        ]

    if doparallel == 1:
        # Run the mosaicing with parallel processing.
        p = Pool(numprocessors)
        result = p.map_async(run_swarp, swarp_commands)
        poolresult = result.get()
    else:
        # Run the mosaicing without parallel processing.
        for command in swarp_commands:
            run_swarp(command)

    # We insert the STRT_CPU of the first triplestack and the STOP_CPU of 
    # the last triplestack used to make the mosaic.

    j_hdulist = pyfits.open(outname,mode='update')

    j_header = j_hdulist[0].header

    j_header.update('STRT_CPU',j_earliest_start)
    j_header.update('STOP_CPU',j_latest_stop)
    
    j_hdulist.flush()
    j_hdulist.close()

    # Remove extraneous text files
    #system('rm ?_long*mosaic*.txt')

    # End program execution timing.
    end_time = time()
    total_time = end_time - start_time
    print "Program finished, execution time %f seconds." % total_time

def MakeDeepStack(path='./'):
    '''Coadd all like-filter images together in a particular folder.  Used
    after SmartStackRefine (or some other coadding scheme) to make a deep stack
    of every available image.  Assumes a format of [j,h,k]_*coadd[0-9].fits or
    [j,h,k]_*coadd.fits
    '''
    filt_list = ['j','h','k']
    for filt in filt_list:
        globstr1 = path + filt + '_*coadd*[0-9].fits'
        globstr2 = path + filt + '*coadd.fits'
        globlist1 = glob.glob(globstr1)
        globlist2 = glob.glob(globstr2)
        filelist = globlist1 + globlist2
        outname = filt + '_long_deepstack.fits'
        Coadd(filelist,outname)


def smartStackRefine(obsidlist, path=None, date='', mins2n=20, minfilter='j', \
                exclude=False, regfile='j.reg', wcs=False, mincoadd=1, \
                maxcoadd=150, fibonacci=True):
    '''Here we continually coadd each observation in the obs list until 
    the minimum signal to noise is reached.  q_phot is needed.
    
    obsidlist: list of observation IDs 
    date: string date of the observation; required for the new pipleine3 naming scheme
    mins2n: the minimum signal to noise of the source to coadd to (as determined by sextractor; might be underestimate)
    minfilter: which filter the minimum s2n must correspond to
    regfile: the region file containing the location of the source
    wcs: apply wcs fitting at the end of the coaddition? 
    mincoadd: the minimum number of observations to add together at each refinement step
    fibonacci: if true, increase the number of observations to add together according to the fibonacci sequence
    
    '''
    from Phot import q_phot
    
    ## PHOTOMETRY PARAMTERS ##
    ap=3
    doupper=False
    
    # Choose which filter to base the minimum s/n off of
    # i.e. require the s/n to be above the threshold for a specific filter, 
    # any of them, or all of them, to continue on.
    minfilteroklist = ['j','h','k','all','any']
    if minfilter.lower() not in minfilteroklist:
        raise ValueError
    elif minfilter.lower() == 'j':
        filtindex = 1
    elif minfilter.lower() == 'h':
        filtindex = 0
    elif minfilter.lower() == 'k':
        filtindex = 2
    elif minfilter.lower() == 'any' or minfilter.lower() == 'all':
        print 'Not implemented yet'
        raise ValueError
    ## END PHOTOMETRY PARAM ##
    
    prep(obsidlist, path=path, exclude=exclude, date=date)
    firstid = obsidlist[0] # the first id is all you need for coadd() and cleanup()
    
    # Assume no modifications have been made between j,h, and k files
    # I.e that j file is the same order, etc as the others.
    
    j_filename_old = "j_long_triplestacks_full.txt"
    j_file = file(j_filename_old,"r")
    j_list_full = j_file.readlines()
    
    total_length = len(j_list_full)
    initial_sum_length = 0
    initial_obs_number = 1
    
    doubling_count=0
    length_count=0
    
    sum_length = initial_sum_length
    obs_num_i = initial_obs_number
    obs_num_f = 0
    
    coaddlist = []
    copyidlist = obsidlist[:]
    
    # Raise error if mins2n is not a number or 
    try:
        mins2n = float(mins2n)
    except:
        raise(TypeError)
    if mins2n < 3.0:
        print 'Unrealistically low S2N'
    elif mins2n < 1.0:
        print 'S2N too low - exiting'
        raise(ValueError)
    
    # Two nested while loops.  One keeps going until all observations are 
    # used up.  The other keeps going until the maximum s/n is reached.    
    # Do the loop until mins2n is reached 
    while obs_num_f < total_length:
        s2n = 0
        coadd_increment = mincoadd
                
        while s2n < mins2n:
        # If the s/n was not reached in the previous iteration,
        # Remove one from the copy list and add it to our current list
            # coaddlist.append(copyidlist.pop(0))
            #         
            # print 'added: ' + str(coaddlist)
            # print 'remaining: ' + str(copyidlist)
                        
            # increment the observation number; either mincoadd value or an
            # increasing value due to fibonacci add
            obs_num_f += coadd_increment
            # Don't try to add more observations than we have: 
            if obs_num_f > total_length:
                obs_num_f = total_length+1
            
            myrange = (obs_num_i, obs_num_f)
            num_of_obs = obs_num_f - obs_num_i
            
            dobreak = False
            # if we're running out of observations, tack the rest on to the end
            if obs_num_f + sum_length > total_length:
                obs_num_f = total_length
                myrange = (obs_num_i, obs_num_f)
                print 'We have run out of run out observations'
                dobreak = True 
            
            # Coadd Everything in coaddlist
            new_coadd_list = coadd(firstid,coadd_range=myrange, dowcs=wcs)
            
            # Do photometry on the resultant coadd 
            print 'Now doing photometry on %s' % (new_coadd_list[filtindex])
            photdict = q_phot.dophot(new_coadd_list[filtindex],regfile,ap=ap,do_upper=doupper)
            
            if not 'targ_s2n' in photdict:
                # if an upper limit found, give a tiny s2n so the loop keeps going
                s2n=0.0
            else:
                s2n = photdict['targ_s2n']
            
            print 'The S/N reached is %s from a total of %s images: %s' % (str(s2n), str(obs_num_f - obs_num_i + 1), str(myrange))
            
            # End here if we've run out of observations to stack.  
            if dobreak:
                break
            
            if s2n >= mins2n or num_of_obs >= maxcoadd:
                if s2n >= mins2n:
                    print 'Minimum S/N Reached! '
                if num_of_obs >= maxcoadd:
                    print 'Maximum allowed number of obs reached.'
                print 'Moving on to the next images set.'
                print 'Renaming images.'
                basename = photdict['FileName'].rstrip('fits')[1:]
                base_split = basename.split('coadd')
                newbasename = base_split[0] + str(obs_num_i) + '-' + str(obs_num_f) + '_coadd' + base_split[1]
                filt_list = ['j','h','k']
                for filt in filt_list:
                    mvname = 'mv %s%sfits %s%sfits' % (filt, basename, filt, newbasename)
                    mvwtname = 'mv %s%sweight.fits %s%sweight.fits' % (filt, basename, filt, newbasename)
                    mvcatname = 'mv %s%sfinalcat.txt %s%sfinalcat.txt' % (filt, basename, filt, newbasename)
                    os.system(mvname)
                    os.system(mvwtname)
                    os.system(mvcatname)
            else:
                # If not reached, add one more to the coaddlist and try again.
                basename = photdict['FileName'].rstrip('fits')[1:]
                # remove the ?_long_basename*.fits files created
                rmname = 'rm -f ?' + basename + '*fits'
                rmcatname = 'rm -f ?' + basename + '*finalcat.txt'
                print 'performing system command ' + rmname
                os.system(rmname)
                os.system(rmcatname)
                
                if fibonacci:
                    coadd_increment += obs_num_f - obs_num_i
            
            # Update the initial obs number for the NEXT observation
        obs_num_i = obs_num_f + 1

        print myrange 

    cleanup(firstid)
    

def smartStackDoubling(obsidlist, doubling_time, path=None):
    '''I need to comment this more.  This is a rough first go at doing 
    smart coaddition using a doubling technique. 
    '''
    prep(obsidlist)
    firstid = obsidlist[0] # the first id is all you need for coadd() and cleanup()
    
    # Assume no modifications have been made between j,h, and k files
    # I.e that j file is the same order, etc as the others.
    j_filename_old = "j_long_triplestacks_full.txt"
    j_file = file(j_filename_old,"r")
    j_list_full = j_file.readlines()
    
    total_length = len(j_list_full)
        
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
        
        # Update the initial obs number for the NEXT observation
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

def prep(obsid, date='', path=None, exclude=False, single=False):
    '''Given a string or list of obsids, combine them into a text file
    containing all of the observations to coadd. Exclude keyword will 
    exclude triplestacks which times are specified 
    (eg:['06h10m32s', '06h11m08s', '06h11m44s']).  
    '''
    makername = '/mosaic_maker.py'
    
    # find path:
    if not path:
        path = './'
    outputdir = os.path.abspath(path) + makername
    if not os.path.exists(outputdir):
        copyit = "cp %s %s" % (mosaic_maker_path, outputdir)
        os.system(copyit)
    # first, remove old text files
    rmstr = 'rm ?_long_triplestacks_full.txt'
    os.system(rmstr)
    # if obsid is a string, convert it into a list
    if isinstance(obsid,str):
        obsid = [obsid]
    # if obsid is not a list now, raise exception
    if not isinstance(obsid,list):
        raise TypeError('obsid is of invalid type; needs to be list or str')
    for oid in obsid:
        globstr = path + oid +'_'+ date + '-reduction_output'
        if not glob.glob(globstr):
            printstr = 'Search directory does not exist: %s' % (globstr)
            sys.exit(printstr)
        if not single:
            prepstr = pypath + ' ' + outputdir + " -o " + oid + " -p"
        if date and not single:
            prepstr = pypath + ' ' + outputdir + " -o " + oid + " -d " + date + " -p"
        if single and not date:
            prepstr = pypath + ' ' + outputdir + " -r -o " + oid + " -p"	
        if single and date:
            prepstr = pypath + ' ' + outputdir + " -r -o " + oid + " -d " + date + " -p"
            
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
            sorted_item = item.replace('stacks.txt','stacks_sorted.txt')
            f.close()
            f=open(sorted_item,'w')
            for line in linelist:
                f.write(line)
            f.close()
            new_item = item.replace('stacks.txt','stacks_full.txt')
            syscmd = 'cat %s/%s >> %s/%s' % (path,sorted_item,path,new_item)
            os.system(syscmd)
            syscmd = 'rm %s/%s' % (path,item)
            os.system(syscmd)
            syscmd = 'rm %s/%s' % (path,sorted_item)
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


def cleanup(obsid,path=None,opt_str=''):
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
        
    
def coadd(obsid, path=None, max_sum=None,dowcs=False,coadd_range=None, single=False):
    
    mosaic_list = []
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
            if not single:
                coaddstr = pypath + " mosaic_maker.py -o " + obsid 
            else:
                coaddstr = pypath + " mosaic_maker.py -r -o " + obsid  
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
                    mosaic_list.append(new_item)
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
    
    # Return an array of the new files created
    return mosaic_list
