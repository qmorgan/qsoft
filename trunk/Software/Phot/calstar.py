from Phot import q_phot
import matplotlib
from scipy import array
from Phot import t_mid
import numpy
import os
import pylab
from MiscBin import qPickle
import glob

storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def magplot(reg, filelist, out_pickle=None, ap=None, triggerid = None, globit = False, noerr=False, magrange=None, caliblimit=True):
    
    '''
    Plot magnitudes of calibration stars as a function of time.
    
    Do after the initial coaddition of triplestacks to plot the magnitudes of 
    calibration stars as a function of time for each science image.
    
    Do once for EACH BAND
    
    Requirements: q_phot and t_mid.
    
    
    Arguments: 
        reg: region file of calibration stars to test and plot
        filelist: EITHER: list of files to run photometry on, or a string
                  to do a glob search of files in a directory (if globit=True)
        out_pickle: override the filename of the pickle file to write to. 
                    Will use default naming convention if not specified.
        triggerid: Swift trigger id of the GRB (if applicable; otherwise None)
        globit: if True, do a glob search to get a list of filenames to run 
                photometry on instead of specifying an explicit list.
        testind: index of image in the filelist from which to run the initial 
                photometry on to get the calibration star keywords. make sure
                all your calib stars are present and viewable in this image.
        noerr: if True, do not plot error bars
    '''

    if globit == True:
        globstr1 = str(filelist) + '*coadd*[0-9].fits'
        globstr2 = str(filelist) + '*coadd.fits'
        globlist1 = glob.glob(globstr1)
        globlist2 = glob.glob(globstr2)
        filelist = globlist1 + globlist2
        print 'globit actiavated'
        print filelist
    
    unique_name = (filelist[0].split('_'))[2]
    filt = filelist[0][0]
    
    calib_star_keys = []
    testind = 0
    caldict = {}
    matplotlib.pyplot.clf()
    regpath = reg
    temppath = storepath + 'temp.reg'
    regfile = open(regpath, 'r')
    reglist = regfile.readlines()
    callist = []
    for line in reglist:
        if 'circle' in line:
            callist += [line]
        else:
            pass
    
    colornumber = len(callist)
    n_stars = len(callist)
    
    while((len(calib_star_keys) < len(callist)) and testind < len(filelist)):
        tempreg = open(temppath, 'w')
        tempreg.write('# Region file format: DS9 version 4.1\n')
        secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
        tempreg.write(secondstr)
        tempreg.write('fk5\n')
        #tmp_str = star_reg
        print callist[0]
        # Add a few arcseconds to the first position to make sure we don't use it as a target
        test_ra = float(callist[0].lstrip('circle(').split(',')[0]) + 0.005
        test_dec = float(callist[0].lstrip('circle(').split(',')[1]) + 0.005

        tmp_str = 'circle(%f,%f,4") # width=2 font="helvetica 16 normal"\n' % (test_ra,test_dec)
    
        tempreg.write(tmp_str)
        tempreg.close()

        # Grab the calib stars we will be looping over:
        calregion =  '/calstarregs/' + os.path.basename(reg) 
        print "Using image #%i to get the calib stars; if not all are present \
                in final plot, try a different image" % (testind)
        testimage = filelist[testind]
        photdict = q_phot.photreturn(os.path.basename(reg), testimage, reg=temppath, calregion=calregion, aper=ap, auto_upper=False, caliblimit=caliblimit)
        
        for key in photdict[testimage]['calib_stars'].keys():
            if not key in calib_star_keys:
                calib_star_keys.append(key)
                
        testind += 1
    print 'length of stuff'
    print len(callist)
    print len(calib_star_keys)
    
    for index, ra_str in enumerate(calib_star_keys):
        # if os.path.exists(temppath):
        #     os.remove(temppath)
        datalist = []
        dataerrlist = []
        timelist = []
        timeerrlist = []
        colorstr = str(float((1/colornumber))*float(index + 1))
        colortuple = (colorstr, 0.5, 0)
        starname = 'star'+str(index)


        precal_dict = {}

        for image in filelist:
            print '**************************************'
            print 'Photometry of star' + str(index) 
            print 'doing image ' + image
            calregion =  '/calstarregs/' + os.path.basename(reg) 
            data = q_phot.photreturn(os.path.basename(reg), image, reg=temppath, calregion=calregion, aper=ap, auto_upper=False, caliblimit=caliblimit)
            image_data = data[image]
            if image[0] != filt:
                raise ValueError('Filter for %s does not match the others') % (image)
            if ra_str in image_data['calib_stars']:
                datalist += [image_data['calib_stars'][ra_str]['new_mag']] 
                dataerrlist += [image_data['calib_stars'][ra_str]['new_e_mag']]
                time = float(t_mid.t_mid(image, trigger=triggerid))
                terr = float(t_mid.t_mid(image, trigger=triggerid,delta = True))/2.
                timetuple = (time, terr)
                image_data.update({'t_mid':timetuple})
                timelist += [time]
                timeerrlist += [terr]
                dec_str = str(image_data['calib_stars'][ra_str]['dec'])[0:7]
                parent_label = image
                precal_dict.update({parent_label:image_data})
            else:
                print 'WARNING: CALIB STAR %s IS NOT USABLE FOR THIS IMAGE' % (ra_str)

        datarr = array(datalist)
        daterrarr = array(dataerrlist)
        timarr = array(timelist)
        timerrarr = array(timeerrlist)
        
        if noerr==True:
            pylab.plot(timarr,datarr,'o',label=str((ra_str,dec_str)))
        else:
            pylab.errorbar(timarr,datarr,yerr=daterrarr,fmt='o',label=str((ra_str,dec_str))) #star_pos_str)
        
        caldict.update({ra_str:precal_dict})

        #matplotlib.pyplot.errorbar(timarr, datarr, yerr = daterrarr, label = starname, fmt='k.', color = colortuple) 
    star_stdv = numpy.std(datarr)
    print 'The standard deviation of the calibration stars is: (DISREGARD THIS, THIS STDV IS PROBABLY WRONG)'
    print star_stdv
    plottitle = '%s Calibration Stars Magnitude vs. t_mid' % (filt)
    plotylabel = '%s Magnitude' % (filt)
    matplotlib.pyplot.title(plottitle)
    matplotlib.pyplot.xlabel('Time After Burst (s)')
    matplotlib.pyplot.ylabel(plotylabel)
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlim((ax.get_xlim()[0]),(ax.get_xlim()[1])*1.2)
    matplotlib.pyplot.legend()

    if magrange:
        ax.set_ylim(magrange)
        
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    DPI = F.get_dpi()
    # F.set_size_inches( (DefaultSize[0]*2.5, DefaultSize[1]*2.5) )
    # was getting incresingly larger with multiple runs
    F.set_size_inches((20, 15))
    n_stars_str = str(n_stars)
    if not out_pickle:
        picklepath = storepath + unique_name + '_' + filt + '_' + 'ap' + str(ap) + '_' + n_stars_str + '_cal_stars.data'
    else:
        picklepath = out_pickle
        
    filepath = storepath + unique_name + '_' + filt + '_' + 'ap' + str(ap) + '_' + n_stars_str + '_cal_stars.png'
    #matplotlib.pyplot.savefig(filepath)

    qPickle.save(caldict, picklepath, clobber=True)

    F.savefig(filepath)
    
    return caldict
    
def star_stdv(magplotdict):

    '''Find the standard deviation of the reference stars using the output 
    picklefile from mag_plot.
    
    Load the pickle file saved from magplot, calculate the standard deviation 
    for each star in here.
    '''
    
    stdv_dict = {}

    for star in magplotdict:
        maglist = []
        magerrlist = []
        
        for image in magplotdict[star]:
            #if 'targ_mag' in magplotdict[star][image]:
             #   print 'targ_mag key found in magplotdict, continuing'
             maglist += [magplotdict[star][image]['calib_stars'][star]['new_mag']]
            # magerrlist += [magplotdict[star][image]['calib_stars'][star]['new_e_mag']]
            #else:
             #   print 'No targ_mag key found in magplotdict, continuing'
        star_stdv = numpy.std(maglist)
        star_stdv_dict = {star:star_stdv}
        stdv_dict.update(star_stdv_dict)
    return stdv_dict


def getstar(reg, out_pickle, filename_h, filename_j, filename_k, \
    ap_h=None,ap_j=None,ap_k=None,triggerid=None, calibration_reg=None, caliblimit=False):
    
    '''
    After creating a calibration deep-stack for all images, use this function
    to perform photomotery of all calibration stars in the calibration region 
    file, and outputs a pickle file which is to be used as a photometry 
    dictionary which is to be used to replace 2mass.
    
    Requirements: q_phot and qPickle. 
    
    note:
    The keyword calibration_reg is for the calibration stars used to 
    calibrate these calibration stars.  For now, just leave as 'none'
    '''

    stardict = {}
    stardict_h = {}
    stardict_j = {}
    stardict_k = {}
    if not ap_h:
        ap_h = raw_input('Enter H aperture: ')
    if not ap_j:
        ap_j = raw_input('Enter J aperture: ')
    if not ap_k:
        ap_k = raw_input('Enter K aperture: ')
        
    regpath = storepath + reg
    regfile = open(regpath, 'r')
    reglist = regfile.readlines()
    temppath = storepath + 'temp.reg'
    star_pos_list = []

    ##################################################################
    #This part is actually not needed, but in case we want to get the star's postition...

    for line in reglist:
        if 'circle' in line:
            star_str = line.strip('circle').strip().strip('")').strip('(').split(',')
            ra_str = star_str[0]
            dec_str = star_str[1]
            star_pos = (float(ra_str), float(dec_str))
            star_pos_list += [star_pos]
        else:
            pass
    
    #End uneeded part of uneededness
    ###################################################################


    callist = []
    for line in reglist:
        if 'circle' in line:
            callist += [line]
        else:
            pass
    
    keylist = []
    
    for index, star_reg in enumerate(callist):
        if os.path.exists(temppath):
            os.remove(temppath)
        starname = 'star'+str(index)
        tempreg = open(temppath, 'w')
        tempreg.write('# Region file format: DS9 version 4.1\n')
        secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
        tempreg.write(secondstr)
        tempreg.write('fk5\n')
        tmp_str = star_reg
        tempreg.write(tmp_str)
        tempreg.close()
        
        star_str = star_reg.strip('circle').strip().strip('")').strip('(').split(',')
        ra_str = star_str[0]
        dec_str = star_str[1]
        ra_round = ra_str[0:8]
        dec_round =dec_str[0:7]
        star_pos = (ra_round, dec_round)
        star_pos_str = str(star_pos)

        data_h = q_phot.dophot(filename_h, temppath, calreg=calibration_reg, ap=ap_h, caliblimit=caliblimit)
        parent_label = star_pos_str
        time = float(t_mid.t_mid(filename_h, trigger=triggerid))
        terr = float(t_mid.t_mid(filename_h, trigger=triggerid,delta = True))/2.
        timetuple = (time, terr)
        data_h.update({'t_mid':timetuple})
        this_star_dict_h = {parent_label:data_h}
        stardict_h.update(this_star_dict_h)
        
        keylist.append(parent_label)
        
        data_j = q_phot.dophot(filename_j, temppath, calreg=calibration_reg, ap=ap_j, caliblimit=caliblimit)
        parent_label = star_pos_str
        time = float(t_mid.t_mid(filename_j, trigger=triggerid))
        terr = float(t_mid.t_mid(filename_j, trigger=triggerid,delta = True))/2.
        timetuple = (time, terr)
        data_j.update({'t_mid':timetuple})
        this_star_dict_j = {parent_label:data_j}
        stardict_j.update(this_star_dict_j)

        data_k = q_phot.dophot(filename_k, temppath, calreg=calibration_reg, ap=ap_k, caliblimit=caliblimit)
        parent_label = star_pos_str
        time = float(t_mid.t_mid(filename_k, trigger=triggerid))
        terr = float(t_mid.t_mid(filename_k, trigger=triggerid,delta = True))/2.
        timetuple = (time, terr)
        data_k.update({'t_mid':timetuple})
        this_star_dict_k = {parent_label:data_k}
        stardict_k.update(this_star_dict_k)

    h_dict = {'h':stardict_h}
    j_dict = {'j':stardict_j}
    k_dict = {'k':stardict_k}
    stardict.update(h_dict)
    stardict.update(j_dict)
    stardict.update(k_dict)

    picklepath = storepath + out_pickle + '.data'
    qPickle.save(stardict, picklepath, clobber = True)
    
    print 'Created a dictionary for the following star locations:'
    print keylist
    
    return stardict

def extract(caldict_h, caldict_j, caldict_k, starRA):
    ''' extract specific stars from caldict '''
    from operator import itemgetter

    matplotlib.pyplot.clf()
    star_h = caldict_h[starRA]
    star_j = caldict_j[starRA]
    star_k = caldict_k[starRA]
    h_list = []
    h_err_list = []
    j_list = []
    j_err_list = []
    k_list = []
    k_err_list = []
    timlist_h = []
    terlist_h = []
    timlist_j = []
    terlist_j = []
    timlist_k = []
    terlist_k = []

    #h
    #sorting w.r.t time
    mosaiclist=[]
    for mosaics in star_h:
        mosaiclist += [star_h[mosaics]]
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)
    
    for mosaics in mosaiclist:
        h_list += [mosaics['calib_stars'][starRA]['new_mag']]
        h_err_list += [mosaics['calib_stars'][starRA]['new_e_mag']]
        timlist_h += [mosaics['t_mid'][0]]
        terlist_h += [mosaics['t_mid'][1]]

    #j
    #sorting w.r.t time
    mosaiclist=[]
    for mosaics in star_j:
        mosaiclist += [star_j[mosaics]]
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)
    
    for mosaics in mosaiclist:
        j_list += [mosaics['calib_stars'][starRA]['new_mag']]
        j_err_list += [mosaics['calib_stars'][starRA]['new_e_mag']]
        timlist_j += [mosaics['t_mid'][0]]
        terlist_j += [mosaics['t_mid'][1]]
    
    #k
    #sorting w.r.t time
    mosaiclist=[]
    for mosaics in star_k:
        mosaiclist += [star_k[mosaics]]
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)
    
    for mosaics in mosaiclist:
        k_list += [mosaics['calib_stars'][starRA]['new_mag']]
        k_err_list += [mosaics['calib_stars'][starRA]['new_e_mag']]
        timlist_k += [mosaics['t_mid'][0]]
        terlist_k += [mosaics['t_mid'][1]]

    matplotlib.pyplot.errorbar(timlist_h, h_list, yerr=h_err_list, xerr=terlist_h, marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green', label = 'h')
    matplotlib.pyplot.errorbar(timlist_j, j_list, yerr=j_err_list, xerr=terlist_j, marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', ecolor = 'blue', label = 'j')
    matplotlib.pyplot.errorbar(timlist_k, k_list, yerr=k_err_list, xerr=terlist_k,  \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'k')

    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    #matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    #ax.set_xscale('log')
    matplotlib.pyplot.legend()
   # uniquename = caldict_h.keys()[0].split('_')[2]
    cname = str(starRA)
    matplotlib.pyplot.title(' Calibration Star '+'('+cname+')')
    savepath = storepath + '_' + cname + '_singlecalib.png'
    
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    DPI = F.get_dpi()
    # F.set_size_inches( (DefaultSize[0]*2.5, DefaultSize[1]*2.5) )
    # was getting incresingly larger with multiple runs
    F.set_size_inches((20, 15))

    print 'single calibration plot saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)    
    matplotlib.pyplot.close()

def extractloop(caldict_h, caldict_j, caldict_k):
    ''' extracts all the calibration star plots from caldicts (done for every h, j, k bands) '''

    for star in caldict_h.keys():
        extract(caldict_h, caldict_j, caldict_k, star)
