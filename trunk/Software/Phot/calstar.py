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

def magplot(reg, filelist, out_pickle, ap=None, triggerid = None, globit = False):
    
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
        out_pickle: filename of the pickle file to write to
        triggerid: Swift trigger id of the GRB (if applicable; otherwise None)
        globit: if True, do a glob search to get a list of filenames to run 
                photometry on instead of specifying an explicit list.
    '''

    if globit == True:
        globstr1 = str(filelist) + '*coadd*[0-9].fits'
        globstr2 = str(filelist) + '*coadd.fits'
        globlist1 = glob.glob(globstr1)
        globlist2 = glob.glob(globstr2)
        filelist = globlist1 + globlist2
        print 'globit actiavated'
        print filelist

    caldict = {}
    matplotlib.pyplot.clf()
    regpath = reg
    temppath = storepath + 'temp.reg'
    picklepath = storepath + out_pickle +'.data'
    regfile = open(regpath, 'r')
    reglist = regfile.readlines()
    callist = []
    for line in reglist:
        if 'circle' in line:
            callist += [line]
        else:
            pass
    
    colornumber = len(callist)
    
    for index, star_reg in enumerate(callist):
        if os.path.exists(temppath):
            os.remove(temppath)
        datalist = []
        dataerrlist = []
        timelist = []
        timeerrlist = []
        colorstr = str(float((1/colornumber))*float(index + 1))
        colortuple = (colorstr, 0.5, 0)
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

        star_str = star_reg.strip('circle').split(',')
        ra_str = star_str[0].strip('(')
        dec_str = star_str[1]
        ra_round = ra_str[0:8]
        dec_round =dec_str[0:7]
        star_pos = (ra_round, dec_round)
        star_pos_str = str(star_pos)

        precal_dict = {}

        for image in filelist:
            print '**************************************'
            print 'Photometry of star' + str(index) 
            print 'doing image ' + image
            data = q_phot.dophot(image, temppath, ap=ap)
            if 'targ_mag' in data:
                datalist += [data['targ_mag'][0]] 
                dataerrlist += [data['targ_mag'][1]]
                time = float(t_mid.t_mid(image, trigger=triggerid))
                terr = float(t_mid.t_mid(image, trigger=triggerid,delta = True))/2.
                timetuple = (time, terr)
                data.update({'t_mid':timetuple})
                timelist += [time]
                timeerrlist += [terr]
                parent_label = image
                precal_dict.update({parent_label:data})
            else:
                pass

        datarr = array(datalist)
        daterrarr = array(dataerrlist)
        timarr = array(timelist)
        timerrarr = array(timeerrlist)
        
        pylab.errorbar(timarr,datarr,yerr=daterrarr,fmt='o',label=star_pos_str)
        
        caldict.update({star_pos_str:precal_dict})

        #matplotlib.pyplot.errorbar(timarr, datarr, yerr = daterrarr, label = starname, fmt='k.', color = colortuple) 
         
    matplotlib.pyplot.title('Calibration Stars Magnitude vs. t_mid')
    matplotlib.pyplot.xlabel('Time After Burst (s)')
    matplotlib.pyplot.ylabel('Magnitude')
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlim((ax.get_xlim()[0]),(ax.get_xlim()[1])*1.2)
    matplotlib.pyplot.legend()
    
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    DPI = F.get_dpi()
    # F.set_size_inches( (DefaultSize[0]*2.5, DefaultSize[1]*2.5) )
    # was getting incresingly larger with multiple runs
    F.set_size_inches((20, 15))

    unique_name = (filelist[0].split('_'))[2]
    filepath = storepath + unique_name + '_calibration_stars.png'
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
            maglist += [magplotdict[star][image]['targ_mag'][0]]
            magerrlist += [magplotdict[star][image]['targ_mag'][1]]
        
        star_stdv = numpy.std(maglist)
        star_stdv_dict = {star:star_stdv}
        stdv_dict.update(star_stdv_dict)
    return stdv_dict


def getstar(reg, out_pickle, filename_h, filename_j, filename_k, \
    ap_h=None,ap_j=None,ap_k=None,triggerid=None, calibration_reg=None):
    
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

        data_h = q_phot.dophot(filename_h, temppath, calreg=calibration_reg, ap=ap_h)
        parent_label = star_pos_str
        time = float(t_mid.t_mid(filename_h, trigger=triggerid))
        terr = float(t_mid.t_mid(filename_h, trigger=triggerid,delta = True))/2.
        timetuple = (time, terr)
        data_h.update({'t_mid':timetuple})
        this_star_dict_h = {parent_label:data_h}
        stardict_h.update(this_star_dict_h)
        
        keylist.append(parent_label)
        
        data_j = q_phot.dophot(filename_j, temppath, calreg=calibration_reg, ap=ap_j)
        parent_label = star_pos_str
        time = float(t_mid.t_mid(filename_j, trigger=triggerid))
        terr = float(t_mid.t_mid(filename_j, trigger=triggerid,delta = True))/2.
        timetuple = (time, terr)
        data_j.update({'t_mid':timetuple})
        this_star_dict_j = {parent_label:data_j}
        stardict_j.update(this_star_dict_j)

        data_k = q_phot.dophot(filename_k, temppath, calreg=calibration_reg, ap=ap_k)
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



    
