import matplotlib
import os
import numpy as np
from numpy import array as arr
import cosmocalc
from Phot import t_mid
from MiscBin import qPickle
import pidly
import glob
import datetime
import pyfits

# autophot, '051008_g.fits', '051008.sdss', 'allhosts.pos', rad=1.2

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'
sextractor_bin = "sex"

idl_path = '/Applications/itt/idl71/bin/idl'
idl = pidly.IDL(idl_path)

class Event():
    """docstring for Event"""
    def __init__(self, eventname):
        self.eventname = eventname
        
class Image():
    """docstring for Image
    
    
    
    """
    def __init__(self, imagefilename,objectfile=None, calfile=None,autocal=True,scope=None,autofilter=True,filt=None):
        self.imagefilename = imagefilename
        
        image_name = self.imagefilename
        # open up the file and read the header
        hdulist = pyfits.open(image_name)
        self.image_data = hdulist[0].data
        self.imagefile_header = hdulist[0].header
        hdulist.close()
        
        # Get scope before filter because some are determinant on scope
        if not scope:
            self.scope = self._get_scope()
            print "Determined scope to be %s" % self.scope
        else:
            print "Scope explicitly specified as %s" % scope
            self.scope = scope
            
        if not filt and autofilter == True:
            self.filt=self._get_filter()
            print 'Setting filter as %s, determined from image' % (self.filt)
        elif filt:
            self.filt=filt
            print 'Forcing to be %s' % (filt)
        else:
            print "could not determine filter; remaining as None"
            self.filt=None
        
            
        # self.objectname = objectname # could potentially look in header for this name
        if autocal and not calfile:
            # grab the 2mass calibration just for kicks; dont automatically use it though
            twomasscalpath = storepath + imagefilename + 'twomass.txt'
            self.get_2mass_calibration(twomasscalpath)
            if self.twomassfield:
                print "2mass Calibration file created"
            
            print "No calibration file specified; attempting to grab SDSS Calibration"
            
            calpath = storepath + imagefilename + 'sdss.txt'
            self.get_sdss_calibration(calpath) 
            if self.sdssfield:
                self.calfile = self.sdsscal
                print "SDSS Calibration file created"
            else:
                self.calfile=None
                print "Warning: No SDSS calibration; Please provide your own calibration field and run again."
                return

            
        elif not calfile:
            print "No calibration file specified and autocal = False; Please provide calibration file and run again"
            return
        else:
            self.calfile = calfile
            
        self.objectfile = objectfile
        # self.ap = ap
    
    
    def _get_filter(self):
        # determine the filter, if possible
        if self.imagefile_header.has_key("FILTERS"):
            filt = str(self.imagefile_header["FILTERS"])
        elif self.imagefile_header.has_key("FILTER"):
            filt = str(self.imagefile_header["FILTER"])
        elif self.imagefile_header.has_key("FILT"):
            filt = str(self.imagefile_header["FILT"])
        else:
            filt = None
        self.header_filter=filt # this is the string of the actual header keyword
        # do auto filter overrides to determine which filter to use for actual calibration
        if self.scope == 'kait':
            if filt == 'clear':
                filt = 'R'
                print 'Reassigning KAIT clear filter as R'
        return filt
        
    def _get_scope(self):
        # determine the telescope, if possible
        if self.imagefile_header.has_key("TELESCOP"):
            TELESCOP = str(self.imagefile_header["TELESCOP"])
            self.header_scope = TELESCOP
        else:
            TELESCOP = None
        
        self.header_scope=TELESCOP
        
        if TELESCOP.strip() == 'K.A.I.T.':
            scope = 'kait'
        elif TELESCOP.find('PAIRITEL') != -1:
            scope = 'pairitel'
        else:
            scope = TELESCOP
            
        return scope
        
    # def extract_objects(self,objectfile):
    #     '''
    #     If given an object file in the format
    #     Objname RA Dec
    #     Extract each line as a position and add to the object dictionary
    #     ignore anything after the 3rd item
    #     '''
    #     f=open(objectfile,'r')
    #     for line in f.readlines():
    #         if line[0] == '#':
    #             pass
    #         else:
    #             linesplit = line.split()
    #         
        
    def do_phot(self,ap,limsigma=3.0,plotcalib=True,offset_calc_type='weighted_mean'):
        '''
        q_phot.do_phot subkeys:
        ['calib_stars',
         'faintest_s2n',
         'HJD_mid',
         'STOP_CPU',x
         'targ_s2n',
         'sex_faintest',
         'N_dither',
         'HJD_start',
         '2mass_abs_avg_dev',
         'STRT_CPU',x
         't_mid',
         'FileName',x
         'filter',x
         'targ_mag',x
         'Aperture',x
         'HJD_stop',
         'targ_flux',
         'zp',x
         'EXPTIME']x
         '''
         
        if not hasattr(self,'objectfile'):
            print "Need to specify object file first, can't do photometry without positions!"
            return

        photdict = {'FileName':self.imagefilename}
        offset_calc_type=offset_calc_type.lower()

        # error checking of individual telescopes
        if self.scope == 'kait':
            #if KAIT, verify that ccdproc has been done
            if not self.imagefile_header.has_key('BIASID'):
                warnmsg = "BIASID doesnt exist for KAIT image"
                print warnmsg
                # raise Exception(warnmsg)
            if not self.imagefile_header.has_key('DARKID'):
                warnmsg = "DARKID doesnt exist for KAIT image"
                print warnmsg
                # raise Exception(warnmsg)
            if not self.imagefile_header.has_key('FLATID'):
                warnmsg = "FLATID doesnt exist for KAIT image"
                print warnmsg
                # raise Exception(warnmsg)


        if self.imagefile_header.has_key("STRT_CPU"): #PAIRITEL
            strt_cpu = str(self.imagefile_header["STRT_CPU"])
            stop_cpu = str(self.imagefile_header["STOP_CPU"])
        elif self.imagefile_header.has_key("DATE-OBS") and self.imagefile_header.has_key("UT"): #KAIT

            bbb = str(self.imagefile_header["DATE-OBS"])
            ccc = str(self.imagefile_header['UT'])
            exp = float(self.imagefile_header["EXPTIME"])
            strt_cpu = bbb[6:10]+'-'+bbb[3:5]+'-'+bbb[0:2] + ' ' + ccc
            # stop_cpu = 'no_stop_cpu' # could just add the exposure time..
            
            start = datetime.datetime.strptime(strt_cpu,"%Y-%m-%d %H:%M:%S")
            # multiply exposure time by a million, convert to integer, treat as microseconds 
            exp_timedelta = datetime.timedelta(0,0,int(exp*1e6)) #assuming exp in seconds
            stop = start + exp_timedelta
            stop_cpu=datetime.datetime.strftime(stop,"%Y-%m-%d %H:%M:%S")
            
            
        else:
            strt_cpu = 'no_strt_cpu'
            stop_cpu = 'no_stop_cpu'


        photdict.update({'STRT_CPU':strt_cpu})
        photdict.update({'STOP_CPU':stop_cpu})    
        photdict.update({'Aperture':ap})


        #Performing Photometry
        if not plotcalib:
            plotcalibidl = '0'
        else:
            plotcalibidl = '1'
            
        IDL_command = "autophot, '" + str(self.imagefilename) + "', '" + str(self.calfile) +\
            "', '"  + str(self.objectfile) + "', rad=" + str(ap)+", filter='"+str(self.filt)+"', limsigma='"+\
            str(limsigma)+"', plotcalib="+plotcalibidl+", offset_calc_type='"+offset_calc_type+"'"
        print IDL_command
        idl(IDL_command)
        # 
        # #Read the filename
        filename = 'tmpphot.txt'
        f=open(filename,'r')
        lines = f.readlines()
        magdict = {} #dictionary of magnitudes; note departure here from q_phot! 
        for line in lines:
            if line[0:4] == "Phot":
                # ['Phot:', '051008', ':', 'g', '=', '24.46', '+/-', '0.06', '(+/-', '0.01)']
                photlist = line.split()
                objstr = photlist[1]
                filtstr = photlist[3] # CHECK IF SAME GIVEN UP EARLEIR 
                equals = photlist[4]
                if equals.strip() == '=':
                    mag = float(photlist[5])
                    magerr = float(photlist[7])
                    syserr = float(photlist[9].strip(')'))
                    # to match qphot targ mag format of tuple of mag and mag err
                    # here add the stat error and sys error in quadrature
                    # not 100% correct since there is some covariance, but roughly right
                    targ_mag = (mag, np.sqrt(magerr**2+syserr**2))
                    magdict.update({objstr:targ_mag})
                    photdict.update({'targ_mag':targ_mag}) # the last item having photometry done will appear here
                    photdict.update({'filter':filtstr})
                    #saving the individual values and errors seperately
                    photdict.update({'mag_value':mag})
                    photdict.update({'mag_err_stat':magerr})
                    photdict.update({'mag_err_sys':syserr})
                    
                elif equals.strip() == '>':
                    print "No detection; grabbing upper limit"
                    magul = float(photlist[5])
                    limsig = float(photlist[6].lstrip('['))
                    assert photlist[7].strip() == 'sig]'
                    assert photlist[8].strip() == '(+/-'
                    syserr = float(photlist[9].strip(')'))
                    print magul
                    print limsig
                    targ_ul = (magul,limsig)
                    magdict.update({objstr:targ_ul})
                    photdict.update({'filter':filtstr}) 
                    #saving the individual values and errors seperately
                    photdict.update({'mag_value':999})
                    photdict.update({'mag_err_stat':999})
                    photdict.update({'mag_err_sys':syserr})
                    
            elif line[0:4] == "Expt":
                explist = line.split()
                exptime = float(explist[1])
                photdict.update({'EXPTIME':exptime})

            elif line[0:18] == 'Measured zeropoint':
                zplist = line.split()
                zp = (float(zplist[-3]),float(zplist[-1]))
                photdict.update({'zp':zp})
            elif line[0:18] == 'Aperture dependant':
                zp_ap_dep_list = line.split()
                zp_ap_dep = float(zp_ap_dep_list[-1])
                photdict.update({'zp_ap_dep':zp_ap_dep})
            elif line[0:20] == 'observed-catalog RMS':
                rmslist= line.split()
                rms_measured = float(rmslist[-4])
                rms_expected = float(rmslist[-2])
                photdict.update({'offset_rms_measured':rms_measured})
                photdict.update({'offset_rms_expected':rms_expected})
            elif line[0:13] == 'Color scatter':
                color_scatter = float(line.split()[-1])
                photdict.update({'color_scatter':color_scatter})
            elif line[0:13] == 'Median offset':
                offset_median = float(line.split()[-1])
                photdict.update({'offset_median':offset_median})
            elif line[0:20] == 'Weighted mean offset':
                offset_weighted_mean = float(line.split()[-1])
                photdict.update({'offset_weighted_mean':offset_weighted_mean})
            elif line[0:20] == 'Error in weighted me':
                offset_weighted_mean_err = float(line.split()[-1])
                photdict.update({'offset_weighted_mean_err':offset_weighted_mean_err})
                
        photdict.update({'magdict':magdict})
        photdict.update({'offset_calc_type':offset_calc_type})
        
        self.imagedict = photdict
        return photdict


    def get_2mass_calibration(self,twomasscal):
        '''Print 2mass calibration text file (if available) to the current 
        directory with filename twomasscal
        '''
        #IDL> printtwomass, 'OBS1_R.fits', outfile="calib.txt"

        #Performing Photometry
        IDL_command = "printcatalog, '" + str(self.imagefilename) + "', '2mass', outfile='" + str(twomasscal) +  "', unc=1"  
        idl(IDL_command)
        self.twomasscal = twomasscal
        if os.path.exists(self.twomasscal):
            self.twomassfield=True
            print "Successfully wrote twomass calibration file to %s" % (self.twomasscal) 
        else:
            self.twomassfield=False
            print "Failed to write twomass calibration file."
            
    def get_sdss_calibration(self,sdsscal):
        '''Print sdss calibration text file (if available) to the current 
        directory with filename sdsscal
        '''
        #IDL> printsdss, 'OBS1_R.fits', outfile="calib.txt"

        #Performing Photometry
        IDL_command = "printsdss, '" + str(self.imagefilename) + "', outfile='" + str(sdsscal) +  "', unc=1"  
        idl(IDL_command)
        self.sdsscal = sdsscal
        if os.path.exists(self.sdsscal):
            self.sdssfield=True
            print "Successfully wrote sdss calibration file to %s" % (self.sdsscal) 
        else:
            self.sdssfield=False
            print "Failed to write sdss calibration file."
    
    
    def p_photreturn(self,outname,ap,limsigma=3.0,plotcalib=True,\
            offset_calc_type='weighted_mean',clobber=False, \
            utburst=None):
        '''
        attempt to build up same structure as the photdict from q_phot
    
        keys: filename
    
    
        '''
        if utburst == None:
            utburst = datetime.datetime(1858, 11, 17) #just use mjd
        
        offset_calc_type=offset_calc_type.lower()
        photdict={}
        newname =  self.imagefilename + '_ap' + str(ap) 
        filepath = storepath + 'phot_'+ outname #outname is the filepath
        # if calregion:
        #     calibration_list = openCalRegion(calregion)
        #     n_calstars = len(calibration_list)
        #     filepath += '_WithCalReg' + str(n_calstars)
        # if stardict:
        #     filepath += '_WithDeepStack'
        filepath += '.data'
        while clobber == False:     # why did i make this a while loop?? lol
            if os.path.isfile(filepath) == True:
                data = qPickle.load(filepath)
                if newname in data:
                    return data
                else:
                    clobber = True
            else:
                clobber = True
    
        while clobber == True:
    
            if os.path.isfile(filepath) == False:
                photdict = {}  
            else:
                #f = file(filepath)
                photdict = qPickle.load(filepath) # This line loads the pickle file, enabling photLoop to work                
    
            # create dictionary for file
            data = self.do_phot(ap=ap,limsigma=limsigma,plotcalib=plotcalib,offset_calc_type=offset_calc_type)
            if not data:
                print "Photometry failed. No data returned."
                return
            #rerun to get upper limit??
            
            label = newname
            # somehow update time here?
            if self.scope == 'pairitel':
                tdict = {'utburst':utburst,'STOP_CPU':data['STOP_CPU'],'STRT_CPU':data['STRT_CPU']}
                time = float(t_mid.t_mid(time_dict=tdict))
                terr = float(t_mid.t_mid(delta = True, time_dict=tdict))/2.
                timetuple = (time, terr)
                data.update({'t_mid':timetuple})
            elif self.scope == 'kait':
                # untested
                tmid = startexp2tmid(utburst,data['STRT_CPU'],data['EXPTIME']) 
                terr = data['EXPTIME']
                timetuple = (time, terr)
                data.update({'t_mid':timetuple})
                
            photdict.update({label:data})
            qPickle.save(photdict, filepath, clobber = True)
            return photdict
    
    
        qPickle.save(photdict, filepath, clobber = True)
    
    
        return photdict


def p_photLoop(outname,ap,objectfile,calfile,clobber=False,forcefilter=False,\
    offset_calc_type='weighted_mean', utburst=None):
    '''Run photreturn on every file in a directory; return a dictionary
    with the keywords as each filename that was observed with photreturn
    '''   
    import glob
    GRBlist = []
    GRBlistwweight = glob.glob('*.fits')
    filepath = storepath + outname + '.data'
    
    if not os.path.exists(calfile):
        raise ValueError("calfile does not exist in specified path")
    if not os.path.exists(objectfile):
        raise ValueError("objectfile does not exist in specified path")
    if clobber == True:
        if os.path.isfile(filepath) == True:
            os.remove(filepath)
    for item in GRBlistwweight:        # Remove the weight images from the list
        if item.find('weight') == -1:
            GRBlist.append(item)
    if not GRBlist:
        raise ValueError('No files to perform photometry on. In right directory?')
    if forcefilter:
        filt=forcefilter
    else:
        filt=None
    for filename in GRBlist:
        img=Image(filename,objectfile=objectfile,calfile=calfile,filt=filt,autofilter=True)
        photdict=img.p_photreturn(outname,ap=ap,plotcalib=False,offset_calc_type=offset_calc_type,clobber=clobber,utburst=utburst) #dont plot for this
    return photdict

def kait_data_check(directory):
    globlist = glob.glob(directory)

    for filestr in globlist:
        if 'fit' in filestr.split('.')[-1]: #confirm that fits file. could be fit, fits, fts.. 
            hdulist = pyfits.open(filestr)
            header = hdulist[0].header
            hdulist.close()
            string = "%s \t %s \t %s \t %s \t %s" % (filestr, header['UT'], header['EXPTIME'], header['FILTERS'], header['AIRMASS'])
            print string


def textoutput(photdict,objname,utburst=None,filt=None, source=None, day=False):
    '''outputs a text file from photdict.  If filt specified, only output that
    particular filter. Requires a utburst string in the form of hh:mm:ss. 
    Days are integers (e.g. second day observation => day=2. WARNING: If the 
    observations consists of multiple days, care needs to be exercised! 
    Current fix: do it manually!
    
    
    Here we provide the optional utburst string to find tmid of the burst 
    if not specified in the dictionary.
    '''

    import datetime
    from operator import itemgetter
    
    # filts = ['j','h','k']
    # if filt:
    #     if not filt in filts:
    #         raise ValueError('Unacceptable value for filt. Needs to be j, h, or k')
            
    # objname = photdict.keys()[0]
    if filt:
        objname = objname + '_' + filt
    savepath = storepath + objname + '_data.txt'
    text = file(savepath, "w")

    text.write('@inunit=days')
    text.write('\n')
    text.write('@expunit=sec')
    text.write('\n')
    if utburst !=  None:
        utburst_text = '@utburst=' + str(utburst)
        text.write(utburst_text)
    text.write('\n')
    scopetext = '@source=%s' % (source)
    text.write(scopetext)
    text.write('\n')
    namelist = ['%', 't_mid', 'ut_start_date', 'ut_start_time', 'exp', 'filt', '=', 'mag', 'emag', 'lim']
    # namelist = ['%', 't_mid', 'ut_start_date', 'ut_start_time', 'ut_end_date', 'ut_end_time', 'exp', 'filt', '=', 'mag', 'emag', 'lim']
    text.write(' '.join(namelist))
    text.write('\n')

    mosaiclist=[]
    for mosaics in photdict:
        mosaiclist += [photdict[mosaics]]
    #sorting w.r.t time
    get = itemgetter('STRT_CPU')
    mosaiclist.sort(key=get)

    h_list = []
    her_list = []
    j_list = []
    jer_list = []
    k_list = []
    ker_list = []
    tstart_list = []
    tstop_list = []
    exp_list = []

    #old - doesnt work with new utburst format
    # burst_time_sec = float(utburst.split(':')[0])*3600 + float(utburst.split(':')[1])*60 + float(utburst.split(':')[2]) 
    
    for mosaics in mosaiclist:
        exp = mosaics['EXPTIME']
        
        strt_cpu = mosaics['STRT_CPU']
        # kind of a hack; i tried to fix this for kait and auto-assign tmid in photreturn
        if 't_mid' not in mosaics:
            if utburst != None:
                print "ATTEMPTING TO AUTO-EXTRACT T-MID FROM EXPOSURE TIME AND TSTART"
                tmid = startexp2tmid(utburst,strt_cpu,exp)
            else:
                print 't_mid not in mosaics and utburst not specified; cannot get t_mid'
                raise Exception
        # stop_cpu = mosaics['STOP_CPU']
        t_mid_days = tmid/86400.
        
        timestart_str = mosaics['STRT_CPU'].split(' ')[1]
        timestart_sec = float(timestart_str.split(':')[0])*3600. + float(timestart_str.split(':')[1])*60. + float(timestart_str.split(':')[2])
        # if day:
        #     start_after_burst_sec = timestart_sec - burst_time_sec + (day-1)*86400
        # else:
        #     start_after_burst_sec = timestart_sec - burst_time_sec
                # timestop_str = mosaics['STOP_CPU'].split(' ')[1]
        # timestop_sec = float(timestop_str.split(':')[0])*3600. + float(timestop_str.split(':')[1])*60. + float(timestop_str.split(':')[2])
        # if day:
        #     stop_after_burst_sec = timestop_sec - burst_time_sec + (day-1)*86400
        # else:
        #     stop_after_burst_sec = timestop_sec - burst_time_sec
        
        if 'magdict' in mosaics: 
            mag = float(mosaics['magdict'][objname][0])
            magerr = float(mosaics['magdict'][objname][1])        
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)    
 
        filt = mosaics['filter']
        
        
        if str(magerr)=="3.0":  #indicative of upperlimit
            magerr=0
            # datalist = [str(t_mid_days), str(strt_cpu), str(stop_cpu), str(exp), filt, '=', str(mag), str(magerr), 'yes']
            datalist = [str(t_mid_days), str(strt_cpu),  str(exp), filt, '=', str(mag), str(magerr), 'yes']

        else:   
            # datalist = [str(t_mid_days), str(strt_cpu), str(stop_cpu), str(exp), filt, '=', str(mag), str(magerr)]
            datalist = [str(t_mid_days), str(strt_cpu), str(exp), filt, '=', str(mag), str(magerr)]

        # old    
        # if 'upper_green' in mosaics:
        #     datalist = [str(start_after_burst_sec), str(stop_after_burst_sec), str(exp), filt, '=', str(mag), str(magerr), 'yes']
        # else:   
        #     datalist = [str(start_after_burst_sec), str(stop_after_burst_sec), str(exp), filt, '=', str(mag), str(magerr)]
        # 

        text.write(' '.join(datalist))
        text.write('\n')
    text.close()


def startexp2tmid(utburst,uttstart,exp):
    '''
     utburst = '2012-01-19 04:04:30.21' for 120119A
    exposure time in seconds
    %Y-%m-%d %H:%M:%S
    '''
    if utburst == None:
        utburst = datetime.datetime(1858, 11, 17) #just use mjd
    
    start = datetime.datetime.strptime(uttstart.split('.')[0], "%Y-%m-%d %H:%M:%S")
    #handle the fractions of a second, if there are any
    if len(uttstart.split('.')) == 2:
        start_microseconds_str = uttstart.split('.')[1]
        start_microseconds = int(start_microseconds_str)
        start = start + datetime.timedelta(microseconds=start_microseconds)
    
    burst = datetime.datetime.strptime(utburst.split('.')[0], "%Y-%m-%d %H:%M:%S")
    #handle the fractions of a second, if there are any
    if len(utburst.split('.')) == 2:
        burst_microseconds_str = utburst.split('.')[1]
        burst_microseconds = int(burst_microseconds_str)
        burst = burst + datetime.timedelta(microseconds=burst_microseconds)
        
    tmid = (start + (datetime.timedelta(seconds=exp))/2 - burst)
    return tmid.days*86400 + tmid.seconds + tmid.microseconds/1e6
    
def qmorgan_test_photometry():
    '''
    For quick testing/example of photometry. Not transferrable.
    '''
    print "Initializing directory and filenames..."
    directory = "/Users/amorgan/testphot/"   
    filename = directory + '051008_g.fits'
    position_file = directory + 'allhosts.pos'
            
