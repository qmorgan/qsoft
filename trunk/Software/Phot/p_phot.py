import matplotlib
import os
import numpy as np
from numpy import array as arr
import cosmocalc
from Phot import t_mid
from MiscBin import qPickle
import pidly
import glob
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
    """docstring for Image"""
    def __init__(self, imagefilename,objectfile=None, calfile=None,autocal=True,scope=None,autofilter=True,filt=None):
        self.imagefilename = imagefilename
        
        image_name = self.imagefilename
        # open up the file and read the header
        hdulist = pyfits.open(image_name)
        self.image_data = hdulist[0].data
        self.imagefile_header = hdulist[0].header
        hdulist.close()
        
        if not filt and autofilter == True:
            self.filt=self._get_filter()
        elif filt:
            self.filt=filt
        else:
            print "could not determine filter; remaining as None"
            self.filt=None
            
        if not scope:
            self.scope = self._get_scope()
            print "Determined scope to be %s" % self.scope
        else:
            print "Scope explicitly specified as %s" % scope
            self.scope = scope
            
        # self.objectname = objectname # could potentially look in header for this name
        if autocal and not calfile:
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
        
    def do_phot(self,ap,limsigma=3.0):
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
            strt_cpu = bbb[6:10]+'-'+bbb[3:5]+'-'+bbb[0:2] + ccc
            stop_cpu = 'no_stop_cpu' # could just add the exposure time..
        else:
            strt_cpu = 'no_strt_cpu'
            stop_cpu = 'no_stop_cpu'


        photdict.update({'STRT_CPU':strt_cpu})
        photdict.update({'STOP_CPU':stop_cpu})    
        photdict.update({'Aperture':ap})


        #Performing Photometry
        IDL_command = "autophot, '" + str(self.imagefilename) + "', '" + str(self.calfile) +\
            "', '"  + str(self.objectfile) + "', rad=" + str(ap)+", filter='"+str(self.filt)+"', limsigma='"+\
            str(limsigma)+"'"
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
                    photdict.update({'filter':filtstr})
                    
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
                    photdict.update({'filter':filtstr}) #why again? 
            elif line[0:4] == "Expt":
                explist = line.split()
                exptime = float(explist[1])
                photdict.update({'EXPTIME':exptime})

            elif line[0:18] == 'Measured zeropoint':
                zplist = line.split()
                zp = (float(zplist[-3]),float(zplist[-1]))
                photdict.update({'zp':zp})
        photdict.update({'magdict':magdict})

        self.imagedict = photdict
        return photdict


    def get_sdss_calibration(self,sdsscal):
        #IDL> printsdss, 'OBS1_R.fits', outfile="calib.txt"

        #Performing Photometry
        IDL_command = "printsdss, '" + str(self.imagefilename) + "', outfile='" + str(sdsscal) +  "', unc=1"  
        idl(IDL_command)
        self.sdsscal = sdsscal
        if os.path.exists(self.sdsscal):
            self.sdssfield=True
        else:
            self.sdssfield=False
    
    
    def p_photreturn(self,ap,clobber=False):
        '''
        attempt to build up same structure as the photdict from q_phot

        keys: filename


        '''
        photdict={}

        filepath = storepath + 'phot_'+ self.imagefilename + 'ap' + str(ap)  
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
                if self.imagefilename in data:
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
            data = self.do_phot(ap=ap)
            if not data:
                print "Photometry failed. No data returned."
                return
            #rerun to get upper limit??
            
            label = data['FileName'] + "ap" + str(ap)
            # somehow update time here?
            # time = float(t_mid.t_mid(filename, trigger = trigger_id))
            # terr = float(t_mid.t_mid(filename, delta = True, trigger = trigger_id))/2.
            # timetuple = (time, terr)
            # data.update({'t_mid':timetuple})
            photdict.update({label:data})
            qPickle.save(photdict, filepath, clobber = True)
            return photdict


        qPickle.save(photdict, filepath, clobber = True)


        return photdict


# def p_phot_loop(directory):
#     for photfile in directiory:
#         if 

def kait_data_check(directory):
    globlist = glob.glob(directory)

    for filestr in globlist:
        if 'fit' in filestr.split('.')[-1]: #confirm that fits file. could be fit, fits, fts.. 
            hdulist = pyfits.open(filestr)
            header = hdulist[0].header
            hdulist.close()
            string = "%s \t %s \t %s \t %s \t %s" % (filestr, header['UT'], header['EXPTIME'], header['FILTERS'], header['AIRMASS'])
            print string





def qmorgan_test_photometry():
    '''
    For quick testing/example of photometry. Not transferrable.
    '''
    print "Initializing directory and filenames..."
    directory = "/Users/amorgan/testphot/"   
    filename = directory + '051008_g.fits'
    position_file = directory + 'allhosts.pos'
            
