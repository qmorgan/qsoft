#!/usr/bin/env python
# encoding: utf-8
"""
q.py
Author: Adam Morgan
Created: Aug 2, 2009
Last Updated: Aug 2, 2009
	Created
	Adding comment
A conglomeration of random definitions that I use!.
"""
import sys
import os
import math
import numpy
import ephem
import time

# Constants - all CGS
c = 2.99792458e10 # cm/s
h = 6.626e-27 # erg s

#import units

class Obs:
    """Observatory Class"""
    def __init__(self,name='Empty Observatory',scopes={},
    latitude='0:00:00.0',longitude='0:00:00.0',elevation=0.0):
        self.name=name
        self.scopes = scopes
        self.Observer = ephem.Observer()
        self.Observer.lat = latitude
        self.Observer.long = longitude
        self.Observer.elevation=elevation
    
    def lst(self):
        '''Given the location of the observatory, calculate the local sidereal
        time.
        '''
        self.Observer.date=ephem.now()
        return str(self.Observer.sidereal_time())
    
    def lstclock(self):
        print "Ctrl-C to exit clock"
        print ">>>>>><<<<<<"
        while True:
            timestr = ">>%s<<\r" % str(self.lst())[:8]
            sys.stdout.write(timestr)
            sys.stdout.flush()
            time.sleep(1)
    
    def Tonight(self):
        sun=ephem.Sun()
        print 'Sunset: %s UT' % (self.Observer.next_setting(sun))
        print 'Sunrise: %s UT' % (self.Observer.next_rising(sun))
        # TODO: Moon rise, moon illumination
    
    def HourAngle(self,pos):
        # Todo!
        pass
        
class Scope:
    """Telescope Class
    name: Telescope Name (str)
    diameter: telescope diameter in meters
    f_ratio: telescope f-ratio
    mount: [Alt/Az,equitorial,transit,fixed] 
    style: 
    filts: dictionary of filt instances
    """
    def __init__(self,name='Empty Telescope',diameter=0.0,f_ratio=0.0,
        mount='Unknown',style='Unknown',filts={}):
        self.name=name
        self.diameter=diameter
        self.f_ratio=f_ratio
        self.mount=mount
        self.style=style
        self.filts=filts
    
class filt:
    """Filter
    filt.wave = wavelength in cm
    filt.freq = frequeny in hertz
    filt.energy = energy in erg
    """
    def __init__(self,val,valtype='wave',fluxconv=0.0,zp=99,comment='None'):
        self.val = val
        self.valtype = valtype
        self.fluxconv = fluxconv
        self.zp = zp
        self.comment = comment
        acceptabletypes = ['wave','freq','energy']
        try:
            acceptabletypes.index(valtype)
        except:
            print "Cannot initialize, needs to be of type: ", acceptabletypes 
        if self.valtype == 'wave':
            self.wave = self.val
            self.freq = c/self.wave
            self.energy = h*c/self.wave
        if self.valtype == 'freq':
            self.freq = self.freq
            self.wave = c/self.freq
            self.energy = h*self.freq
        if self.valtype == 'energy':
            self.energy = self.energy
            self.freq = self.energy/h
            self.wave = h*c/self.energy

# UVOT - source: Poole et al. 2008
uvotfilts={
    'vv':filt(5402e-8,valtype='wave',fluxconv=2.614e-16,zp=17.89, comment='UVOT V Filter'),
    'bb':filt(4329e-8,valtype='wave',fluxconv=1.472e-16,zp=19.11, comment='UVOT B Filter'),
    'uu':filt(3501e-8,valtype='wave',fluxconv=1.63e-16, zp=18.34, comment='UVOT U Filter'),
    'w1':filt(2634e-8,valtype='wave',fluxconv=4.00e-16, zp=17.49, comment='UVOT UVW1 Filter'),
    'm2':filt(2231e-8,valtype='wave',fluxconv=8.50e-16, zp=16.82, comment='UVOT UVM2 Filter'),
    'w2':filt(2030e-8,valtype='wave',fluxconv=6.2e-16,  zp=17.35, comment='UVOT UVW2 Filter'),
    'wh':filt(3471e-8,valtype='wave',fluxconv=3.7e-17,  zp=20.29, comment='UVOT UVW2 Filter')} 

# Cousin RI 
cousinfilts={
    'Rc':filt(6470e-8,valtype='wave',comment='Cousins R'),
    'Ic':filt(7865e-8,valtype='wave',comment='Cousins I')}

# 2mass
twomassfilts={
    'J':filt(1.25e-4,valtype='wave',comment='PAIRITEL J Filter'),
    'H':filt(1.65e-4,valtype='wave',comment='PAIRITEL H Filter'),
    'Ks':filt(2.15e-4,valtype='wave',comment='PAIRITEL Ks Filter')}

# Sloan
sloanfilts = {
    'u':filt(3540e-8,valtype='wave',comment="Sloan u'"),
    'g':filt(4750e-8,valtype='wave',comment="Sloan g'"),
    'r':filt(6220e-8,valtype='wave',comment="Sloan r'"),
    'i':filt(7630e-8,valtype='wave',comment="Sloan i'"),
    'z':filt(9050e-8,valtype='wave',comment="Sloan z'")}

# Johnson
johnsonfilts = {
    'U':filt(3640e-8,valtype='wave',comment="Johnson U"),
    'B':filt(4420e-8,valtype='wave',comment="Johnson B"),
    'V':filt(5400e-8,valtype='wave',comment="Johnson V")}

# UKIDSS (WFCAM)
ukidssfilts = {
    'Z':filt(8820e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey Z'),
    'Y':filt(10310e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey Y'),
    'J':filt(12480e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey J'),
    'H':filt(16310e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey H'),
    'K':filt(22010e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey K')}

# Initialize telescope instances
# Observatory" instances, which contain "Scopes" that use "Filts"
# Maybe also add "Instruments" class
whipple_scopes = {
    'PAIRITEL':Scope(
        name='PAIRITEL',
        diameter=1.3,
        f_ratio=13.4,
        mount='Equatorial',
        style='Cassegrain',
        filts=twomassfilts),
    'MMTO':Scope(
        name='MMTO',
        diameter=6.5)}
lick_scopes = {
    'KAIT':Scope(
        name='KAIT',
        diameter=0.76,
        mount='Equatorial'),
    'Shane':Scope(
        name='Shane',
        diameter=3.0),
    'Nickel':Scope(
        name='Nickel',
        diameter=1.0)}

whipple = Obs(name='Fred Lawrence Whipple Observatory',
    scopes=whipple_scopes,
    latitude='31:40:52',
    longitude='-110:52:42',
    elevation=2606.0)

lick = Obs(name="Lick Observatory",
    scopes=lick_scopes,
    latitude='37:20:29',
    longitude='-121:38:34',
    elevation=1283.0)

palomar = Obs(name="Palomar Observatory",
    latitude='33:21:21',
    longitude='-116:51:50',
    elevation=1713.0)

keck = Obs(name="W. M. Keck Observatory",
    latitude='19:49:35',
    longitude='-155:28:27',
    elevation=4145.0)

gemini_n = Obs(name="Gemini North Observatory",
    latitude='19:49:26',
    longitude='-155:28:09',
    elevation=2722.0)

gemini_s = Obs(name="Gemini South Observatory",
    latitude='-30:14:27',
    longitude='-70:44:12',
    elevation=2722.0)

gtc = Obs(name="Gran Telescopio Canarias",
    latitude='28:45:24',
    longitude='-17:53:31.3',
    elevation=2267.0)

uvot = Scope(name='UVOT',
    filts=uvotfilts)


def uvotcr2flux(countrate,filt,red_corr=0.0):
    '''
    Convert Swift UVOT countrates into fluxes.  Uses the 
    
    Arguments: 
        countrate: Corrected UVOT Count Rate (c/s)
        filt: q.filt instance
        red_corr: A_wave in magnitudes 
        
    '''
    uvotflux = countrate * filt.fluxconv
    print filt.comment
    uvotmag = filt.zp - 2.5*math.log10(countrate)
    print '  mag:  ' , uvotmag
    print '  flux: ' , uvotflux, 'erg/s/cm^2/A'
    print '  flux: ' , 1e26 * (filt.wave**2 / c) * uvotflux * 1e11, 'uJy'
    if red_corr != 0.0:
        print " Reddening Corrected Values:"
        corr_uvotmag = uvotmag - red_corr
        corr_rate = 10**(0.4*(filt.zp - corr_uvotmag))
        corr_uvotflux = corr_rate * filt.fluxconv
        print '  mag:  ' , corr_uvotmag
        print '  flux: ' , corr_uvotflux, 'erg/s/cm^2/A'
        print '  flux: ' , 1e26 * (filt.wave**2 / c) * corr_uvotflux * 1e11, 'uJy'
        
def flux2mag():
    """docstring for flux2mag"""
    pass

def mag2flux(mag_1=99,mag_2=99,flux_2=0):
    """Given a magnitude and a flux """
    if mag_1 == 99:
        mag_1 = raw_input('Please enter mag_1: ')
    if mag_2 == 99:
        mag_1 = raw_input('Please enter mag_2: ')
    if flux_2 == 0:
        flux_2 = raw_input('Please enter flux_2: ')
    flux_1 = flux_2 * 10**(-0.4*(mag_1 - mag_2))
    return flux_1

# Zeropoint: 1 count/second 
# Frequency: 5647 Angstroms hc/lambda

def mag2alpha(mag_1=None,mag_2=None,t_1=None,t_2=None):
    """Given a magnitude and a flux """
    if not mag_1:
        mag_1 = raw_input('Please enter mag_1: ')
    if not mag_2:
        mag_2 = raw_input('Please enter mag_2: ')
    if not t_1:
        t_1 = raw_input('Please enter t_1: ')
    if not t_2:
        t_2 = raw_input('Please enter t_2: ')
    try: 
        mag_1=float(mag_1) 
        mag_2=float(mag_2) 
        t_1=float(t_1)
        t_2=float(t_2) 
    except:
        return None           
    alpha = 0.4 * (mag_1 - mag_2)/(math.log10(t_1/t_2))
    return alpha

def tstart_exp_2_tmid_arr(burst_time=None,tstartarr=None,exparr=None,fmt='%y-%m-%d %H:%M:%S'):
    '''Takes the t_start and exposure time and calculates the t_mid
    
    Input: 
    burst_time = '09-11-27 23:25:45' 
    array of tstarts = ['09-11-28 01:02:20','09-11-28 01:06:20']
    array of exposure times = [123, 241]
    
    BE WARY THIS ASSUMES CONSTANT EXPOSURE i.e. tstart + exp = tstop
    
    Retuns array in seconds since burst
    
    '''
    
    import time 
    
    if not burst_time:
        burst_time = raw_input('Enter GRB Time: ')
    if not tstartarr:
        tstartarr = raw_input('Enter array of start times: ')
    if not exparr:
        exparr = raw_input('Enter array of exposure times (seconds): ')
    if not fmt:
        y_n = raw_input('Time Format of "%y-%m-%d %H:%M:%S" OK?  y/n: ')
        if y_n.lower() == 'y' or y_n.lower() == 'yes':
            fmt = "%y-%m-%d %H:%M:%S"
        elif y_n.lower() == 'n' or y_n.lower() == 'no':
            fmt = raw_input('Enter Desired Format: ')
        else:
            print 'Need y or n answer. Try again'
            return
    
    burst_time_sse = time.mktime(time.strptime(burst_time,fmt))
    if len(tstartarr) != len(exparr):
        raise ValueError('Length of Arrays does not match')
    
    zippedtime = zip(tstartarr,exparr)
    tmid_since_burst_arr = []
    
    for tupl in zippedtime:
        tmid = time.mktime(time.strptime(tupl[0],fmt)) - burst_time_sse + tupl[1]
        tmid_since_burst_arr.append(tmid)
        print tmid
    
    return tmid_since_burst_arr

def dhms2h(d1=None,h1=None,m1=None,s1=None,d2=None,h2=None,m2=None,s2=None):
    if d1 == None:
        d2 = float(raw_input("GRB Day: "))
        h2 = float(raw_input("GRB Hour: "))
        m2 = float(raw_input("GRB Minute: "))
        s2 = float(raw_input("GRB Second: "))
        d1 = float(raw_input("PTEL Day: "))
        h1 = float(raw_input("PTEL hour: "))
        m1 = float(raw_input("PTEL minute: "))
        s1 = float(raw_input("PTEL second: "))
    hours1 =  d1*24.0 + float(h1) + m1/60.0 + s1/3600.0
    hours2 =  d2*24.0 + float(h2) + m2/60.0 + s2/3600.0
    print hours1-hours2
    
def dec2sex(dec_pos):
    '''
    Convert decimal degrees position float tuple (ra,dec) into 
    sexagesimal string tuple
    
    '''
    import math
    
    if type(dec_pos).__name__ != 'tuple':
        print 'Was Epecting Tuple position'
        return
    elif len(dec_pos) != 2:
        print 'Was Expecting tuple of length two'
        return
    # If integer values given, convert to float
    if type(dec_pos[0]).__name__ == 'int': dec_pos = (float(dec_pos[0]), dec_pos[1])
    if type(dec_pos[1]).__name__ == 'int': dec_pos = (dec_pos[0], float(dec_pos[1]))
    if type(dec_pos[0]).__name__ != 'float' and type(dec_pos[1]).__name__ != 'float':
        print 'Decimal Degrees entries need to be of type float'
        return
    rangeflag = 0
    if not 0.0 <= dec_pos[0] < 360.0:
        print 'RA out of range.  0.0 <= RA < 360.0'
        rangeflag = 1
    if not -90.0 <= dec_pos[1] <= 90.0:
        print 'Dec out of range.  -90.0 <= Dec <= 90.0'
        rangeflag = 1
    if rangeflag == 1:
        print 'Returning due to out of range error'
        return
    
    abs_dec = abs(dec_pos[1])
    dec = [0.,0.,0.]
    str_dec = ['','','']
    
    ra_hours = dec_pos[0]/15.0
    ra = [0.,0.,0.]
    str_ra = ['','','']
    
    dec[0] = math.floor(abs_dec) ; remain = abs_dec - dec[0]
    dec[1] = math.floor(remain * 60.0) ; remain = abs_dec - dec[0] - dec[1]/60.0
    dec[2] = remain * 3600.0
    if dec_pos[1] < 0.0: dec[0] *= -1
    
    ra[0] = math.floor(ra_hours) ; remain = ra_hours - ra[0]
    ra[1] = math.floor(remain * 60.0) ; remain = ra_hours - ra[0] - ra[1]/60.0
    ra[2] = remain * 3600.0
    
    index = 0
    # Added the split to deal with decimal seconds e.g. 17:51:8.568 -> 17:51:08.56
    while index < 3:
        str_dec[index] = str(dec[index]).rstrip('0').rstrip('.')
        if len(str_dec[index].split('.')[0]) == 1: 
            str_dec[index] = '0' + str_dec[index]
        index += 1
    index = 0
    while index < 3:
        str_ra[index] = str(ra[index]).rstrip('0').rstrip('.')
        if len(str_ra[index].split('.')[0]) == 1:
            str_ra[index] = '0' + str_ra[index]
        index += 1
    
    if dec_pos[1] >= 0.0: str_dec[0] = '+' + str_dec[0]
    
    # If the arcseconds/seconds is tiny enough, just round to zero
    if dec[2] < 1e-10: str_dec[2] = '00'
    if ra[2] < 1e-10: str_ra[2] = '00'
    
    # If a rounding error caused the seconds/arsceconds to be 60, 
    # take care of it.  Comment out and test dex2sex((123.0,24))
    # or (30.5, -80.5) for an example
    if str_ra[2] == '60':
        str_ra[1] = str(int(str_ra[1])+1) # add to minutes
        if len(str_ra[1]) == 1: # add extra zero if too short
            str_ra[1] = '0' + str_ra[1] 
        str_ra[2] = '00' # subtract from seconds
    if str_dec[2] == '60':
        str_dec[1] = str(int(str_dec[1])+1) # add to arcmin
        if len(str_dec[1]) == 1: # add extra zero if too short
            str_dec[1] = '0' + str_dec[1]
        str_dec[2] = '00' # subtract from arcseconds
        
    str_dec_str = "%s:%s:%s" % (str_dec[0],str_dec[1],str_dec[2])
    str_ra_str = "%s:%s:%s" % (str_ra[0],str_ra[1],str_ra[2])
    
    sex_pos = (str_ra_str,str_dec_str)
    return sex_pos

def sex2dec(sex_pos):
    '''
    Convert sexagesimal position string tuple into decimal degree tuple
    
    '''
    if type(sex_pos).__name__ != 'tuple':
        print 'Was Epecting Tuple position'
        return
    elif len(sex_pos) != 2:
        print 'Was Expecting tuple of length two'
        return
    elif type(sex_pos[0]).__name__ != 'str' and type(sex_pos[1]).__name__ != 'str':
        print 'Sexagesimal entries need to be of type string'
        return
    
    if sex_pos[0].find(':') != -1:
        ra_list = sex_pos[0].split(':')
        dec_list = sex_pos[1].split(':')
    elif sex_pos[0].find(' ') != -1:
        ra_list = sex_pos[0].split(' ')
        dec_list = sex_pos[1].split(' ')
    else:
        print "Positions not formatted correctly! ('12:34:56.7','-65:43:21.0')"
        print "Cannot split based on ':' or ' '"
        print sex_pos
        return
        
    try:
        ras = float(ra_list[2]); ram = int(ra_list[1]); rah = int(ra_list[0])
        decs = float(dec_list[2]); decm = int(dec_list[1]); decd = int(dec_list[0])
    except IndexError:
        print "Positions not formatted correctly! ('12:34:56.7','-65:43:21.0')"
        print "Cannot split into 3 positions"
        print sex_pos
        return
    except ValueError:
        print "Positions not formatted correctly! ('12:34:56.7','-65:43:21.0')"
        print "Cannot convert 3 split positions into ints and floats"
        print sex_pos
        return
    
    rangeflag = 0
    if not 0 <= rah <= 24:
        print 'RA hours out of range: 0 <= rah <= 24'
        rangeflag = 1
    if not 0 <= ram < 60: 
        print 'RA minutes out of range: 0 <= ram < 60'
        rangeflag = 1    
    if not 0.0 <= ras < 60.0:
        print 'RA seconds out of range: 0.0 <= ras < 60.0'
        rangeflag = 1
    if not -90 <= decd <= 90:
        print 'Dec degrees out of range: -90 <= decd <= 90'
        rangeflag = 1
    if not 0 <= decm < 60:
        print 'Dec arcminutes out of range: 0 <= decm < 60'
        rangeflag = 1
    if not 0.0 <= decs < 60.0:
        print 'Dec arcseconds out of range: 0.0 <= decs < 60.0'
        rangeflag = 1
    if rangeflag == 1:
        print 'Returning due to out of range error'
        return
    
    # IS THIS MATH RIGHT???  9/25/09 - fixed math
    
    ra_ddeg = (float(rah) + float(ram)/60. + ras/3600.)*15.0
    if dec_list[0][0] == '-':  #if it is a negative number
        dec_ddeg = float(decd) - float(decm)/60. - decs/3600.
    else:
        dec_ddeg = float(decd) + float(decm)/60. + decs/3600.
    
    ddeg_pos = (ra_ddeg,dec_ddeg)
    return ddeg_pos


def sphere_dist(ra1,dec1,ra2,dec2):
    cosdistance = math.cos(ra1-ra2)*math.cos(dec1)*math.cos(dec2) + \
        math.sin(dec1)*math.sin(dec2)
    distance_deg = math.acos(cosdistance)
    distance_asec = distance_deg*3600.0
    return distance_asec


def Standardize(mylist):
    '''Given a list or array, returns the "standardized" version of the list
    with a mean of zero and a standard deviation of one.
    UNTESTED
    
    '''
    arr = scipy.array(mylist)
    return (arr-arr.mean())/(arr.std())

def RemoveNaN(myarr):
    '''Returns a shortened numpy array with all the NaN values removed'''
    return(myarr[numpy.isfinite(myarr).nonzero()])

def where(a,val,cond='==',wherenot=False):
    """
    Analogous to the 'where' function in IDL
    See thread:
    http://mail.python.org/pipermail/python-list/2007-March/431322.html
    
    Returns a list of indices in the list 'a' where a[i] == val
    If wherenot=True, returns an index where a[i]!= val
    e.g. 
    >>> a=[0,1,2,5,42,2]
    >>> where(a,42)
    [4]
    >>> where(a,2)
    [2,5]
    >>> where(a,2,wherenot=True)
    [0,1,3,4]
    >>> where(a,999)
    []
    
    This is similar to the numpy array np.nonzero() function
    In [103]: a
    Out[103]: [1, 2, 3]

    In [104]: q.where(a,2,cond='<=')
    Out[104]: [0, 1]

    In [105]: npa=np.array(a)

    In [106]: np.nonzero(npa<=2)
    Out[106]: (array([0, 1]),)
    
    """
    if wherenot == True:
        cond='!='
    if cond == '==':
    	return [i for i in xrange(len(a)) if a[i]==val]
    elif cond == '!=':
    	return [i for i in xrange(len(a)) if a[i]!=val]
    elif cond == '>':
        return [i for i in xrange(len(a)) if a[i]>val]
    elif cond == '<':
        return [i for i in xrange(len(a)) if a[i]<val]
    elif cond == '<=':
        return [i for i in xrange(len(a)) if a[i]<=val]
    elif cond == '>=':
        return [i for i in xrange(len(a)) if a[i]>=val]
    else:
        raise ValueError('Unsupported Condition.  Choose ==,!=,>,<,>=,<=.')    



def object2dict(obj,include=[],force=True):
    '''Given an object with attributes, return a dictionary with the attribute
    names as keys and attribute values as values.  Default to include all 
    attributes; specify a list of strings of attributes to include.
    
    If force=True, an exception will not be raised if an item in include is 
    not an attribute of obj.
    
    Usage:
    objdict = object2dict(obj,include=['attr1','attr2'])
    '''
    objdict = {}
    if not include:
        include = dir(obj)
    for attr in include:
        if force:
            try: objdict.update({attr:getattr(obj,attr)})
            except: pass
        else: objdict.update({attr:getattr(obj,attr)})
            
    return objdict
    

def SaveSortedDictTxt(indict,outpath='prettydict.txt'):
    '''Given a dictionary, save a prettified version of it to a text file
    for human consumption.
    '''
    import pprint
    f = open(outpath,'w')
    keys = indict.keys()
    keys.sort()
    sorted_vals = map(indict.get,keys)
    for item in keys:
        myindex=keys.index(item)
        mystr = '\n\n** %s **\n' % (item)
        f.write(mystr)
        prettydict = pprint.pformat(sorted_vals[myindex])
        f.write(prettydict)
    f.close

