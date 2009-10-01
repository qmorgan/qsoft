#!/usr/bin/env python
# encoding: utf-8
"""
q.py
Author: Adam Morgan
Created: Aug 2, 2009
Last Updated: Aug 2, 2009
	Created
	
A conglomeration of random definitions that I use!.
"""
import sys
import os
import math

# Constants - all CGS
c = 2.99792458e10 # cm/s
h = 6.626e-27 # erg s

#import units

class filt:
    """Filter"""
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
    
vv = filt(5402e-8,valtype='wave',fluxconv=2.614e-16,zp=17.89, comment='UVOT V Filter')
bb = filt(4329e-8,valtype='wave',fluxconv=1.472e-16,zp=19.11, comment='UVOT B Filter')
uu = filt(3501e-8,valtype='wave',fluxconv=1.63e-16, zp=18.34, comment='UVOT U Filter')
w1 = filt(2634e-8,valtype='wave',fluxconv=4.00e-16, zp=17.49, comment='UVOT UVW1 Filter')
m2 = filt(2231e-8,valtype='wave',fluxconv=8.50e-16, zp=16.82, comment='UVOT UVM2 Filter')
w2 = filt(2030e-8,valtype='wave',fluxconv=6.2e-16,  zp=17.35, comment='UVOT UVW2 Filter')
wh = filt(3471e-8,valtype='wave',fluxconv=3.7e-17,  zp=20.29, comment='UVOT UVW2 Filter')

def uvotcr2flux(countrate,filt,red_corr=0.0):
    '''Use http://wwwmacho.mcmaster.ca/JAVA/Acurve.html to find absoprtion
    corrections.
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

