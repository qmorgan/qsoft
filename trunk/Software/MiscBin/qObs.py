#!/usr/bin/env python
# encoding: utf-8
"""
qObs.py
Author: Adam Morgan

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
    
    def GetMoon(self,silent=False):
        '''Calculate the relative position of the moon'''
        self.Observer.date=ephem.now()
        self.moon = ephem.Moon()
        self.moon.compute(self.Observer)
        if not silent:
            print '*******   Moon Position   *******'
            print 'RA, Dec: (%s,%s)' % (self.moon.ra,self.moon.dec)
            print 'Alt, Az: (%s,%s)' % (self.moon.alt,self.moon.az)
            print 'Illumination: %1.4f' % (self.moon.moon_phase)
            print 'Current time: %s UT' % (self.Observer.date)
            print 'Moonset: %s UT' % (self.Observer.next_setting(self.moon))
            print 'Moonrise: %s UT' % (self.Observer.next_rising(self.moon))
    
    def GetSun(self,silent=False):
        '''Calculate the relative position of the Sun'''
        self.Observer.date=ephem.now()
        self.sun = ephem.Sun()
        self.sun.compute(self.Observer)
        if not silent:
            print '*******    Sun Position    *******'
            print 'RA, Dec: (%s,%s)' % (self.sun.ra,self.sun.dec)
            print 'Alt, Az: (%s,%s)' % (self.sun.alt,self.sun.az)
            print 'Current time: %s UT' % (self.Observer.date)
            print 'Sunset: %s UT' % (self.Observer.next_setting(self.sun))
            print 'Sunrise: %s UT' % (self.Observer.next_rising(self.sun))
    
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
    def __init__(self,val,valtype='wave',fwhm=0.0,zpflux=0.0,kval=0.0,fluxconv=0.0,zp=99,name='none',comment='None'):
        '''Given a val and a valtype, convert the relevant value into wavelength,
        energy, and frequency units.  
        
        wave: wavelength in centimeters
        freq: frequency in hz
        energy: energy in erg
        
        fwhm: same units as valtype
        zpflux: flux corresponding to a magnitude of 0 (in microjanskies)
        kval: such that F(uJy) = 10^(kval - mag/2.5)
        fluxconv: flux density conversion factor for uvot in ergs cm^-2 s^-1 angstrom^-1
        zp: zero point magnitude (currently only for uvot)
        
        name: name of the filter for plot annotations
        comment: more info about the filter
        
        Given a zpflux, determine the kval for flux conversions, where
        F = 10^(kval - mag/2.5)
        log(F) = kval - mag/2.5
        2.5log(F) = 2.5kval - mag
        and if zeropoint, mag=0, so if zpflux is in microjanskies, 
        log(F) = log(zpflux) = kval 
        
        If zpflux and kval are both given, raise an exception. Only one should 
        be given, and the other can be converted.
        
        For convinience, energy conversions are done:
        wave_A: wavelength in angstroms
        wave_nm: wavelength in nanometers
        wave_um: wavelength in microns 
        wave_m: wavelength in meters
        
        energy_eV: energy in electron volts
        energy_J: energy in Joules
        
        
        '''
        self.zpflux = zpflux
        self.kval = kval
        
        if zpflux and kval:
            print 'Cannot initialize; only either zpflux or kval can be specified.'
            return
        if zpflux and not kval:
            self.kval = numpy.log10(zpflux) #assuming zpflux in uJy         
        if kval and not zpflux:            
            self.zpflux = 10**kval # gives zpflux in uJy
        
        self.zpflux_uJy = self.zpflux #no conversion
        self.zpflux_Jy = self.zpflux*1E-6
        self.zpflux_cgs = self.zpflux_Jy*1E-23 #erg/cm^2/s/Hz
        self.zpflux_mks = self.zpflux_Jy*1E-26 #W/m^2/Hz
            
        self.val = val
        self.valtype = valtype
        self.fluxconv = fluxconv
        self.zp = zp
        self.fwhm = fwhm
        self.comment = comment
        self.name = name
        acceptabletypes = ['wave','freq','energy']
        try:
            acceptabletypes.index(valtype)
        except:
            print "Cannot initialize, needs to be of type: ", acceptabletypes 
            return
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
        
        # do unit conversions for convinience
        self.wave_cm = self.wave #no conversion
        self.wave_A = self.wave * 1E8
        self.wave_nm = self.wave * 1E7
        self.wave_um = self.wave * 1E4
        self.wave_m = self.wave * 1E-2
        
        self.freq_Hz = self.freq #no conversion
        
        self.energy_erg = self.energy #no conversion
        self.energy_eV = self.energy * 6.24150974*1E11
        self.energy_J = self.energy * 1E-7

# UVOT - source: Poole et al. 2008 and heasarc.nasa.gov/docs/swift/analysis/ for fwhm
# white values from roming et al 2005 (verified in arxiv/0809.4193)
VV=filt(5468e-8,valtype='wave',fwhm=769e-8,fluxconv=2.614e-16,zp=17.89,name='UVOT V', comment='UVOT V Filter')
BB=filt(4392e-8,valtype='wave',fwhm=975e-8,fluxconv=1.472e-16,zp=19.11,name='UVOT B', comment='UVOT B Filter')
UU=filt(3465e-8,valtype='wave',fwhm=785e-8,fluxconv=1.63e-16, zp=18.34,name='UVOT U', comment='UVOT U Filter')
W1=filt(2600e-8,valtype='wave',fwhm=700e-8,fluxconv=4.00e-16, zp=17.49,name='UVW1', comment='UVOT UVW1 Filter')
M2=filt(2246e-8,valtype='wave',fwhm=510e-8,fluxconv=8.50e-16, zp=16.82,name='UVM2', comment='UVOT UVM2 Filter')
W2=filt(1928e-8,valtype='wave',fwhm=760e-8,fluxconv=6.2e-16,  zp=17.35,name='UVW2', comment='UVOT UVW2 Filter')
WH=filt(3850e-8,valtype='wave',fwhm=2600e-8,fluxconv=3.7e-17,  zp=20.29,name='White', comment='UVOT White Filter')
uvotfilts={'vv':VV,'bb':BB,'uu':UU,'w1':W1,'m2':M2,'w2':W2,'wh':WH} 

# Cousin RI 
# Fukujita et al. 1995
Rc=filt(6588e-8,valtype='wave',fwhm=1568e-8,zpflux=3104e6,name='R',comment='Cousins R')
Ic=filt(8060e-8,valtype='wave',fwhm=1542e-8,zpflux=2432e6,name='I',comment='Cousins I')
cousinfilts={'Rc':Rc,'Ic':Ic}

# 2mass
# cohen et al. 2003
# zp (jy) J:1594, H:1024, Ks:666.7
J=filt(1.235e-4,valtype='wave',fwhm=1620e-8,zpflux=1594e6,name='J',comment='PAIRITEL J Filter')
H=filt(1.662e-4,valtype='wave',fwhm=2510e-8,zpflux=1024e6,name='H',comment='PAIRITEL H Filter')
Ks=filt(2.159e-4,valtype='wave',fwhm=2620e-8,zpflux=666.7e6,name='Ks',comment='PAIRITEL Ks Filter')
twomassfilts={'J':J,'H':H,'Ks':Ks}

# Sloan
# Fukujita et al. 1995
# http://www.sdss3.org/dr8/algorithms/fluxcal.php#SDSStoAB
# The SDSS photometry is intended to be on the AB system (Oke & Gunn 1983), 
# by which a magnitude 0 object should have the same counts as a source of 
# Fν = 3631 Jy. However, this is known not to be exactly true, such that the 
# photometric zeropoints are slightly off the AB standard. We continue to work 
# to pin down these shifts. Our present estimate, based on comparison to the 
# STIS standards of Bohlin, Dickinson, & Calzetti (2001) and confirmed by SDSS 
# photometry and spectroscopy of fainter hot white dwarfs, is that the u band 
# zeropoint is in error by 0.04 mag, uAB = uSDSS - 0.04 mag, and that g, r, and 
# i are close to AB. These statements are certainly not precise to better than 
# 0.01 mag; in addition, they depend critically on the system response of the
# SDSS 2.5-meter, which was measured by Doi et al. (2010). The z band zeropoint 
# is not as certain at this time, but there is mild evidence that it may be s
# hifted by about 0.02 mag in the sense zAB = zSDSS + 0.02 mag. The large shift
#  in the u band was expected because the adopted magnitude of the SDSS standard 
#  BD+17 in Fukugita et al. (1996) was computed at zero airmass, thereby making 
# the assumed u response bluer than that of the USNO system response.
# kval = log(3631e6ujy) = 9.56
# u: 9.56 - 0.04/2.5 = 9.544 => zpflux = 3500
# z: 9.56 + 0.02/2.5 = 9.568 => zpflux = 3698
u=filt(3585e-8,valtype='wave',fwhm=556e-8,zpflux=3500e6,name='u',comment="Sloan u'")
g=filt(4858e-8,valtype='wave',fwhm=1297e-8,zpflux=3631e6,name='g',comment="Sloan g'")
r=filt(6290e-8,valtype='wave',fwhm=1358e-8,zpflux=3631e6,name='r',comment="Sloan r'")
i=filt(7706e-8,valtype='wave',fwhm=1547e-8,zpflux=3631e6,name='i',comment="Sloan i'")
z=filt(9222e-8,valtype='wave',fwhm=1530e-8,zpflux=3698e6,name='z',comment="Sloan z'")
sloanfilts = {'u':u,'g':g,'r':r,'i':i,'z':z}

# Johnson
# Fukujita et al. 1995
U=filt(3652e-8,valtype='wave',fwhm=526e-8,zpflux=1923e6,name='U',comment="Johnson U")
B=filt(4458e-8,valtype='wave',fwhm=1008e-8,zpflux=4130e6,name='B',comment="Johnson B")
V=filt(5505e-8,valtype='wave',fwhm=827e-8,zpflux=3689e6,name='V',comment="Johnson V")
johnsonfilts = {'U':U,'B':B,'V':V}

# UKIDSS (WFCAM)
# NOT VERIFIED VALUES; be warned
ukZ=filt(8820e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey Z')
ukY=filt(10310e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey Y')
ukJ=filt(12480e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey J')
ukH=filt(16310e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey H')
ukK=filt(22010e-8,valtype='wave',comment='UKIRT IR Deep Sky Survey K')
ukidssfilts = {'Z':ukZ,'Y':ukY,'J':ukJ,'H':ukH,'K':ukK}

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

uvot = Scope(name='UVOT',
    filts=uvotfilts)


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
    
