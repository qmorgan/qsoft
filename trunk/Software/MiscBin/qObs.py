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
