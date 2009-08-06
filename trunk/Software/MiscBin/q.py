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