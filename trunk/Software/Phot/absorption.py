#!/usr/bin/env python

'''
Absorption Calculator

Python port of Doug Welch's Excellent Absorption Law calculator
http://wwwmacho.mcmaster.ca/JAVA/Acurve.html

This calculator determines the total absorption at wavelengths 
between 0.10 and 3.33 um using the algorithm determined by 
Cardelli, Clayton, and Mathis (1989) [ApJ 345 245].

inputs: 
rv the ratio of total to selective absorption at V
av the total absorption in mags at V
wv the target wavelength (in um)!

'''
__author__ = 'Doug Welch, ported to Python by Adam N Morgan'
__email__ = 'qmorgan@gmail.com'

import sys

def usage():

    print """
usage: absorption.py <A_V> <LAMBDA> <R_V[optional]>

  where <A_V> is the total absorption in mags at V 
  <LAMBDA> is the wavelength in MICROMETERS (um)
  <R_V> is the ratio of total to selective absorption at V
     (this is an optional parameter; if no value is specified, 
     the commonly accepted value of 3.1 for the Milky Way is 
     assumed)
  
  """
    sys.exit()

def AbsCalc(wv,av,rv=3.1):
    x = 1.0/wv
    
    if x >= 0.3 and x < 1.1:
        ax = 0.574*x**1.61
        bx = -0.527*x**1.61
    
    elif x >=1.1 and x < 3.3:  # 0.9091 > wv > 0.3030
        y = x - 1.82
        ax = 1.0 + (0.17699 - 0.50447*y)*y
        ax -= 0.02427*y**3
        ax += 0.72085*y**4
        ax += 0.01979*y**5
        ax -= 0.77530*y**6
        ax += 0.32999*y**7
        
        bx = 1.41338*y
        bx += 2.28305*y**2
        bx += 1.07233*y**3
        bx -= 5.38434*y**4
        bx -= 0.62251*y**5
        bx += 5.30260*y**6
        bx -= 2.09002*y**7
        
    elif x >= 3.3 and x < 8.0:  # 0.3030 > wv > 0.125
        if x >= 5.9 and x <8.0:
            xx = x-5.9
            fa = (xx**2)*(-0.04473 - 0.009779*xx)
            fb = (xx**2)*( 0.2130  + 0.1207*xx)
        else:
            fa = 0.0
            fb = 0.0
        
        ax = 1.752 - 0.316*x 
        ax -= 0.104/((x-4.67)**2 + 0.341)
        ax += fa
        
        bx = -3.090 + 1.825*x
        bx += 1.206/((x-4.62)**2 + 0.263)
        bx += fb
    
    elif x >= 8.0 and x <= 10.0:  # 0.125 > wv > 0.1
        
        xx = (x - 8.0)
        
        ax = -1.073 - 0.628*xx + 0.137*xx**2
        ax -= 0.070*xx**3
        
        bx = 13.670 + 4.257*xx - 0.420*xx**2
        bx += 0.374*xx**3
        
    
    else:
        raise ValueError('Wavelength out of Range')
        
    al = (ax + bx/rv) * av
    return al
    
if __name__ == "__main__":
    
    # invoked from the command line
    if len(sys.argv) == 1 or len(sys.argv) > 4:
        print usage()
        sys.exit(0)
        
    if ((sys.argv[1] == "-h") or (sys.argv[1] == "--h")):
        print usage()
        sys.exit(0)
        
    wv = float(sys.argv[1])
    av = float(sys.argv[2])
    if len(sys.argv) > 3:
        rv = float(sys.argv[3])
    else:
        rv = 3.1
    
    a_wav = AbsCalc(wv,av,rv)
    print a_wav
    sys.exit()