import matplotlib
import os
import numpy
from pylab import *
from numpy import array as arr
import cosmocalc

# autophot, '051008_g.fits', '051008.sdss', 'allhosts.pos', rad=1.2

def p_phot(burst_image,calfile,objectfile,ap=1.2):

    # band = band.upper() #making sure band in upper case
    import pidly
    idl_path = '/Applications/itt/idl71/bin/idl'
    idl = pidly.IDL(idl_path)
    #Performing Photometry
    IDL_command = "autophot, '" + str(burst_image) + "', '" + str(calfile) +  "', '"  + str(objectfile) + "', rad=" + str(ap)
    idl(IDL_command)
    # 
    # #Read the filename
    filename = 'tmpphot.txt'
    f=open(filename,'r')
    lines = f.readlines()
    for line in lines:
        print line
    #     if band in line:
    #         flux = (line.split()[1], line.split()[2])
    # # 
    # return flux
