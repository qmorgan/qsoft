"""
scamp_wcs_fit.py fits pipeline3 mosaic images with WCS header info.

Run the program in a directory with your mosiacs.

You will also need the following files:
    gauss_2.0_3x3.conv          -   SExtractor filter for alignment catalogs
    pairitel_align.param        -   SExtractor output catalog parameters file
                                    for alignment catalog production
    pairitel_align.sex          -   SExtractor configuration file for alignment
                                    catalog production
    pairitel_align.scamp        -   SCAMP configuration file for WCS fitting
    scamp_wcs_fit.py            -   This file, obviously
    
Example program call:
    python scamp_wcs_fit.py [obj_string]
    python scamp_wcs_fit.py MFS.1.2

Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/11/18
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
from os import system
import sys
#-------------------------------------------------------------------------------
# START TIMING
from time import time
# Begin program execution timing.
start_time = time()
#-------------------------------------------------------------------------------
# USER-SPECIFIC PATHS AND VALUES
# We'll be using SExtractor. Define the proper bin executable here. 
# This is particularly important if you have different installations and want to
# use a specific version in the reduction. Also, note that the SExtractor 
# packaged with Scisoft probably won't work properly on Intel Macs. It is
# strongly advised that you recompile and install SExtractor from 
# source if SExtractor's check images are null images.
sextractor_bin = "/opt/local/bin/sex"
# Same for SCAMP and MissFITS.
scamp_bin = "/opt/local/bin/scamp"
missfits_bin = "python /Users/amorgan/q_soft/trunk/Software/MiscBin/missfits.py"
#-------------------------------------------------------------------------------
# BEGIN MAIN PROGRAM
obj_string = sys.argv[1]
# Next, run SExtractor on each 
system(sextractor_bin + " j_long_" + obj_string + "_coadd.fits " + 
    "-c pairitel_align.sex " + 
    "-CATALOG_NAME j_long_" + obj_string + "_coadd.ldac " + 
    "-WEIGHT_IMAGE j_long_" + obj_string + "_coadd.weight.fits")
system(sextractor_bin + " h_long_" + obj_string + "_coadd.fits " + 
    "-c pairitel_align.sex " + 
    "-CATALOG_NAME h_long_" + obj_string + "_coadd.ldac " + 
    "-WEIGHT_IMAGE j_long_" + obj_string + "_coadd.weight.fits")
system(sextractor_bin + " k_long_" + obj_string + "_coadd.fits " + 
    "-c pairitel_align.sex " + 
    "-CATALOG_NAME k_long_" + obj_string + "_coadd.ldac " + 
    "-WEIGHT_IMAGE j_long_" + obj_string + "_coadd.weight.fits")
# And then, run Scamp.
system(scamp_bin + " j_long_" + obj_string + "_coadd.ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
system(scamp_bin + " h_long_" + obj_string + "_coadd.ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
system(scamp_bin + " k_long_" + obj_string + "_coadd.ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
# Finally, use missfits to apply the new WCS header info to the image.
system(missfits_bin + " j_long_" + obj_string + "_coadd.fits")
system(missfits_bin + " h_long_" + obj_string + "_coadd.fits") 
system(missfits_bin + " k_long_" + obj_string + "_coadd.fits") 
# Finish up by copying the WCS header info from the ?_long*coadd.fits mosaics 
# into the corresponding weightmaps and ?_short*coadd.fits mosaics and 
# weightmaps.
system("cp j_long_" + obj_string + "_coadd.head j_long_" + obj_string + 
    "_coadd.weight.head")
system("cp j_long_" + obj_string + "_coadd.head j_short_" + obj_string + 
    "_coadd.head")
system("cp j_long_" + obj_string + "_coadd.head j_short_" + obj_string + 
    "_coadd.weight.head")
system("cp h_long_" + obj_string + "_coadd.head h_long_" + obj_string + 
    "_coadd.weight.head")
system("cp h_long_" + obj_string + "_coadd.head h_short_" + obj_string + 
    "_coadd.head")
system("cp h_long_" + obj_string + "_coadd.head h_short_" + obj_string + 
    "_coadd.weight.head")
system("cp k_long_" + obj_string + "_coadd.head k_long_" + obj_string + 
    "_coadd.weight.head")
system("cp k_long_" + obj_string + "_coadd.head k_short_" + obj_string + 
    "_coadd.head")
system("cp k_long_" + obj_string + "_coadd.head k_short_" + obj_string + 
    "_coadd.weight.head")
system(missfits_bin + " j_long_" + obj_string + "_coadd.weight.fits")
system(missfits_bin + " h_long_" + obj_string + "_coadd.weight.fits") 
system(missfits_bin + " k_long_" + obj_string + "_coadd.weight.fits")
system(missfits_bin + " j_short_" + obj_string + "_coadd.fits")
system(missfits_bin + " h_short_" + obj_string + "_coadd.fits") 
system(missfits_bin + " k_short_" + obj_string + "_coadd.fits")
system(missfits_bin + " j_short_" + obj_string + "_coadd.weight.fits")
system(missfits_bin + " h_short_" + obj_string + "_coadd.weight.fits") 
system(missfits_bin + " k_short_" + obj_string + "_coadd.weight.fits")
# Clean up the directory of all the intermediate files.
system("rm missfits.xml")
#system("rm *.back")
system("rm *.ldac")
system("rm *.head")
# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time
