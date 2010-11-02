"""
scamp_wcs_fit.py fits pipeline3 mosaic images with WCS header info.

You will also need the following files:
    gauss_2.0_3x3.conv          -   SExtractor filter for alignment catalogs
    pairitel_align.param        -   SExtractor output catalog parameters file
                                    for alignment catalog production
    pairitel_align.sex          -   SExtractor configuration file for alignment
                                    catalog production
    pairitel_align.scamp        -   SCAMP configuration file for WCS fitting
    scamp_wcs_fit_dir.py        -   This file, obviously
    
Example program call:
    python scamp_wcs_fit_dir.py [image_directory]
    python scamp_wcs_fit_dir.py final_mosaics

Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/12/03
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
from os import system, listdir
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
image_directory = sys.argv[1]
raw_file_list = listdir(image_directory)
file_list = []
weight_file_list = []
for filename in raw_file_list:
    if (filename[-5:] == ".fits") and (filename[-10:] == "coadd.fits"):
        file_list.append(image_directory + "/" + filename)
        weight_file_list.append(image_directory + "/" + filename[:-5] + ".weight.fits")

print file_list
for n in range(len(file_list)):
    # Run SExtractor on each 
    system(sextractor_bin + " " + file_list[n] + " -c pairitel_align.sex " + 
        "-CATALOG_NAME " + file_list[n][:-5] + ".ldac " + 
        "-WEIGHT_IMAGE " + weight_file_list[n])
    # And then, run Scamp.
    system(scamp_bin + " " + file_list[n][:-5] + ".ldac " + 
        "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
    # Finally, use missfits to apply the new WCS header info to the image.
    system(missfits_bin + " " + file_list[n])
    # Finish up by copying the WCS header info into the weightmap.
    system("cp " + file_list[n][:-5] + ".head " + 
        weight_file_list[n][:-5] + ".head")
    system(missfits_bin + " " + weight_file_list[n])

# Clean up the directory of all the intermediate files.
# system("rm " + image_directory + "/missfits.xml")
system("rm " + image_directory + "/*.back")
system("rm " + image_directory + "/*.ldac")
system("rm " + image_directory + "/*.head")

# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time
