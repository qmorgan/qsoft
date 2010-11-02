"""
big_swarp.py swarps pipeline3 mosaic images into a deep mosaic.

You will also need the following files:
    pairitel_redux.swarp        -   SCAMP configuration file for WCS fitting
    big_swarp.py                -   This file, obviously
    
Example program call:
    python big_swarp.py [image_directory]
    python big_swarp.py final_mosaics_fit

Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/12/03
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
from os import system, listdir
import pyfits
import sys
#-------------------------------------------------------------------------------
# START TIMING
from time import time
# Begin program execution timing.
start_time = time()
#-------------------------------------------------------------------------------
# USER-SPECIFIC PATHS AND VALUES
swarp_bin = "/Applications/scisoft//i386/bin/swarp"
sethead_bin = "/Applications/scisoft//i386/bin/sethead"
#-------------------------------------------------------------------------------
# BEGIN MAIN PROGRAM
image_directory = sys.argv[1]
raw_file_list = listdir(image_directory)
j_files = file("j_input_images.txt", "w")
j_weights = file("j_input_weights.txt", "w")
h_files = file("h_input_images.txt", "w")
h_weights = file("h_input_weights.txt", "w")
k_files = file("k_input_images.txt", "w")
k_weights = file("k_input_weights.txt", "w")

for filename in raw_file_list:
    if (filename[-5:] == ".fits") and (filename[-10:] == "coadd.fits"):
        if filename[0] == "j":
            j_files.write(image_directory + "/" + filename + "\n")
            j_weights.write(image_directory + "/" + filename[:-5] + 
                ".weight.fits\n")
            
        if filename[0] == "h":
            h_files.write(image_directory + "/" + filename + "\n")
            h_weights.write(image_directory + "/" + filename[:-5] + 
                ".weight.fits\n")
        if filename[0] == "k":
            k_files.write(image_directory + "/" + filename + "\n")
            k_weights.write(image_directory + "/" + filename[:-5] + 
                ".weight.fits\n")
j_files.close()
j_weights.close()
h_files.close()
h_weights.close()
k_files.close()
k_weights.close()

system(swarp_bin + " @j_input_images.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @j_input_weights.txt " + 
    "-IMAGEOUT_NAME j_mosaic_coadd.fits " + 
    "-WEIGHTOUT_NAME j_mosaic_coadd.weight.fits")
system(swarp_bin + " @h_input_images.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @h_input_weights.txt " + 
    "-IMAGEOUT_NAME h_mosaic_coadd.fits " + 
    "-WEIGHTOUT_NAME h_mosaic_coadd.weight.fits")
system(swarp_bin + " @k_input_images.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @k_input_weights.txt " + 
    "-IMAGEOUT_NAME k_mosaic_coadd.fits " + 
    "-WEIGHTOUT_NAME k_mosaic_coadd.weight.fits")

system("rm j_input_images.txt")
system("rm j_input_weights.txt")
system("rm h_input_images.txt")
system("rm h_input_weights.txt")
system("rm k_input_images.txt")
system("rm k_input_weights.txt")

# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time