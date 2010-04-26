"""
big_mosaic_maker.py mosaics your PAIRITEL data after full pipeline3 reduction 
with manual control over input mosaic files. You have to run big_mosaic_maker.py
in the same directory that contains the reduction_output directory from 
pipeline3.

You will also need the following files:
    big_mosaic_maker.py         -   This file, obviously
    pairitel_redux.swarp        -   SWarp configuration file for mosacing

Dependencies:
    SWarp
    WCSTools (sethead)
    Python 2.5 or better
    pyfits (python module)

Example Usage: 
    python mosaic_maker.py -o PULSE.33.
        [Then user edits the newly created j_long_triplestacks.txt, 
        h_long_triplestacks.txt, and k_long_triplestacks.txt files to exclude 
        unwanted triplestacks from the final mosaics]
    python mosaic_maker.py -o PULSE.33.1 
        [This will use the recently edited text files to make mosaics]
        
Author: Christopher Klein
        Adam Morgan
Contact: cklein@astro.berkeley.edu
Date: began 2009/09/04
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
import pyfits
import glob
import os
from os import system
import sys
try:
    from multiprocessing import Pool
    from multiprocessing import cpu_count as cpuCount
    doparallel = 1
except:
    print "Parallel processing library not installed. Not paralellizing."
    doparallel = 0
from optparse import OptionParser
#-------------------------------------------------------------------------------
# START TIMING
from time import time
# Begin program execution timing.
start_time = time()
#-------------------------------------------------------------------------------
# USER-SPECIFIC PATHS AND VALUES
# We'll be using SWarp. Define the proper bin executable here. This is 
# particularly important if you have different installations and want to use a 
# specific version in the reduction. Also, note that the SWarp packaged with 
# Scisoft probably won't work properly on Intel 64-bit Macs. It is strongly 
# advised that you recompile and install SWarp from source if SWarp's coadd 
# mosaics are null images.
swarp_bin = "/Applications/scisoft//i386/bin/swarp"
# We'll also be using some WCSTools. There should be no problem with the Scisoft
# version of thses.
sethead_bin = "/Applications/scisoft//i386/bin/sethead"
# For parallel processing we need the number of available CPUs. You could hard-
# code this value to something else, but it would result in non-optimal 
# performance.
if doparallel == 1:
    numprocessors = cpuCount()
else:
    numprocessors = 1
# Debugging mode outputs intermediate image files to help track down problems.
# It requires more time to run because of the additional hard disk reads/writes.
DEBUG = False
#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# Simple function allowing parallelization of mosaicing with SWarp.
def run_swarp(command):
    system(command)
    return
#-------------------------------------------------------------------------------
# BEGIN MAIN PROGRAM

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'


parser = OptionParser()

parser.add_option("-o", "--obs", "--obs-id", "--obs-string", action="store",     
    type="string", dest="obs_string", help=("observation ID string of epoch " + 
    "to be reduced"))
parser.add_option("-p", "--prep",
                  action="store_true", dest="do_prep", default=False,
                  help=("generate SWarp input lists to edit manually before " + 
                  "running in mosaicing mode"))
# AMorgan adds option
parser.add_option("-w", "--wcs", "--do-wcs",
                  action="store_true", dest="do_wcs", default=False,
                  help=("Attempt to do wcs fitting after coaddition"))
(options, args) = parser.parse_args()
obs_string = options.obs_string
if not obs_string:
    print "No obs_string specified, exiting."
    sys.exit()
do_prep = options.do_prep
do_wcs = options.do_wcs

# OBTAIN WORKING DIRECTORIES 
reduction_output_directory = str(obs_string) + "-reduction_output"
triplestacks_path = (reduction_output_directory + "/" + str(obs_string) + 
    "_triplestacks")
triplestackweights_path = (reduction_output_directory + "/" + str(obs_string) + 
    "_triplestackweights")

## Obtain Obs String
globlist = glob.glob(obs_string+'*-reduction_output')
j_long_list = file("j_long_mosaics.txt", "w")
h_long_list = file("h_long_mosaics.txt", "w")
k_long_list = file("k_long_mosaics.txt", "w")

j_long_list_weights = file("j_long_mosaic_weights.txt", "w")
h_long_list_weights = file("h_long_mosaic_weights.txt", "w")
k_long_list_weights = file("k_long_mosaic_weights.txt", "w")

j_stop_list = []
j_start_list = []
h_stop_list = []
h_start_list = []
k_stop_list = []
k_start_list = []

for item in globlist:
    obsid = item.split('-reduction_output')[-2]
    realpath = os.path.realpath(item)
    j_path = realpath + '/' + obsid + '_mosaics/j_long_' + obsid + '_coadd.fits'
    j_w_path = realpath + '/' + obsid + '_mosaics/j_long_' + obsid + '_coadd.weight.fits'
    h_path = realpath + '/' + obsid + '_mosaics/h_long_' + obsid + '_coadd.fits'
    h_w_path = realpath + '/' + obsid + '_mosaics/h_long_' + obsid + '_coadd.weight.fits'
    k_path = realpath + '/' + obsid + '_mosaics/k_long_' + obsid + '_coadd.fits'
    k_w_path = realpath + '/' + obsid + '_mosaics/k_long_' + obsid + '_coadd.weight.fits'
    j_long_list.write(j_path+'\n')
    h_long_list.write(h_path+'\n')
    k_long_list.write(k_path+'\n')
    j_long_list_weights.write(j_w_path+'\n')
    h_long_list_weights.write(h_w_path+'\n')
    k_long_list_weights.write(k_w_path+'\n')
    
    # Obtain the start and stop times of the image
    j_hdulist = pyfits.open(j_path)
    j_header = j_hdulist[0].header
    j_stop_list.append(str(j_header["STOP_CPU"]))
    j_start_list.append(str(j_header["STRT_CPU"]))
    j_hdulist.close()
    
    h_hdulist = pyfits.open(h_path)
    h_header = h_hdulist[0].header
    h_stop_list.append(str(h_header["STOP_CPU"]))
    h_start_list.append(str(h_header["STRT_CPU"]))
    h_hdulist.close()
    
    k_hdulist = pyfits.open(k_path)
    k_header = k_hdulist[0].header
    k_stop_list.append(str(k_header["STOP_CPU"]))
    k_start_list.append(str(k_header["STRT_CPU"]))
    k_hdulist.close()

## Sort the lists of start and stop times
j_start_list.sort()
j_stop_list.sort()
h_start_list.sort()
h_stop_list.sort()
k_start_list.sort()
k_stop_list.sort()
## Take the first start time and the last stop time
j_earliest_start = j_start_list[0]
j_latest_stop = j_stop_list[-1]
h_earliest_start = j_start_list[0]
h_latest_stop = j_stop_list[-1]
k_earliest_start = j_start_list[0]
k_latest_stop = j_stop_list[-1]


# Close the relevant files
j_long_list.close()
j_long_list_weights.close()
h_long_list.close()
h_long_list_weights.close()
k_long_list.close()
k_long_list_weights.close()

# Run the mosaicing. 
# Make list of swarp_commands.
swarp_commands = [
    swarp_bin + " @j_long_mosaics.txt " + 
    "-c " + loadpath + "pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @j_long_mosaic_weights.txt " + 
    "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
    "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits",
    swarp_bin + " @h_long_mosaics.txt " + 
    "-c " + loadpath + "pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @h_long_mosaic_weights.txt " + 
    "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
    "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits",
    swarp_bin + " @k_long_mosaics.txt " + 
    "-c " + loadpath + "pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE @k_long_mosaic_weights.txt " + 
    "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
    "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits"]

if doparallel == 1:
    # Run the mosaicing with parallel processing.
    p = Pool(numprocessors)
    result = p.map_async(run_swarp, swarp_commands)
    poolresult = result.get()
else:
    # Run the mosaicing without parallel processing.
    for command in swarp_commands:
        run_swarp(command)

# We insert the STRT_CPU of the first triplestack and the STOP_CPU of 
# the last triplestack used to make the mosaic.

j_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits",mode='update')
h_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits",mode='update')
k_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits",mode='update')
j_header = j_hdulist[0].header
h_header = h_hdulist[0].header
k_header = k_hdulist[0].header
j_header.update('STRT_CPU',j_earliest_start)
j_header.update('STOP_CPU',j_latest_stop)
h_header.update('STRT_CPU',h_earliest_start)
h_header.update('STOP_CPU',h_latest_stop)
k_header.update('STRT_CPU',k_earliest_start)
k_header.update('STOP_CPU',k_latest_stop)

j_hdulist.flush()
h_hdulist.flush()
k_hdulist.flush()

j_hdulist.close()
h_hdulist.close()
k_hdulist.close()

# Remove extraneous text files
system('rm ?_long*mosaic*.txt')

# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time