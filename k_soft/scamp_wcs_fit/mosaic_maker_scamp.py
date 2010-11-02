"""
mosaic_maker.py will make mosaics over the specified interval (between p_start
and p_end, inclusive). You need to specify the obj_string, p_start, and p_end
in the program call.
Usage: python mosaic_maker.py obj_string p_start p_stop

obj_string is the PAIRITEL database id for the observation (Ex: SN.175).
p_start and p_stop are integers. They define the range of p# triplestacks
to mosaic.

The progenitor triplestacks need to be in a subdirectory 
"./" + obj_string + "_triplestacks". Also, the weightmaps need to be in a 
subdirectory "./" + obj_string + "_weightmaps". 

You will also need the following files:
    gauss_2.0_3x3.conv          -   SExtractor filter for alignment catalogs
    pairitel_align.param        -   SExtractor output catalog parameters file
                                    for alignment catalog production
    pairitel_align.sex          -   SExtractor configuration file for alignment
                                    catalog production
    pairitel_align.scamp        -   SCAMP configuration file for WCS fitting
    mosaic_maker2.py            -   This file, obviously

If you are calling mosaic_maker.py after running the full reduction in the 
current working directory, then all the dependencies should already be 
satisfied.

Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/07/01
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
import pyfits
from os import listdir
from os import system
import sys
from multiprocessing import Pool
from multiprocessing import cpu_count
#-------------------------------------------------------------------------------
# START TIMING
from time import time
# Begin program execution timing.
start_time = time()
#-------------------------------------------------------------------------------
# USER-SPECIFIC PATHS AND VALUES
# We'll be using SExtractor and SWarp. Define the proper bin executable here. 
# This is particularly important if you have different installations and want to
# use a specific version in the reduction. Also, note that the SExtractor and
# SWarp packaged with Scisoft probably won't work properly on Intel Macs. It is
# strongly advised that you recompile and install SExtractor and SWarp from 
# source if SExtractor's check images or SWarp's coadd mosaics are null images.
sextractor_bin = "/opt/local/bin/sex"
swarp_bin = "/Applications/scisoft//i386/bin/swarp"
# Same for SCAMP and MissFITS.
scamp_bin = "/opt/local/bin/scamp"
missfits_bin = "python /Users/amorgan/q_soft/trunk/Software/MiscBin/missfits.py"
# We'll also be using some WCSTools. There should be no problem with the Scisoft
# version of thses.
imwcs_bin = "/Applications/scisoft/i386/bin/imwcs"
cphead_bin = "/Applications/scisoft/i386/bin/cphead"
sethead_bin = "/Applications/scisoft/i386/bin/sethead"
# For parallel processing we need the number of available CPUs. You could hard-
# code this value to something else, but it would result in non-optimal 
# performance.
numprocessors = cpu_count()
#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# Simple function allowing parallelization of mosaicing with SWarp.
def run_swarp(command):
    system(command)
    return
#-------------------------------------------------------------------------------
# BEGIN MAIN PROGRAM
if len(sys.argv) != 4:
    print "Usage: python mosaic_maker.py obj_string p_start p_stop"
    print "obj_string is the PAIRITEL database object id (Ex: SN.175)."
    print "p_start is the p# of the first triplestack to mosaic."
    print "p_end is the p# of the last triplestack to mosaic."
    sys.exit()
# Read in the user's arguments.
obj_string = sys.argv[1]
p_start = int(sys.argv[2])
p_end = int(sys.argv[3])
if (p_start <= 0) or (p_start >= p_end):
    print "p_start must be >= 1 and p_start < p_end."
    sys.exit()
# Create intermediate id string.
id_string = str(p_start) + "-to-" + str(p_end)
# First, generate a list of all the j_long triplestacks.
triplestackfilelist_raw = listdir("./" + obj_string + "_triplestacks")
triplestackfilelist = []
for n in range(len(triplestackfilelist_raw)):
    if ((triplestackfilelist_raw[n][-5:] == ".fits") and 
        (triplestackfilelist_raw[n][:6] == "j_long")):
        triplestackfilelist.append(triplestackfilelist_raw[n])
# Then, clip the list of all p#'s outside the user-specified range.
clippedtriplestackfilelist = []
for filename in triplestackfilelist:
    p_current = int(filename.split("-p")[1].split(".")[0])
    if (p_start <= p_current) and (p_current <= p_end):
        clippedtriplestackfilelist.append("_triplestack" + 
            filename.split("triplestack")[1])
    if p_current == p_start:
        first_filename = filename
    if p_current == p_end:
        last_filename = filename
# Open files to which we'll write the location-paths of the images to be 
# mosaiced.
j_long_progenitors = file("j_long_progenitors.txt", "w")
j_short_progenitors = file("j_short_progenitors.txt", "w")
h_long_progenitors = file("h_long_progenitors.txt", "w")
h_short_progenitors = file("h_short_progenitors.txt", "w")
k_long_progenitors = file("k_long_progenitors.txt", "w")
k_short_progenitors = file("k_short_progenitors.txt", "w")
for filename in clippedtriplestackfilelist:
    j_long_progenitors.write("./" + obj_string + "_triplestacks/j_long" + 
        filename + "\n")
    j_short_progenitors.write("./" + obj_string + "_triplestacks/j_short" + 
        filename + "\n")
    h_long_progenitors.write("./" + obj_string + "_triplestacks/h_long" + 
        filename + "\n")
    h_short_progenitors.write("./" + obj_string + "_triplestacks/h_short" + 
        filename + "\n")
    k_long_progenitors.write("./" + obj_string + "_triplestacks/k_long" + 
        filename + "\n")
    k_short_progenitors.write("./" + obj_string + "_triplestacks/k_short" + 
        filename + "\n")
j_long_progenitors.close()
j_short_progenitors.close()
h_long_progenitors.close()
h_short_progenitors.close()
k_long_progenitors.close()
k_short_progenitors.close()
# We want to propagate the correct exposure time.
num_dither_positions = len(clippedtriplestackfilelist)
mosaic_long_exptime = num_dither_positions * 23.400
mosaic_short_exptime = num_dither_positions * 0.153
# Make list of swarp_commands.
swarp_commands = [
    swarp_bin + " @j_long_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/j_weightmap.fits " + 
    "-IMAGEOUT_NAME j_long_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME j_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits",
    swarp_bin + " @j_short_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/j_weightmap.fits " + 
    "-IMAGEOUT_NAME j_short_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME j_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits",
    swarp_bin + " @h_long_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/h_weightmap.fits " + 
    "-IMAGEOUT_NAME h_long_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME h_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits",
    swarp_bin + " @h_short_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/h_weightmap.fits " + 
    "-IMAGEOUT_NAME h_short_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME h_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits",
    swarp_bin + " @k_long_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/k_weightmap.fits " + 
    "-IMAGEOUT_NAME k_long_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME k_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits",
    swarp_bin + " @k_short_progenitors.txt " + 
    "-c pairitel_redux.swarp " + 
    "-WEIGHT_IMAGE ./" + obj_string + "_weightmaps/k_weightmap.fits " + 
    "-IMAGEOUT_NAME k_short_" + obj_string + "_coadd." + id_string + ".fits " + 
    "-WEIGHTOUT_NAME k_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits"]
# Run the mosaicing with parallel processing.
p = Pool(numprocessors)
result = p.map_async(run_swarp, swarp_commands)
poolresult = result.get()
# Set the exposure time in the mosaics.
system(sethead_bin + " j_long_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_long_exptime))
system(sethead_bin + " j_short_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_short_exptime))
system(sethead_bin + " h_long_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_long_exptime))
system(sethead_bin + " h_short_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_short_exptime))
system(sethead_bin + " k_long_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_long_exptime))
system(sethead_bin + " k_short_" + obj_string + "_coadd." + id_string + 
    ".fits EXPTIME=" + str(mosaic_short_exptime))
# We insert the STRT_CPU of the first triplestack and the STOP_CPU of 
# the last triplestack used to make the mosaic.
system(cphead_bin + " ./" + obj_string + "_triplestacks/" + first_filename + 
    " *coadd*.fits STRT_CPU")
system(cphead_bin + " ./" + obj_string + "_triplestacks/" + last_filename + 
    " *coadd*.fits STOP_CPU")
# And, finally, we insert the STRT_CPU and STOP_CPU for each rawfile which was
# used in the mosaic.
for filename in clippedtriplestackfilelist:
    p_num = int(filename.split("-p")[1].split(".")[0])
    num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
    num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
    num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
    strt_0 = "STRT" + num_0
    strt_1 = "STRT" + num_1
    strt_2 = "STRT" + num_2
    stop_0 = "STOP" + num_0
    stop_1 = "STOP" + num_1
    stop_2 = "STOP" + num_2
    trip_hdulist = pyfits.open("./" + obj_string + 
        "_triplestacks/j_long" + filename)
    trip_header = trip_hdulist[0].header
    strt_0_val = str(trip_header["STRT0000"])
    strt_1_val = str(trip_header["STRT0001"])
    strt_2_val = str(trip_header["STRT0002"])
    stop_0_val = str(trip_header["STOP0000"])
    stop_1_val = str(trip_header["STOP0001"])
    stop_2_val = str(trip_header["STOP0002"])
    trip_hdulist.close()
    system(sethead_bin + " *" + obj_string + "_coadd*.fits "
        "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, 
        strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, 
        stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
# Run the wcs fitting.
print "Now wcs fitting coadd images."
t1 = time()
# Recover the CRVAL1 and CRVAL2 from the fitted triplestack.
p1_hdulist = pyfits.open("./" + obj_string + 
    "_triplestacks/wcsfit-j_long-p1.fits")
p1_header = p1_hdulist[0].header
p1_crval1 = float(p1_header["CRVAL1"])
p1_crval2 = float(p1_header["CRVAL2"])
p1_hdulist.close()
# Copy the CRVAL1 and CRVAL2 which were fit for the p1 triplestacks (via imwcs)
# into the ?_long*coadd.fits mosaics. These will provide the rough astrometric
# solution which scamp will refine.
system(sethead_bin + " j_long_" + obj_string + "_coadd." + id_string + 
    ".fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
system(sethead_bin + " h_long_" + obj_string + "_coadd." + id_string + 
    ".fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
system(sethead_bin + " k_long_" + obj_string + "_coadd." + id_string + 
    ".fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
system(sethead_bin + " j_long_" + obj_string + "_coadd." + id_string + 
    ".weight.fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
system(sethead_bin + " h_long_" + obj_string + "_coadd." + id_string + 
    ".weight.fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
system(sethead_bin + " k_long_" + obj_string + "_coadd." + id_string + 
    ".weight.fits " + "CRVAL1=%f CRVAL2=%f" % (p1_crval1, p1_crval2))
# Next, run SExtractor on each 
system(sextractor_bin + " j_long_" + obj_string + "_coadd." + id_string + 
    ".fits -c pairitel_align.sex -CATALOG_NAME j_long_" + obj_string + 
    "_coadd." + id_string + ".ldac -WEIGHT_IMAGE j_long_" + obj_string + 
    "_coadd." + id_string + ".weight.fits")
system(sextractor_bin + " h_long_" + obj_string + "_coadd." + id_string + 
    ".fits -c pairitel_align.sex -CATALOG_NAME h_long_" + obj_string + 
    "_coadd." + id_string + ".ldac -WEIGHT_IMAGE j_long_" + obj_string + 
    "_coadd." + id_string + ".weight.fits")
system(sextractor_bin + " k_long_" + obj_string + "_coadd." + id_string + 
    ".fits -c pairitel_align.sex -CATALOG_NAME k_long_" + obj_string + 
    "_coadd." + id_string + ".ldac -WEIGHT_IMAGE j_long_" + obj_string + 
    "_coadd." + id_string + ".weight.fits")
# And then, run Scamp.
system(scamp_bin + " j_long_" + obj_string + "_coadd." + id_string + ".ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
system(scamp_bin + " h_long_" + obj_string + "_coadd." + id_string + ".ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
system(scamp_bin + " k_long_" + obj_string + "_coadd." + id_string + ".ldac " + 
    "-c pairitel_align.scamp -ASTREF_CATALOG USNO-B1")
# Finally, use missfits to apply the new WCS header info to the image.
system(missfits_bin + " j_long_" + obj_string + "_coadd." + id_string + ".fits")
system(missfits_bin + " h_long_" + obj_string + "_coadd." + id_string + ".fits") 
system(missfits_bin + " k_long_" + obj_string + "_coadd." + id_string + ".fits") 
# Finish up by copying the WCS header info from the ?_long*coadd.fits mosaics 
# into the corresponding weightmaps and ?_short*coadd.fits mosaics and 
# weightmaps.
system("cp j_long_" + obj_string + "_coadd." + id_string + ".head j_long_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system("cp j_long_" + obj_string + "_coadd." + id_string + ".head j_short_" + 
    obj_string + "_coadd." + id_string + ".head")
system("cp j_long_" + obj_string + "_coadd." + id_string + ".head j_short_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system("cp h_long_" + obj_string + "_coadd." + id_string + ".head h_long_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system("cp h_long_" + obj_string + "_coadd." + id_string + ".head h_short_" + 
    obj_string + "_coadd." + id_string + ".head")
system("cp h_long_" + obj_string + "_coadd." + id_string + ".head h_short_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system("cp k_long_" + obj_string + "_coadd." + id_string + ".head k_long_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system("cp k_long_" + obj_string + "_coadd." + id_string + ".head k_short_" + 
    obj_string + "_coadd." + id_string + ".head")
system("cp k_long_" + obj_string + "_coadd." + id_string + ".head k_short_" + 
    obj_string + "_coadd." + id_string + ".weight.head")
system(missfits_bin + " j_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits")
system(missfits_bin + " h_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits") 
system(missfits_bin + " k_long_" + obj_string + "_coadd." + 
    id_string + ".weight.fits")
system(missfits_bin + " j_short_" + obj_string + "_coadd." + 
    id_string + ".fits")
system(missfits_bin + " h_short_" + obj_string + "_coadd." + 
    id_string + ".fits") 
system(missfits_bin + " k_short_" + obj_string + "_coadd." + 
    id_string + ".fits")
system(missfits_bin + " j_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits")
system(missfits_bin + " h_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits") 
system(missfits_bin + " k_short_" + obj_string + "_coadd." + 
    id_string + ".weight.fits")
# Clean up the directory of all the intermediate files.
system("rm missfits.xml")
system("rm *.back")
system("rm *.ldac")
system("rm *.head")
system("rm *_progenitors.txt")
system("rm -rf ./" + obj_string + "_intermediate_mosaics_" + id_string)
system("mkdir ./" + obj_string + "_intermediate_mosaics_" + id_string)
system("mv -f *" + obj_string + "_coadd*fits ./" + obj_string + 
    "_intermediate_mosaics_" + id_string)
# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time