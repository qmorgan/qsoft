"""
mosaic_maker.py mosaics your PAIRITEL data after full pipeline3 reduction with 
manual control over input triplestack files. You have to run mosaic_maker.py
in the same directory that contains the reduction_output directory from 
pipeline3.

You will also need the following files:
    mosaic_maker.py             -   This file, obviously
    anet.py                     -   Python code to use astrometry.net
    pairitel_redux.swarp        -   SWarp configuration file for mosacing

Dependencies:
    SWarp
    WCSTools (sethead)
    Python 2.5 or better
    pyfits (python module)

Example Usage: 
    python mosaic_maker.py -o PULSE.33.1 -p
        [Then user edits the newly created j_long_triplestacks.txt, 
        h_long_triplestacks.txt, and k_long_triplestacks.txt files to exclude 
        unwanted triplestacks from the final mosaics]
    python mosaic_maker.py -o PULSE.33.1 
        [This will use the recently edited text files to make mosaics]
        
Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/09/04
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
import pyfits
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
swarp_bin = "/usr/local/bin/swarp"
# We'll also be using some WCSTools. There should be no problem with the Scisoft
# version of thses.
sethead_bin = "/usr/bin/sethead"
# For parallel processing we need the number of available CPUs. You could hard-
# code this value to something else, but it would result in non-optimal 
# performance.

# AMORGAN ADDS for WCS FITTING  Change if your python version is different.
# Will want to avoid having to load this in future versions.
pypath = "~/Programs/epd-6.1-1-rh5-x86/bin/"

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
parser = OptionParser()

parser.add_option("-o", "--obs", "--obs-id", "--obs-string", action="store",     
    type="string", dest="obs_string", help=("observation ID string of epoch " + 
    "to be reduced")) 
parser.add_option("-s", "--short", "--shorts",
                  action="store_true", dest="do_short", default=False,
                  help="reduce short read data")
parser.add_option("-p", "--prep",
                  action="store_true", dest="do_prep", default=False,
                  help=("generate SWarp input lists to edit manually before " + 
                  "running in mosaicing mode"))
# AMorgan adds option
parser.add_option("-w", "--wcs", "--do-wcs",
                  action="store_true", dest="do_wcs", default=False,
                  help=("Attempt to do wcs fitting after coaddition"))

# Option for use in single reduced image
parser.add_option("-r", "--single",
                  action="store_true", dest="do_single", default=False,
                  help=("Use this for single reduced image"))

(options, args) = parser.parse_args()

obs_string = options.obs_string
if not obs_string:
    print "No obs_string specified, exiting."
    sys.exit()
do_short = options.do_short
do_prep = options.do_prep
do_wcs = options.do_wcs
do_single = options.do_single

reduction_output_directory = str(obs_string) + "-reduction_output"
triplestacks_path = (reduction_output_directory + "/" + str(obs_string) + 
    "_triplestacks")
triplestackweights_path = (reduction_output_directory + "/" + str(obs_string) + 
    "_triplestackweights")

if do_prep:
    ancillary_path = (reduction_output_directory + "/" + str(obs_string) + 
        "_ancillary")
    j_long_list = file(ancillary_path + "/j_long_triplestacks.txt", "r")
    j_long_list_new = file("j_long_triplestacks.txt", "w")
    for line in j_long_list:
        j_long_list_new.write(line.replace("triplestacks/", triplestacks_path +
            "/"))
    j_long_list.close()
    j_long_list_new.close()
    h_long_list = file(ancillary_path + "/h_long_triplestacks.txt", "r")
    h_long_list_new = file("h_long_triplestacks.txt", "w")
    for line in h_long_list:
        h_long_list_new.write(line.replace("triplestacks/", triplestacks_path +
            "/"))
    h_long_list.close()
    h_long_list_new.close()
    k_long_list = file(ancillary_path + "/k_long_triplestacks.txt", "r")
    k_long_list_new = file("k_long_triplestacks.txt", "w")
    for line in k_long_list:
        k_long_list_new.write(line.replace("triplestacks/", triplestacks_path +
            "/"))
    k_long_list.close()
    k_long_list_new.close()    
    if do_short:
        j_short_list = file(ancillary_path + "/j_short_triplestacks.txt", "r")
        j_short_list_new = file("j_short_triplestacks.txt", "w")
        for line in j_short_list:
            j_short_list_new.write(line.replace("triplestacks/", 
                triplestacks_path + "/"))
        j_short_list.close()
        j_short_list_new.close()
        h_short_list = file(ancillary_path + "/h_short_triplestacks.txt", "r")
        h_short_list_new = file("h_short_triplestacks.txt", "w")
        for line in h_short_list:
            h_short_list_new.write(line.replace("triplestacks/", 
                triplestacks_path + "/"))
        h_short_list.close()
        h_short_list_new.close()
        k_short_list = file(ancillary_path + "/k_short_triplestacks.txt", "r")
        k_short_list_new = file("k_short_triplestacks.txt", "w")
        for line in k_short_list:
            k_short_list_new.write(line.replace("triplestacks/", 
                triplestacks_path + "/"))
        k_short_list.close()
        k_short_list_new.close()
    sys.exit()

j_long_triplestacks_list = []
j_long_list = file("j_long_triplestacks.txt", "r")
j_long_list_weights = file("j_long_triplestackweights.txt", "w")

if not do_single:
    replaced_str = '_triplestack'
else:
    replaced_str = '_reduced'

for line in j_long_list:
    j_long_triplestacks_list.append(line.rstrip())
    j_long_list_weights.write(line.replace(replaced_str,
        "_triplestackweightmap").replace("_triplestacks/", 
        "_triplestackweights/").replace("weightmaps",
        "weights"))
j_long_list.close()
j_long_list_weights.close()
h_long_triplestacks_list = []
h_long_list = file("h_long_triplestacks.txt", "r")
h_long_list_weights = file("h_long_triplestackweights.txt", "w")
for line in h_long_list:
    h_long_triplestacks_list.append(line.rstrip())
    h_long_list_weights.write(line.replace(replaced_str,
        "_triplestackweightmap").replace("_triplestacks/", 
        "_triplestackweights/").replace("weightmaps",
        "weights"))
h_long_list.close()
h_long_list_weights.close()
k_long_triplestacks_list = []
k_long_list = file("k_long_triplestacks.txt", "r")
k_long_list_weights = file("k_long_triplestackweights.txt", "w")
for line in k_long_list:
    k_long_triplestacks_list.append(line.rstrip())
    k_long_list_weights.write(line.replace(replaced_str,
        "_triplestackweightmap").replace("_triplestacks/", 
        "_triplestackweights/").replace("weightmaps",
        "weights"))
k_long_list.close()
k_long_list_weights.close()

if do_short:
    j_short_triplestacks_list = []
    j_short_list = file("j_short_triplestacks.txt", "r")
    j_short_list_weights = file("j_short_triplestackweights.txt", "w")
    for line in j_short_list:
        j_short_triplestacks_list.append(line.rstrip())
        j_short_list_weights.write(line.replace("_triplestack",
            "_triplestackweightmap").replace("_triplestacks/", 
            "_triplestackweights/").replace("weightmaps",
            "weights"))
    j_short_list.close()
    j_short_list_weights.close()
    h_short_triplestacks_list = []
    h_short_list = file("h_short_triplestacks.txt", "r")
    h_short_list_weights = file("h_short_triplestackweights.txt", "w")
    for line in h_short_list:
        h_short_triplestacks_list.append(line.rstrip())
        h_short_list_weights.write(line.replace("_triplestack",
            "_triplestackweightmap").replace("_triplestacks/", 
            "_triplestackweights/").replace("weightmaps",
            "weights"))
    h_short_list.close()
    h_short_list_weights.close()
    k_short_triplestacks_list = []
    k_short_list = file("k_short_triplestacks.txt", "r")
    k_short_list_weights = file("k_short_triplestackweights.txt", "w")
    for line in k_short_list:
        k_short_triplestacks_list.append(line.rstrip())
        k_short_list_weights.write(line.replace("_triplestack",
            "_triplestackweightmap").replace("_triplestacks/", 
            "_triplestackweights/").replace("weightmaps",
            "weights"))
    k_short_list.close()
    k_short_list_weights.close()


# Run the mosaicing. 
# Make list of swarp_commands.
if do_short:
    swarp_commands = [
        swarp_bin + " @j_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
        "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @j_short_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @j_short_triplestackweights.txt " + 
        "-IMAGEOUT_NAME j_short_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME j_short_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @h_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
        "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @h_short_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @h_short_triplestackweights.txt " + 
        "-IMAGEOUT_NAME h_short_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME h_short_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @k_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
        "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @k_short_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @k_short_triplestackweights.txt " + 
        "-IMAGEOUT_NAME k_short_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME k_short_" + obs_string + "_coadd.weight.fits"]
if not do_short:
    swarp_commands = [
        swarp_bin + " @j_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
        "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @h_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
        "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
        "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits",
        swarp_bin + " @k_long_triplestacks.txt " + 
        "-c pairitel_redux.swarp " + 
        "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
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
trip_hdulist = pyfits.open(j_long_triplestacks_list[0])
trip_header = trip_hdulist[0].header
strt_cpu = str(trip_header["STRT_CPU"])
trip_hdulist.close()
trip_hdulist = pyfits.open(j_long_triplestacks_list[-1])
trip_header = trip_hdulist[0].header
stop_cpu = str(trip_header["STOP_CPU"])
trip_hdulist.close()
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits " + 
    "j_long_" + obs_string + "_coadd.weight.fits " + 
    "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
trip_hdulist = pyfits.open(h_long_triplestacks_list[0])
trip_header = trip_hdulist[0].header
strt_cpu = str(trip_header["STRT_CPU"])
trip_hdulist.close()
trip_hdulist = pyfits.open(h_long_triplestacks_list[-1])
trip_header = trip_hdulist[0].header
stop_cpu = str(trip_header["STOP_CPU"])
trip_hdulist.close()
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits " + 
    "h_long_" + obs_string + "_coadd.weight.fits " + 
    "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
trip_hdulist = pyfits.open(k_long_triplestacks_list[0])
trip_header = trip_hdulist[0].header
strt_cpu = str(trip_header["STRT_CPU"])
trip_hdulist.close()
trip_hdulist = pyfits.open(k_long_triplestacks_list[-1])
trip_header = trip_hdulist[0].header
stop_cpu = str(trip_header["STOP_CPU"])
trip_hdulist.close()
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits " + 
    "k_long_" + obs_string + "_coadd.weight.fits " + 
    "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
if do_short:
    trip_hdulist = pyfits.open(j_short_triplestacks_list[0])
    trip_header = trip_hdulist[0].header
    strt_cpu = str(trip_header["STRT_CPU"])
    trip_hdulist.close()
    trip_hdulist = pyfits.open(j_short_triplestacks_list[-1])
    trip_header = trip_hdulist[0].header
    stop_cpu = str(trip_header["STOP_CPU"])
    trip_hdulist.close()
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits " + 
        "j_short_" + obs_string + "_coadd.weight.fits " + 
        "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
    trip_hdulist = pyfits.open(h_short_triplestacks_list[0])
    trip_header = trip_hdulist[0].header
    strt_cpu = str(trip_header["STRT_CPU"])
    trip_hdulist.close()
    trip_hdulist = pyfits.open(h_short_triplestacks_list[-1])
    trip_header = trip_hdulist[0].header
    stop_cpu = str(trip_header["STOP_CPU"])
    trip_hdulist.close()
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits " + 
        "h_short_" + obs_string + "_coadd.weight.fits " + 
        "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
    trip_hdulist = pyfits.open(k_short_triplestacks_list[0])
    trip_header = trip_hdulist[0].header
    strt_cpu = str(trip_header["STRT_CPU"])
    trip_hdulist.close()
    trip_hdulist = pyfits.open(k_short_triplestacks_list[-1])
    trip_header = trip_hdulist[0].header
    stop_cpu = str(trip_header["STOP_CPU"])
    trip_hdulist.close()
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits " + 
        "k_short_" + obs_string + "_coadd.weight.fits " + 
        "STOP_CPU='%s'" % (strt_cpu, stop_cpu))

# And, finally, we insert the STRT_CPU and STOP_CPU for each rawfile which was
# used in the mosaic.

def strt_stop(j_long_triplestacks_list, h_long_triplestacks_list, k_long_triplestacks_list):

    for triplestackfile in j_long_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " j_long_" + obs_string + "_coadd.fits " + "j_long_" + obs_string + "_coadd.weight.fits " + "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0,strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
    for triplestackfile in h_long_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " h_long_" + obs_string + "_coadd.fits " + "h_long_" + obs_string + "_coadd.weight.fits " + "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
    for triplestackfile in k_long_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " k_long_" + obs_string + "_coadd.fits " + "k_long_" + obs_string + "_coadd.weight.fits " + "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))

#Temporarily make it so this code only work with single files

if not do_single:
    strt_stop(j_long_triplestacks_list, h_long_triplestacks_list, k_long_triplestacks_list)
else:
    pass

if do_short:
    for triplestackfile in j_short_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " j_short_" + obs_string + "_coadd.fits " + 
            "j_short_" + obs_string + "_coadd.weight.fits " + 
            "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, 
            strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, 
            stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
    for triplestackfile in h_short_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " h_short_" + obs_string + "_coadd.fits " + 
            "h_short_" + obs_string + "_coadd.weight.fits " + 
            "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, 
            strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, 
            stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
    for triplestackfile in k_short_triplestacks_list:
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open(triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " k_short_" + obs_string + "_coadd.fits " + 
            "k_short_" + obs_string + "_coadd.weight.fits " + 
            "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, 
            strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, 
            stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))
# AMORGAN ADDS OPTION
if do_wcs:            
    # Run the wcs fitting.
    # NOTE TO SELF.  Change this so that pypath will not need to be known.
    # But the simple hack should work for now.
    system(pypath + 
        "python anet.py *_long_" + obs_string + "_coadd.fits")
    j_long_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits", 
        "readonly")
    j_long_weights_hdulist = pyfits.open("j_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        j_short_hdulist = pyfits.open("j_short_" + obs_string + "_coadd.fits", 
            "update")
        j_short_weights_hdulist = pyfits.open("j_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp0 = j_long_hdulist[0].header
    try:
        j_long_an_jobid = tmp0["AN_JOBID"]
    except:
        j_long_an_jobid = False
    del tmp0["SIMPLE"]
    del tmp0["BITPIX"]
    del tmp0["NAXIS"]
    j_long_hdulist.close()
    h_long_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits", 
        "readonly")
    h_long_weights_hdulist = pyfits.open("h_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        h_short_hdulist = pyfits.open("h_short_" + obs_string + "_coadd.fits", 
            "update")
        h_short_weights_hdulist = pyfits.open("h_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp1 = h_long_hdulist[0].header
    try:
        h_long_an_jobid = tmp1["AN_JOBID"]
    except:
        h_long_an_jobid = False
    del tmp1["SIMPLE"]
    del tmp1["BITPIX"]
    del tmp1["NAXIS"]
    h_long_hdulist.close()
    k_long_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits", 
        "readonly")
    k_long_weights_hdulist = pyfits.open("k_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        k_short_hdulist = pyfits.open("k_short_" + obs_string + "_coadd.fits", 
            "update")
        k_short_weights_hdulist = pyfits.open("k_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp2 = k_long_hdulist[0].header
    try:
        k_long_an_jobid = tmp2["AN_JOBID"]
    except:
        k_long_an_jobid = False
    del tmp2["SIMPLE"]
    del tmp2["BITPIX"]
    del tmp2["NAXIS"]
    k_long_hdulist.close()
    if j_long_an_jobid:
        for c in tmp0.ascardlist():
            j_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                j_short_hdulist[0].header.update(c.key,c.value,c.comment)
                j_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if h_long_an_jobid:
        for c in tmp1.ascardlist():
            h_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                h_short_hdulist[0].header.update(c.key,c.value,c.comment)
                h_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if k_long_an_jobid:
        for c in tmp2.ascardlist():
            k_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                k_short_hdulist[0].header.update(c.key,c.value,c.comment)
                k_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if j_long_an_jobid and not h_long_an_jobid:
        h_long_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits", 
        "update")
        for c in tmp0.ascardlist():
            if c.key == "FILTER":
                continue
            if c.key == "CRPIX1":
                h_long_hdulist[0].header.update(c.key,float(c.value) - 
                    2.492148,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    2.492148,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,float(c.value) - 
                        2.492148,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        2.492148,c.comment)
            elif c.key == "CRPIX2":
                h_long_hdulist[0].header.update(c.key,float(c.value) - 
                    0.0969426,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    0.0969426,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,float(c.value) - 
                        0.0969426,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        0.0969426,c.comment)
            else:
                h_long_hdulist[0].header.update(c.key,c.value,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,c.value,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,c.value,
                        c.comment)
        h_long_hdulist.verify("silentfix")
        h_long_hdulist.close(output_verify='warn')  
    if j_long_an_jobid and not k_long_an_jobid:
        k_long_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits", 
        "update")
        for c in tmp0.ascardlist():
            if c.key == "FILTER":
                continue
            elif c.key == "CRPIX1":
                k_long_hdulist[0].header.update(c.key,float(c.value) - 
                    0.5982118,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    0.5982118,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,float(c.value) - 
                        0.5982118,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        0.5982118,c.comment)
            elif c.key == "CRPIX2":
                k_long_hdulist[0].header.update(c.key,float(c.value) + 
                    6.07977,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,float(c.value) + 
                    6.07977,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,float(c.value) + 
                        6.07977,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,float(c.value) + 
                        6.07977,c.comment)
            elif c.key == "CD1_2":
                k_long_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                    c.comment)
                k_long_weights_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                    c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                        c.comment)
                    k_short_weights_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                        c.comment)
            elif c.key == "CD2_1":
                k_long_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                    c.comment)
                k_long_weights_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                    c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                        c.comment)
                    k_short_weights_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                        c.comment)
            else:
                k_long_hdulist[0].header.update(c.key,c.value,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,c.value,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,c.value,
                        c.comment)
        k_long_hdulist.verify("silentfix")
        k_long_hdulist.close(output_verify='warn')
    j_long_weights_hdulist.verify("silentfix")
    j_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        j_short_hdulist.verify("silentfix")
        j_short_hdulist.close(output_verify='warn')
        j_short_weights_hdulist.verify("silentfix")
        j_short_weights_hdulist.close(output_verify='warn')
    h_long_weights_hdulist.verify("silentfix")
    h_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        h_short_hdulist.verify("silentfix")
        h_short_hdulist.close(output_verify='warn')
        h_short_weights_hdulist.verify("silentfix")
        h_short_weights_hdulist.close(output_verify='warn')
    k_long_weights_hdulist.verify("silentfix")
    k_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        k_short_hdulist.verify("silentfix")
        k_short_hdulist.close(output_verify='warn')
        k_short_weights_hdulist.verify("silentfix")
        k_short_weights_hdulist.close(output_verify='warn')

# We want to propagate the correct exposure time.
j_long_num_dither_positions = len(j_long_triplestacks_list)
j_long_mosaic_exptime = j_long_num_dither_positions * 23.400
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(j_long_mosaic_exptime))
h_long_num_dither_positions = len(h_long_triplestacks_list)
h_long_mosaic_exptime = h_long_num_dither_positions * 23.400
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(h_long_mosaic_exptime))
k_long_num_dither_positions = len(k_long_triplestacks_list)
k_long_mosaic_exptime = k_long_num_dither_positions * 23.400
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(k_long_mosaic_exptime))
if do_short:
    j_short_num_dither_positions = len(j_short_triplestacks_list)
    j_short_mosaic_exptime = j_short_num_dither_positions * 23.400
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(j_short_mosaic_exptime))
    h_short_num_dither_positions = len(h_short_triplestacks_list)
    h_short_mosaic_exptime = h_short_num_dither_positions * 23.400
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(h_short_mosaic_exptime))
    k_short_num_dither_positions = len(k_short_triplestacks_list)
    k_short_mosaic_exptime = k_short_num_dither_positions * 23.400
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(k_short_mosaic_exptime))

# Clean up the directory of all the intermediate files.
system("rm *_triplestackweights.txt")

# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time
