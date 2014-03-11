"""
Photometry Package for PAIRITEL Data

The following files are required in the $Q_DIR/load directory:
    pairitel_photo.sex
    pairitel_photo.param
    default.nnw
    star_find.param

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Authors: Christopher Klein & Adam Morgan
Contact: cklein@berkeley.edu, amorgan@berkeley.edu
"""


# --------------------------    LIBRARY and MODULE IMPORTS  --------------------

import sys
import ephem
import pylab
from pylab import *
import os
from os import system
from math import log
from scipy import average, std, median, random, zeros, clip, array, append
from scipy import histogram, optimize, sqrt, pi, exp, where, split
from scipy import sum, std, resize, size, isnan
from scipy import optimize, shape, indices, arange, ravel, sqrt
from scipy.ndimage import gaussian_filter, binary_dilation
from scipy import loadtxt, savetxt, sort
import pyfits
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import operator
from Phot import t_mid
from MiscBin import qPickle
from pylab import close
from time import time
import numpy

# --------------------------    USER INPUT PARAMETERS   --------------------

# progenitor_image_name = "j_long_laptop_coadd.fits"
# region_file = "ds9.reg"
# progenitor_image_name = "j_long_GRB.10928.1_coadd_1-1.fits"
# region_file = "PTEL.reg"
# 
# sextractor_bin = "/opt/local/bin/sex"
# weight_image_name = progenitor_image_name.replace("fits", "weight.fits")
# 
# 

# --------------------------    FUNCTION DEFINITIONS    ------------------------

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'
sextractor_bin = "sex"

#############
# NOTE: Make sure the path is set correctly in the keyword
# PARAMETERS_NAME in $Q_DIR/load/pairitel_photo.sex
#############

# Calculate the residual for the gaussian fit.
def gaussian_residuals(params, y, x):
    central_height, center_x, x_width = params
    err = y - ((central_height/sqrt(2*pi*x_width**2))*
        exp((-1*(x-center_x)**2)/(2*x_width**2)))
    return err
# Fit the input data with a gaussian curve.
def fit_gaussian(y, x):
    input_params = [5000.0, 0.0, 30.0]
    fit_params, success = optimize.leastsq(gaussian_residuals, input_params, 
        args=(y, x), maxfev=10000)
    return fit_params
# Call Source Extractor in two-image mode.
def make_sex_cat_fake_stars(fake_image, progenitor_image_name, 
    weight_image_name, aperture_size):
    min_pix = pi*(aperture_size*aperture_size/4)/4
    sexcat_file = progenitor_image_name.replace(".fits", ".sex")
    system(sextractor_bin + " " + fake_image + ", " + progenitor_image_name + 
        " -c " + loadpath + "pairitel_photo.sex " + 
        "-DETECT_MINAREA " + str(min_pix) +
        " -BACK_SIZE 64 " + 
        "-BACK_FILTERSIZE 3 " + 
        "-BACK_TYPE AUTO " + 
        "-BACK_VALUE 0.0 " +
        "-FILTER N " + 
        "-PHOT_APERTURES " + str(aperture_size) + 
        " -WEIGHT_TYPE NONE" +
        " -WEIGHT_IMAGE " + weight_image_name + 
        " -CATALOG_NAME " + sexcat_file)
    return sexcat_file
# Call Source Extractor in two-image mode with no weight image.
def make_sex_cat_fake_stars_no_weight(fake_image, progenitor_image_name, aperture_size):
    sexcat_file = progenitor_image_name.replace(".fits", ".sex")
    min_pix = pi*(aperture_size*aperture_size/4)/4
    system(sextractor_bin + " " + fake_image + ", " + progenitor_image_name + 
        " -c " + loadpath + "pairitel_photo.sex " + 
        "-DETECT_MINAREA " + str(min_pix) +
        " -BACK_SIZE 64 " + 
        "-BACK_FILTERSIZE 3 " + 
        "-BACK_TYPE AUTO " + 
        "-BACK_VALUE 0.0 " +
        "-FILTER N " + 
        "-BACK_TYPE MANUAL " +
        "-PHOT_APERTURES " + str(aperture_size) + 
        " -WEIGHT_TYPE NONE" + 
        " -CATALOG_NAME " + sexcat_file)
    return sexcat_file
# Call Source Extractor for photometry on science image.
def make_sex_cat(progenitor_image_name, weight_image_name, aperture_size):
    sexcat_file = progenitor_image_name.replace(".fits", ".sex")
    min_pix = pi*(aperture_size*aperture_size/4)/4
    system(sextractor_bin + " " + progenitor_image_name + " -c " + loadpath + "pairitel_photo.sex " + 
        "-DETECT_MINAREA " + str(min_pix) +
        " -BACK_SIZE 64 " + 
        "-BACK_FILTERSIZE 3 " + 
        "-BACK_TYPE AUTO " + 
        "-BACK_VALUE 0.0 " +
        "-FILTER N " +
        "-PHOT_APERTURES " + str(aperture_size) + 
        " -WEIGHT_IMAGE " + weight_image_name + 
        " -CATALOG_NAME " + sexcat_file)
    return sexcat_file
# Calculate the heliocentric julian date.
def heliocentric_julian_date(observation_time, observation_ra_radians, 
    observation_dec_radians):
    # Compute the observation time in Heliocentric Julian Date. First convert to 
    # julian date (jd) by adding the offset to ephem.date format.
    observation_time_jd = float(observation_time) + 2415020.0
    # Calculate the Sun's position in the sky.
    sun = ephem.Sun()
    sun.compute(observation_time)
    sun_ra_radians = float(sun.a_ra)
    sun_dec_radians = float(sun.a_dec)
    # Calculate the Earth-Sun light travel time in days.
    earth_sun_light_travel_time = sun.earth_distance * 0.00577551833
    # Calculate the observation time in heliocentric julian date.
    observation_time_hjd = (observation_time_jd - earth_sun_light_travel_time * 
        (sin(observation_dec_radians) * sin(sun_ra_radians) + 
        cos(observation_dec_radians) * cos(sun_ra_radians) * 
        cos(observation_dec_radians - sun_dec_radians)))
    return observation_time_hjd
# Median Absolute Deviation clipping for input list of numbers.
def mad_clipping(input_data, sigma_clip_level):
    medval = median(input_data)
    sigma = 1.48 * median(abs(medval - input_data))
    high_sigma_clip_limit = medval + sigma_clip_level * sigma
    low_sigma_clip_limit = medval - sigma_clip_level * sigma
    clipped_data = []
    for value in input_data:
        if (value > low_sigma_clip_limit) and (value < high_sigma_clip_limit):
            clipped_data.append(value)
    clipped_data_array = array(clipped_data)
    new_medval = median(clipped_data_array)
    new_sigma = 1.48 * median(abs(medval - clipped_data_array))
    return clipped_data_array, new_medval, new_sigma

def moffat_residuals(params, z, x, y):
    central_height, center_x, center_y, alpha, beta = params
    r2 = (x - center_x)**2 + (y - center_y)**2
    err = z - (central_height * (1 + (r2/(alpha**2)))**(-beta))
    return err
    

def fit_moffat(z, x, y):
    """Returns (central_height, center_x, center_y, alpha, beta
    the moffat parameters of a 2D distribution found by a fit"""
    input_params = [15000.0, 5.0, 5.0, 5.0, 5.0]
    fit_params, success = optimize.leastsq(moffat_residuals, input_params, 
        args=(z, x, y), maxfev=10000)
    return fit_params
    
def fit_fwhm(sat_locations, objects_data, fwhm, fwhm_stdev):
    fwhm_list = []
    for location in sat_locations:
        m = location[0] # vertical from bottom
        n = location[1] # horizontal from left
        submask_data = objects_data[m-5:m+5,n-5:n+5]
        z = []
        x = []
        y = []
        x_y = []
        r = []
        for i in range(len(submask_data)):
            i_pixel = i + 0
            for j in range(len(submask_data[i])):
                if submask_data[i][j] != 0:
                    j_pixel = j + 0
                    x.append(j_pixel)
                    y.append(i_pixel)
                    x_y.append([j_pixel, i_pixel])
                    z.append(submask_data[i][j])
        try:
            fit_params = fit_moffat(z, x, y)
        except(TypeError):
            continue
        central_height, center_x, center_y, alpha, beta = fit_params
        for coordinate in x_y:
            r.append(sqrt((coordinate[0]-center_x)**2 + 
                (coordinate[1]-center_y)**2))
        fit_data = []
        residual_data = []
        abs_residual_data = []
        for k in range(len(r)):
            fit_point = (central_height * (1 + (r[k]/alpha)**2)**(-beta))
            fit_data.append(fit_point)
            residual_data.append(z[k] - fit_point)
            abs_residual_data.append(abs(z[k] - fit_point))
    
        metric = mean(abs_residual_data) / central_height # Lower is better
        fwhm_unqiue = 2*alpha*sqrt(2**(1/beta) - 1)
        if ((metric < 0.02) and (fwhm_unqiue > (fwhm - fwhm_stdev)) and 
            (fwhm_unqiue < (fwhm + fwhm_stdev))):
#             print ("Finished fit for saturated star at " + str(m) + " " 
#                 + str(n) + ". FWHM: " + str(fwhm_unqiue) + " with metric: " + str(metric))
            fwhm_list.append(fwhm_unqiue)
    fwhm = median(fwhm_list)
    fwhm_stdev = std(fwhm_list)
    return fwhm, fwhm_stdev
    

def find_number_of_dither_positions(progenitor_image_name,fallback_n_dither=0):
    '''Finding the number of triplestacks (dither positions) used in the image
    for use when calulating the uncertainty (an addl source of error occurs
    which reduces as the sqrt of the number of dither positions):
    
    This works by counting the number of images used to make each file, 
    indicated by the number of STRT0001 type instances in the header.
    
    For some manually construceted mosaics these keywords are lost, so the number
    of dither positions will need to be provided either manually (through
    raw input), or via the fallback_n_dither keyword which will assign this 
    value.
    '''
    file_header = pyfits.open(progenitor_image_name)
    STRT_count = 0
    for info in file_header[0].header:
        if 'STRT' in info:
            STRT_count += 1
    n_ind_images = STRT_count - 1 # accounting for the STRT_CPU
    print 'Found ' + str(n_ind_images) + ' individual images'
    # if n_ind_images == 0:
    #     num_triplestacks = 1
    # elif n_ind_images == 1:
    #     num_triplestacks = 1
    # else:
    #     num_triplestacks = n_ind_images/3
    # if n_ind_images == 0:
    #     n_ind_images = 1
    if n_ind_images > 2:
        print "Assuming the number of dither positions (triplestacks) is 1/3rd the number of images."
        num_triplestacks = n_ind_images/3
    elif fallback_n_dither > 0:
        print "\n\nCould not auto-detect number of dither positions! Using fallback value of %s\n\n" % (str(fallback_n_dither))
        num_triplestacks = int(fallback_n_dither)
    else:
        print """Unknown number of dither positions! Needed to estimate PAIRITEL 
uncertainty. (This is a negligible source of error for large mosaics, 
so just enter a large number of order 100/hour if this is a mosaic)."""
        print progenitor_image_name
        num_triplestacks_str = raw_input("Enter number of dither positions: ")
        num_triplestacks = int(num_triplestacks_str)
    if n_ind_images%3 != 0:
        warning = '\n\n\nWarning: Number of individual images is not divisible by 3! Was this mosaic made of triplestacks?\n\n\n'
        print warning
    print "Determined %i input images and %i dither positions" % (n_ind_images,num_triplestacks)
    ditherdict = {'N_dither_pos':num_triplestacks,'N_input_images':n_ind_images}
    return ditherdict
    

def return_ptel_uncertainty(instrumental_error,zeropoint_error,num_triplestacks,
                            base_dither_error=0.03):
    '''
    Dan tested PAIRITEL uncertainty sources extensively for his work on 
    GRB 071025 (Perley et al. 2010, Sec 2.2).
    
    In these mosaics we are resampling (with the Swarp program) from
    PAIRITEL's native pixel scale of 2"/pix to 1"/pix; since we're
    drizzling from 1 into 4 fake pixels, the nominal statistical
    uncertainty is under-sampled by a factor of 2 (actually found to be
    about 2.4).   The instrumental error is thus multiplied by this factor.
    
    The large plate scale of PAIRITEL (2.0arcsec/pix) and the variable sub-pixel 
    response function of the NICMOS3 arrays creates a significant additional 
    uncertainty in each position beyond ordinary photometric errors, estimated 
    at ~3 percent by Blake et al. (2008). To quantify this uncertainty as 
    accurately as possible, we constructed light curves for standard stars
    of different magnitudes in regions of the image free of defects by 
    measuring the image-to-image magnitude variations of bright 
    (source-dominated) stars. An additional uncertainty of approximately 
    ~0.02 mag per position was required to incorporate the observed scatter 
    in the photometry of these objects.
    
    email 5/18/11 to pierre
    ////Grand! Very promising results.  This clearly shows that we've been
    overestimating the per-dither uncertainty in the photometry code
    (currently set to 0.1 mag).  At some point we should verify these
    values with well-detected sources in other fields, but it appears that
    we're looking at values of 0.02 for J, 0.025 for H, and 0.03 for K.
    ////
    
    The error in zeropoint is added in quadrature with these other two sources
    of uncertainty to determine the final uncertainty.
    
    '''
    uncertainty = float(sqrt(zeropoint_error**2 + (instrumental_error*2.4)**2 
                    + (base_dither_error/sqrt(num_triplestacks))**2))
    return uncertainty

# --------------------------    BEGIN PROGRAM   --------------------------------

def dophot(progenitor_image_name,region_file, ap=None, find_fwhm = False, \
        do_upper = False, calstar_reg_output = True, calreg = None, stardict = None,
        cleanup=True, caliblimit=True, verbose=False, autocull=False, for_wcs_coadd=False,
        fallback_n_dither=0, reuse_old_vizcat=False):
    '''
    Do photometry on a given image file given a region file with the coordinates.
    
    Arguments:
        progenitor_image_name: the path to the image file
        region_file: the region file containing a single fk5 coordinate position
            on which to perform photometry.  The associated positional uncertainty
            within this file does not matter.
        ap: desired aperature in arcseconds.  If not specified, we will try 
            to find_fwhm to get the optimal aperture.  If this does not work,
            a default of 4.5 arcseconds will be used
        find_fwhm: find the fwhm of the image.  If ap is not set, then this will
            also be used to attempt to find the optimal aperture, set as
            1.5*FWHM
        do_upper: If False, skip the steps to attempt to find the limiting 
            magnitude of the field.  As this is a very time consuming step,
            it is beneficial to skip if a detection is assured.
        calstar_reg_output: if True, save calibration stars that are found by 
            the code to be usable into a region file, which can be trimmed 
            and fed to calreg
        calreg: If the path to a calibration star region file is provided 
            (created and modified from calstar_reg_output), use only those 
            targets which are BOTH found by the code to be useable AND in
            the calreg file.  Currently no implementation to override what
            the code thinks is usuable calibration stars.
        stardict: If set, use a dictionary of calibration star photometries 
            from a deep stack of the location rather than the 2mass values.  
            CURRENTLY NOT WORKING.
        caliblimit: if True, then if a calibration star has a higher stdev than
            our threshold, then stop the program
        fallback_n_dither: fallback value if the number of dither positions 
            cannot be determined
        reuse_old_vizcat: if True, then do not make the call to download a new
            vizcat and just use the old one in the storepath.
        
    '''
    # Begin timing
    t1 = time()
    
    if not calreg:
        caliblimit = False 
    
    # Store the original image name.
    image_name = progenitor_image_name
    hdulist = pyfits.open(image_name)
    image_data = hdulist[0].data
    hdulist.close()
    photdict = {'FileName':image_name}
    
    weight_image_name = progenitor_image_name.replace("fits", "weight.fits")
    
    # Run source extractor to find the pixel locations of stars.
    print "Running Source Extractor to find locations of stars."
    syscmd = sextractor_bin + " -c " + loadpath + "pairitel_photo.sex " + image_name + \
        " -WEIGHT_IMAGE " + weight_image_name + \
        " -PARAMETERS_NAME " + loadpath + "star_find.param " + \
        "-CATALOG_NAME " + storepath + "star_cat.txt " + \
        "-DETECT_MINAREA 8 " + \
        "-BACK_SIZE 64 " + \
        "-BACK_FILTERSIZE 3 " + \
        "-BACK_TYPE AUTO " + \
        "-BACK_VALUE 0.0 " + \
        "-FILTER N " + \
        "-CHECKIMAGE_TYPE OBJECTS " + \
        "-CHECKIMAGE_NAME " + storepath + "check.fits " + \
        "-PHOT_APERTURES 5.8"
    
    system(syscmd)


    objects_image = storepath + "check.fits"
    hdulist = pyfits.open(objects_image)
    objects_data = hdulist[0].data
    hdulist.close()
    #system("rm "+ storepath + "check.fits")

    sex_file = file(storepath + "star_cat.txt", "r")
    sat_locations = []
    radec_locations = []
    for line in sex_file:
        if (line[0] == "#"):
            continue
        else:
            if float(line.split()[2]) == 0:
                x_index = int(round(float(line.split()[0]) - 1))
                y_index = int(round(float(line.split()[1]) - 1))
                RA_index = float(line.split()[3])
                Dec_index = float(line.split()[4])
                sat_locations.append([y_index, x_index])
                radec_locations.append([RA_index, Dec_index])
               # print str(RA_index) + ' ' + str(Dec_index)
    # print 'sat_locations before culling:'
    # print sat_locations

    
    if calreg:
        if find_fwhm or not ap:
            
            calreg_radec_list = openCalRegion(calreg)

            keep_it = False
            indexlist_fwhm = []

           # print 'calreg_radec_list is:'
           # print calreg_radec_list

            for index, old_radec in enumerate(radec_locations):
                keep_it = False
                #print 'old_radec is'
                #print old_radec
                for new_radec in calreg_radec_list:
                    calstar_radec_list = new_radec.split(',')
                    #print 'new_radec is'
                    #print calstar_radec_list
                    caldist = 206264.806247*(float(ephem.separation(((numpy.pi)*\
                            (1/180.)*(float(old_radec[0])), (numpy.pi)*(1/180.)*\
                            (float(old_radec[1]))),((numpy.pi)*(1/180.)*\
                            (float(calstar_radec_list[0])), (numpy.pi)*(1/180.)*\
                            (float(calstar_radec_list[1]))))))
                        # only remove the calstar entry if it is not our target;
                        # i.e. it is not identified as 'target' above
                        # Remove the calstar entry if distance is small
                   # print caldist
                    if caldist < 5:
                        keep_it = True
                if not keep_it:
                    indexlist_fwhm += [index]
                   # print 'removed index is'
                   # print index

            indexlist_fwhm.sort()
           # print 'here is the indexlist_fwhm of calibration stars which are not used in calculating FWHM:'
           # print indexlist_fwhm
           # print 'length of indexlist_fwhm is'
           # print len(indexlist_fwhm)
            revindex = indexlist_fwhm[::-1]
            for ind in revindex:
                # print "REMOVING %s" % str(sat_locations[ind])
                del sat_locations[ind]

           # print 'sat_locations after culling:'
           # print sat_locations
    
    if find_fwhm or not ap:
        if not ap:
            print 'Aperture is not specified; fitting FWHM to determine'
        fwhm, fwhm_stdev = fit_fwhm(sat_locations, objects_data, 2.8, 0.8)
        print "Initial fwhm fit complete."
        print ("FWHM mean: " + str(fwhm) + " stdev: " + str(fwhm_stdev) + ".")
        fwhm, fwhm_stdev = fit_fwhm(sat_locations, objects_data, fwhm, fwhm_stdev)
        print "Sigma clipping of fwhm complete."
        print ("FWHM mean: " + str(fwhm) + " stdev: " + str(fwhm_stdev) + ".")
        photdict.update({'FWHM':(fwhm,fwhm_stdev)})        
    
    if not ap:
        # If aperture is not specified, fit for ideal size based on the fwhm.
        aperture_size = (1.5)*fwhm
        print "Aperture size: ", aperture_size
        if isnan(aperture_size):
            print "Aperture size not found, using default 4.5."
            aperture_size = 4.5
    else:
        aperture_size = ap
    photdict.update({'Aperture':aperture_size})
    print 'Using aperture size: %f' % (aperture_size)

    '''Read in the image and convert data to array. Also store STRT_CPU, STOP_CPU, 
    and FILTER header keyword values for future use. If the start and stop times
    are not in the header, use placeholder values.'''
    imagefile_hdulist = pyfits.open(progenitor_image_name)
    imagefile_header = imagefile_hdulist[0].header
    image_data = imagefile_hdulist[0].data
    try:
        strt_cpu = str(imagefile_header["STRT_CPU"])
        stop_cpu = str(imagefile_header["STOP_CPU"])
        exptime = str(imagefile_header["EXPTIME"])
    except:
        raise Exception()
        strt_cpu = 999
        stop_cpu = 999
        exptime = 1
        
    band = str(imagefile_header["FILTER"])
    imagefile_hdulist.close()
    
    photdict.update({'STRT_CPU':strt_cpu})
    photdict.update({'STOP_CPU':stop_cpu})
    photdict.update({'EXPTIME':exptime})
    photdict.update({'filter':band})
    
    # Obtain the number of dither positions
    ditherdict = find_number_of_dither_positions(progenitor_image_name,fallback_n_dither=fallback_n_dither)
        
    if do_upper:
        # Store the image's dimensions.
        height, width = image_data.shape[0], image_data.shape[1]
        '''We will drop apertures of 3 pixel radius everywhere on the image where there
        exists science data. Since apertures are circular there will be gaps between
        adjacent apertures. Science data is determined as any pixel with a weightmap
        value > 0. To make sure we drop all possible apertures we need to create an 
        expanded weight mask. This is a mask array with value 1 where aperatures can 
        be centered and 0 elsewhere. It is formed by shrinking slightly the science
        data area from the original image (by about 4 pixels).'''
        imagefile_hdulist = pyfits.open(weight_image_name)
        imagefile_header = imagefile_hdulist[0].header
        weight_data = imagefile_hdulist[0].data
        imagefile_hdulist.close()
        weight_mask = where(weight_data > 0, 0, 1)
        weight_mask_expanded = zeros([height,width])
        binary_dilation(weight_mask, iterations = int(aperture_size), output = weight_mask_expanded)
        weight_mask_expanded = where(weight_mask_expanded == 0, 1, 0)
        # SourceExtractor likes to have a non-zero background, so fill it in with 
        # random data with a mean of 0.
        random_background_data = (-0.5 + 
            random.random_sample(width*height).reshape(height, width))
        '''To run SourceExtractor in two-image mode we need to first generate the fake
        stars image which dictates where the apertures will be dropped on the science
        image. We make fake stars at evenly spaced (6 pixels) intervals, blur them a 
        bit to make them look like stars, and then use the expanded weight mask to 
        enforce the boundaries.'''
        peak = 1000
        aperture_data = array([[0,0,0,0,0,0,0,0]])
        flux_array = array([])
        weight_aperture_data = array([[0,0,0,0,0,0,0,0]])
        order = 3
        spacing = aperture_size*order
        for g in range(order):
            for h in range(order):
                fake_star_data = zeros([height,width])
                for i in range(int(width/aperture_size/order)-1):
                    for j in range(int(height/aperture_size/order)-1):
                        fake_star_data[int((j+1)*spacing+aperture_size*h)][int((i+1)*spacing+aperture_size*g)] = peak
                blurred_fake_star_data = zeros([height,width])
                gaussian_filter(fake_star_data*weight_mask_expanded, 0.7, 
                    output=blurred_fake_star_data)
                output_hdu = pyfits.PrimaryHDU(blurred_fake_star_data + random_background_data)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto(progenitor_image_name.replace("coadd", "fakestars1_" + str(h) + "_" + str(g)))
                # Call Source Extractor in two-image mode for both sky apertures and weight
                # apertures.
                make_sex_cat_fake_stars(progenitor_image_name.replace("coadd", "fakestars1_" + str(h) + "_" + str(g)), 
                    progenitor_image_name, weight_image_name, aperture_size)
                make_sex_cat_fake_stars_no_weight(
                    progenitor_image_name.replace("coadd", "fakestars1_" + str(h) + "_" + str(g)), weight_image_name,
                    aperture_size)
                # Delete the fakestars image, we don't need it anymore.
                system("rm " + progenitor_image_name.replace("coadd", "fakestars1_" + str(h) + "_" + str(g)))
                # Read the SEX catalogs into data lists. The format of the sex_file is: 
                # ra, dec, inst_mag, e_inst_mag, flags, fwhm, signal, noise
    
                temp_array = loadtxt(file(progenitor_image_name.replace(".fits", ".sex")))
                aperture_data = append(aperture_data, temp_array)
                aperture_data = split(aperture_data, len(aperture_data)/8)

                flux_array = append(flux_array, (temp_array[:,6]))

                temp_array = loadtxt(file(weight_image_name.replace(".fits", ".sex")))
                weight_aperture_data = append(weight_aperture_data, temp_array)
                weight_aperture_data = split(weight_aperture_data, len(weight_aperture_data)/8)

        print len(aperture_data), len(weight_aperture_data)

        for g in range(order):
            for h in range(order):
                fake_star_data = zeros([height,width])
                for i in range(int(width/aperture_size/order)-1):
                    for j in range(int(height/aperture_size/order)-1):
                        fake_star_data[int((j+1-1./(order*2))*spacing+aperture_size*h)]\
                            [int((i+1-1./(order*2))*spacing+aperture_size*g)] = peak
                blurred_fake_star_data = zeros([height,width])
                gaussian_filter(fake_star_data*weight_mask_expanded, 0.7, 
                    output=blurred_fake_star_data)
                output_hdu = pyfits.PrimaryHDU(blurred_fake_star_data + random_background_data)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto(progenitor_image_name.replace("coadd", "fakestars2_" + str(h) + "_" + str(g)))
                # Call Source Extractor in two-image mode for both sky apertures and weight
                # apertures.
                make_sex_cat_fake_stars(progenitor_image_name.replace("coadd", "fakestars2_" + str(h) + "_" + str(g)), 
                    progenitor_image_name, weight_image_name, aperture_size)
                make_sex_cat_fake_stars_no_weight(
                    progenitor_image_name.replace("coadd", "fakestars2_" + str(h) + "_" + str(g)), weight_image_name,
                    aperture_size)
                # Delete the fakestars image, we don't need it anymore.
                system("rm " + progenitor_image_name.replace("coadd", "fakestars2_" + str(h) + "_" + str(g)))
                # Read the SEX catalogs into data lists. The format of the sex_file is: 
                # ra, dec, inst_mag, e_inst_mag, flags, fwhm, signal, noise
    
                temp_array = loadtxt(file(progenitor_image_name.replace(".fits", ".sex")))
                aperture_data = append(aperture_data, temp_array)
                aperture_data = split(aperture_data, len(aperture_data)/8)
        
                flux_array = append(flux_array, (temp_array[:,6]))
        
                temp_array = loadtxt(file(weight_image_name.replace(".fits", ".sex")))
                weight_aperture_data = append(weight_aperture_data, temp_array)
                weight_aperture_data = split(weight_aperture_data, len(weight_aperture_data)/8)
    
        print len(aperture_data), len(weight_aperture_data)
    
        dtype = [("ra", float), ("dec", float), ("mag", float), ("e_mag", float), 
            ("flags", float), ("fwhm", float), ("signal", float), ("noise", float)]
        list_o_tuples = []
        for entry in aperture_data[1:]:
            list_o_tuples.append((entry[0], entry[1], entry[2], entry[3], entry[4], 
                entry[5], entry[6], entry[7]))
        aperture_data = list_o_tuples
        list_o_tuples = []
        for entry in weight_aperture_data[1:]:
            list_o_tuples.append((entry[0], entry[1], entry[2], entry[3], entry[4], 
                entry[5], entry[6], entry[7]))
        weight_aperture_data = list_o_tuples

        aperture_data = array(aperture_data, dtype=dtype)
        weight_aperture_data = array(weight_aperture_data, dtype=dtype)

        aperture_data = sort(aperture_data, order=["ra", "dec"])
        weight_aperture_data = sort(weight_aperture_data, order=["ra", "dec"])

        # savetxt("aperture_data.txt", aperture_data)
        # savetxt("weight_aperture_data.txt", weight_aperture_data)
        
        if len(aperture_data) != len (weight_aperture_data):
            print "Aperture Data list length mismatch:"
            print len(aperture_data), len(weight_aperture_data)
            raise(ValueError)
    
        combined_aperture_data = []
        flux_list = []
        mismatch_num = 0
    
        for n in range(len(aperture_data)):
            sky_ra = aperture_data[n][0]
            sky_dec = aperture_data[n][1]
            sky_signal = aperture_data[n][6]
            weight_ra = weight_aperture_data[n][0]
            weight_dec = weight_aperture_data[n][1]
            weight_signal = weight_aperture_data[n][6]
            if str(sky_ra)[:7] == str(weight_ra)[:7] and str(sky_dec)[:7] == str(weight_dec)[:7]:
                combined_aperture_data.append([sky_ra, sky_dec, sky_signal, weight_signal])
                flux_list.append(sky_signal)
            else:
                mismatch_num += 1
        #         print (n, "WARNING: sky and weight aperture miss-match at sky location:" + 
        #             str(sky_ra) + ", " + str(sky_dec) + " and weight location:" + 
        #             str(weight_ra) + ", " + str(weight_dec))
        print "Sky and Weight apperture matching resulted in ", mismatch_num, " fails."
        photdict.update({'SkyWeightMatchFails':mismatch_num})
    
        # Convert the flux_list to flux_array and then use median absolute deviation
        # clipping to generate a new median and sigma.
        flux_array = array(flux_list)
        c_array, new_medval, new_sigma = mad_clipping(flux_array, 2)
        # Filter the entries in aperture_data with the new_medval and new_sigma. Store
        # the surviving data in aperture_data_2 and the flux values in flux_list_2.
        aperture_data_2 = []
        flux_list_2 = []
        weight_flux_list = []
        for aperture in combined_aperture_data:
            if ((aperture[2] > new_medval - 3.0*new_sigma) and 
                (aperture[2] < new_medval + 1.0*new_sigma)):
                    aperture_data_2.append(aperture)
                    flux_list_2.append(aperture[2])
                    weight_flux_list.append(aperture[3])
            

        # Replace flux_array with the new values from flux_list_2. Also, store the
        # number of apertures.
        flux_array = array(flux_list_2)
        num_apertures = len(flux_array)
        # Divide the weight flux list into quartiles. This is the basis of the 
        # color/quality coding.
        median_weight_flux = median(weight_flux_list)
        high_half = []
        low_half = []
        for val in weight_flux_list:
            if val >= median_weight_flux:
                high_half.append(val)
            else:
                low_half.append(val)
        high_median_weight_flux = median(high_half)
        low_median_weight_flux = median(low_half)
        # Use the entries in aperture_data_2 to create a region file with all the 
        # apertures. Very useful for diagnostics.
        flux_list_green = []
        flux_list_yellow = []
        flux_list_orange = []
        flux_list_red = []
        weight_flux_list_green = []
        weight_flux_list_yellow = []
        weight_flux_list_orange = []
        weight_flux_list_red = []
        output_region_file = file(storepath + "aperature_regions.reg", "w")
        output_region_file.write("global color=green dashlist=8 3 width=1 " + 
            "font='helvetica 10 normal' select=1 highlite=1 dash=0 fixed=0 edit=1 " + 
            "move=1 delete=1 include=1 source=1\nfk5\n")
        for aperture in aperture_data_2:
            if aperture[3] >= high_median_weight_flux:
                output_region_file.write('circle(' + str(aperture[0]) + ',' + 
                    str(aperture[1]) + ',' + str(aperture_size/2) + '") # color=green\n')
                flux_list_green.append(aperture[2])
                weight_flux_list_green.append(aperture[3])
            if (aperture[3] < high_median_weight_flux) and (aperture[3] >= median_weight_flux):
                output_region_file.write('circle(' + str(aperture[0]) + ',' + 
                    str(aperture[1]) + ',' + str(aperture_size/2) + '") # color=yellow\n')
                flux_list_yellow.append(aperture[2])
                weight_flux_list_yellow.append(aperture[3])
            if (aperture[3] < median_weight_flux) and (aperture[3] >= low_median_weight_flux):
                output_region_file.write('circle(' + str(aperture[0]) + ',' + 
                    str(aperture[1]) + ',' + str(aperture_size/2) + '") # color=orange\n')
                flux_list_orange.append(aperture[2])
                weight_flux_list_orange.append(aperture[3])
            if aperture[3] < low_median_weight_flux:
                output_region_file.write('circle(' + str(aperture[0]) + ',' + 
                    str(aperture[1]) + ',' + str(aperture_size/2) + '") # color=red\n')
                flux_list_red.append(aperture[2])
                weight_flux_list_red.append(aperture[3])
        output_region_file.close()
        flux_array_green = array(flux_list_green)
        flux_array_yellow = array(flux_list_yellow)
        yellow = False
        if len(flux_array_yellow) > 50:
            yellow = True
        flux_array_orange = array(flux_list_orange)
        orange = False
        if len(flux_array_orange) > 50:
            orange = True
        flux_array_red = array(flux_list_red)
        red = False
        if len(flux_array_red) > 50:
            red = True
        '''Create a histogram from the flux_array. Convert the flux_histogram_data into a
        list for use with the gaussian fitting. Convert the flux_histogram_bins into a
        list of length one less with values replaced with the average value of each 
        bin. This is necessary because the output bins array is actually the bin
        boundaries. Repeat this code block 4 times, once for each color/quality subset.'''
        # GREEN
        flux_histogram_data_green, flux_histogram_bins_green = histogram(
            flux_array_green, int(num_apertures/160))
        bin_half_height_green = (flux_histogram_bins_green[1] - 
            flux_histogram_bins_green[0])
        flux_histogram_data_list_green = list(flux_histogram_data_green)
        flux_histogram_bins_list_green = (list(bin_half_height_green + 
            flux_histogram_bins_green[0:-1]))
        # Fit a gaussian curve to the histogram of the flux data.
        central_height_green, center_x_green, x_width_green = (fit_gaussian(
            flux_histogram_data_list_green, flux_histogram_bins_list_green))
        # In case the fit returns a negative sigma, take the absolute value.
        x_width_green = abs(x_width_green)
        # Create data from the gaussian fit to plot on the flux histogram.
        fit_data_list_green = []
        for bin in flux_histogram_bins_list_green:
            fit_point = ((central_height_green/sqrt(2*pi*x_width_green**2)) * 
                exp((-1*(bin-center_x_green)**2)/(2*x_width_green**2)))
            fit_data_list_green.append(fit_point)
        # YELLOW
        if yellow:
            flux_histogram_data_yellow, flux_histogram_bins_yellow = histogram(
                flux_array_yellow, int(num_apertures/160))
            bin_half_height_yellow = (flux_histogram_bins_yellow[1] - 
                flux_histogram_bins_yellow[0])
            flux_histogram_data_list_yellow = list(flux_histogram_data_yellow)
            flux_histogram_bins_list_yellow = (list(bin_half_height_yellow + 
                flux_histogram_bins_yellow[0:-1]))
            # Fit a gaussian curve to the histogram of the flux data.
            central_height_yellow, center_x_yellow, x_width_yellow = (fit_gaussian(
                flux_histogram_data_list_yellow, flux_histogram_bins_list_yellow))
            # In case the fit returns a negative sigma, take the absolute value.
            x_width_yellow = abs(x_width_yellow)
            # Create data from the gaussian fit to plot on the flux histogram.
            fit_data_list_yellow = []
            for bin in flux_histogram_bins_list_yellow:
                fit_point = ((central_height_yellow/sqrt(2*pi*x_width_yellow**2)) * 
                    exp((-1*(bin-center_x_yellow)**2)/(2*x_width_yellow**2)))
                fit_data_list_yellow.append(fit_point)
        # ORANGE
        if orange:
            flux_histogram_data_orange, flux_histogram_bins_orange = histogram(
                flux_array_orange, int(num_apertures/160))
            bin_half_height_orange = (flux_histogram_bins_orange[1] - 
                flux_histogram_bins_orange[0])
            flux_histogram_data_list_orange = list(flux_histogram_data_orange)
            flux_histogram_bins_list_orange = (list(bin_half_height_orange + 
                flux_histogram_bins_orange[0:-1]))
            # Fit a gaussian curve to the histogram of the flux data.
            central_height_orange, center_x_orange, x_width_orange = (fit_gaussian(
                flux_histogram_data_list_orange, flux_histogram_bins_list_orange))
            # In case the fit returns a negative sigma, take the absolute value.
            x_width_orange = abs(x_width_orange)
            # Create data from the gaussian fit to plot on the flux histogram.
            fit_data_list_orange = []
            for bin in flux_histogram_bins_list_orange:
                fit_point = ((central_height_orange/sqrt(2*pi*x_width_orange**2)) * 
                    exp((-1*(bin-center_x_orange)**2)/(2*x_width_orange**2)))
                fit_data_list_orange.append(fit_point)
        # RED
        if red:
            flux_histogram_data_red, flux_histogram_bins_red = histogram(flux_array_red, 
                int(num_apertures/160))
            bin_half_height_red = flux_histogram_bins_red[1] - flux_histogram_bins_red[0]
            flux_histogram_data_list_red = list(flux_histogram_data_red)
            flux_histogram_bins_list_red = list(bin_half_height_red + 
                flux_histogram_bins_red[0:-1])
            # Fit a gaussian curve to the histogram of the flux data.
            central_height_red, center_x_red, x_width_red = (fit_gaussian(
                flux_histogram_data_list_red, flux_histogram_bins_list_red))
            # In case the fit returns a negative sigma, take the absolute value.
            x_width_red = abs(x_width_red)
            # Create data from the gaussian fit to plot on the flux histogram.
            fit_data_list_red = []
            for bin in flux_histogram_bins_list_red:
                fit_point = ((central_height_red/sqrt(2*pi*x_width_red**2)) * 
                    exp((-1*(bin-center_x_red)**2)/(2*x_width_red**2)))
                fit_data_list_red.append(fit_point)
        # Generate the flux histogram plot.
        plot_title = ("Sky Flux Histogram")
        fig = plt.figure(figsize=(6, 6))
        ax1 = fig.add_subplot(1,1,1)
        ax1.plot(flux_histogram_bins_list_green, flux_histogram_data_list_green, 
            marker = "o", color = "green", linestyle="none", label = "Flux [ADU]")
        ax1.plot(flux_histogram_bins_list_green, fit_data_list_green, marker = ".", 
            color = "green", linestyle="solid", label = "Gaussian Fit")
        if yellow:
            ax1.plot(flux_histogram_bins_list_yellow, flux_histogram_data_list_yellow, 
                marker = "o", color = "yellow", linestyle="none", label = "Flux [ADU]")
            ax1.plot(flux_histogram_bins_list_yellow, fit_data_list_yellow, marker = ".", 
                color = "yellow", linestyle="solid", label = "Gaussian Fit")
        if orange:
            ax1.plot(flux_histogram_bins_list_orange, flux_histogram_data_list_orange, 
                marker = "o", color = "orange", linestyle="none", label = "Flux [ADU]")
            ax1.plot(flux_histogram_bins_list_orange, fit_data_list_orange, marker = ".", 
                color = "orange", linestyle="solid", label = "Gaussian Fit")
        if red:
            ax1.plot(flux_histogram_bins_list_red, flux_histogram_data_list_red, 
                marker = "o", color = "red", linestyle="none", label = "Flux [ADU]")
            ax1.plot(flux_histogram_bins_list_red, fit_data_list_red, marker = ".", 
                color = "red", linestyle="solid", label = "Gaussian Fit")
        ax1.set_ylabel("Number")
        ax1.set_xlabel("<- Fainter     Bin     Brighter ->")
        ax1.set_title(plot_title)
        canvas = FigureCanvas(fig)
        canvas.print_figure("flux_histogram.png", dpi=144)
    # Now we do photometry. We use the region file to extract the target ra and dec.
    region_file = file(region_file, "r")
    for line in region_file:
        if line[:6] == "circle":
            target_ra_deg = float(line.split("(")[1].split(",")[0])
            target_dec_deg = float(line.split("(")[1].split(",")[1])
    target_ra_rad = target_ra_deg * 0.01745329252
    target_dec_rad = target_dec_deg * 0.01745329252

    # We collect the 2MASS photometry of field sources with a query to vizier.
    # Most of this is forming the query and then formatting the returned text file.
    if not reuse_old_vizcat:
        viz_input_file = file(storepath + "viz_input.txt", "w")
        viz_input_file.write(
            "-source=2MASS/PSC\n" + 
            "-c=" + str(target_ra_deg) + " " + str(target_dec_deg) + "\n" +
            "-c.rm=15\n" +
            "-out=_r RAJ2000 DEJ2000 " + 
                "Jmag e_Jmag Jsnr " + 
                "Hmag e_Hmag Hsnr " + 
                "Kmag e_Kmag Ksnr " + 
                "Qflg\n" +
            "-sort=_r\n" + 
            "Qflg==AAA\n" + 
            "Jsnr=>10"
            + "\nJmag=>9"
            )
        viz_input_file.close()
        syscmd = "vizquery -mime=csv -site=cfa %sviz_input.txt > %sviz_output.txt" % (storepath, storepath)
        system(syscmd)
        system("rm -f viz_input.txt")
        viz_output_cropped_file = file(storepath + "viz_output_cropped.txt", "w")
        viz_output_file = file(storepath + "viz_output.txt", "r")
        line_num = 0
        for line in viz_output_file:
            line_num += 1
            if line_num > 46 and line != "\n":
                viz_output_cropped_file.write(line)
        viz_output_cropped_file.close()
        viz_output_file.close()
    else:
        print "Reusing old vizcat in the store directory!"
    # Define the obj_string and band_type from the progenitor_image_name. If the 
    # progenitor_image_name is not in the format used by PARITIEL reduction
    # pipeline 3, these will not work properly.
    try:
        obj_string = progenitor_image_name.split("_")[2]
        band_type = (progenitor_image_name.split("_")[0] + "_" + 
            progenitor_image_name.split("_")[1])
    except:
        raise ValueError('Image name misformatted; needs to be in standard Pipeline3 file format')
    # Run source extractor on the science image.
    sexcat = make_sex_cat(progenitor_image_name, weight_image_name, aperture_size)
    # Store placehold values for the target photometry.
    target_mag = 999
    target_e_mag = 999
    # Create the catalog of 2MASS stars from the vizier data which corrresponds to
    # the filter/band used for the science image.
    if not stardict:
        vizcat_file = file(storepath + "viz_output_cropped.txt", "r")
        vizcat_starlist = []
        for line in vizcat_file:
            data_list = line.rstrip().lstrip().split(";")
            if len(data_list) == 13: # Check formatting 
                try:
                    ra = float(data_list[1]) # degrees
                    dec = float(data_list[2]) # degrees
                    if band == "j":
                        mag = float(data_list[3])
                        e_mag = float(data_list[4])
                        snr = float(data_list[5])
                    if band == "h":
                        mag = float(data_list[6])
                        e_mag = float(data_list[7])
                        snr = float(data_list[8])
                    if band == "k":
                        mag = float(data_list[9])
                        e_mag = float(data_list[10])
                        snr = float(data_list[11])
                    vizcat_starlist.append([ra, dec, mag, e_mag, snr])
                except:
                    print "Problem parsing line in viz_output.txt:"
                    print line
        vizcat_file.close()
        
    # Using deepstack of calibration stars as the catalog, instead of 2MASS
    else:
        vizcat_starlist = []
        caldict = stardict[band] 
        for star in caldict:
            if 'targ_mag' not in caldict[star]:
                pass
            else:
                star_ra = float(star.strip('(').strip(')').split(',')[0].lstrip("'").rstrip("'"))
                star_dec = float(star.strip('(').strip(')').split(',')[1].lstrip(" ").lstrip("'").rstrip("'"))
                star_mag = caldict[star]['targ_mag'][0]
                star_e_mag = caldict[star]['targ_mag'][1]
                star_snr = caldict[star]['targ_s2n']
                vizcat_starlist.append([star_ra, star_dec, star_mag, star_e_mag, star_snr])
    
    if vizcat_starlist:
        print 'first entry in vizquery starlist [ra, dec, mag, e_mag, snr]:'
        print vizcat_starlist[0]
    else:
        print '\n\n\n****No Vizquery starlist!!****\n\n\n'
    
    # Create the sexcat_starlist from the Source Extractor output catalog. Also fill
    # in the sex_inst_mag_list which will be used as a diagnostic check on the 
    # computed upper limit.
    sexcat_starlist = []
    sex_inst_mag_list = []
    sex_inst_flux_list = []
    sexcat_file = file(sexcat, "r")
    for line in sexcat_file:
        ra = float(line.split()[0])
        dec = float(line.split()[1])
        mag = float(line.split()[2])
        mag_err = float(line.split()[3])
        flags = int(line.split()[4])
        fwhm_image = float(line.split()[5])
        flux = float(line.split()[6])
        flux_err = float(line.split()[7])
        sexcat_starlist.append([ra, dec, mag, mag_err, flags,fwhm_image,flux,flux_err])
        if flags == 0:
            sex_inst_mag_list.append([mag, mag_err, ra, dec])
            sex_inst_flux_list.append([flux,flux_err, ra, dec])
    sexcat_file.close()
    
    print "Length of pre-truncated sextractor starlist is %s" % str(len(sexcat_starlist))
    print "Length of vizquery 2mass starlist is %s" % str(len(vizcat_starlist))
    # Compare the entries in sexcat_starlist and vizcat_starlist to create a 
    # combined_starlist which has entries for sources with both archival and new
    # instrumental data. The target need not be included in the archive.
    # [ra, dec, 2massmag, 2masserr instmag, insterr, sexflags, 'calib'/'target']

    # debug printing
    # import pprint
    # print 'sexcat_starlist is'
    # pprint.pprint(sexcat_starlist)

    # print '================================================================='
    # print '================================================================='
    # print '================================================================='
    # print 'vizcat_starlist is'
    # pprint.pprint(vizcat_starlist)
    
    combined_starlist = []
    for sexcat_star in sexcat_starlist:
        sexcat_star_ra_rad = sexcat_star[0] * 0.01745329252
        sexcat_star_dec_rad = sexcat_star[1] * 0.01745329252
        # First check to see if our target location was detected by sextractor
        # If it was, add it to the starlist.  
        target_separation_arcsec = 206264.806247*(float(ephem.separation(
            (target_ra_rad, target_dec_rad), 
            (sexcat_star_ra_rad, sexcat_star_dec_rad))))
        if target_separation_arcsec < 5:
            combined_starlist.append([sexcat_star[0], sexcat_star[1], 
                999, 999, 
                sexcat_star[2], sexcat_star[3], 
                sexcat_star[4],'target'])
            src_flux = sexcat_star[6]
            src_flux_err = sexcat_star[7]
            src_s2n = src_flux/src_flux_err
            print "SOURCE FOUND WITH S/N OF %f (%f +/- %f)" % (src_s2n, src_flux, src_flux_err)
            continue
        # Now, cross correlate the calibration catalog with our instrumental
        # Magnitudes, and store any matches in the starlist
        for calibcat_star in vizcat_starlist:
            calibcat_star_ra_rad = calibcat_star[0] * 0.01745329252
            calibcat_star_dec_rad = calibcat_star[1] * 0.01745329252
            separation_arcsec = 206264.806247*(float(ephem.separation(
                (calibcat_star_ra_rad, calibcat_star_dec_rad), 
                (sexcat_star_ra_rad, sexcat_star_dec_rad))))
            if separation_arcsec < 5:
                combined_starlist.append([calibcat_star[0], calibcat_star[1], 
                    calibcat_star[2], calibcat_star[3], 
                    sexcat_star[2], sexcat_star[3], 
                    sexcat_star[4],'calib',separation_arcsec])

    
    print 'length of combined_starlist is %s' % str(len(combined_starlist))
    if len(combined_starlist) == 0:
        if not for_wcs_coadd:
            raise Exception('Combined Starlist is empty; no cross correlations found between sextractor catalog and 2mass.')
        else:
            photdict = {'FileName':image_name}
            return photdict
    # Deleting stars in combined_starlist that is not in the inputted region file (if not calreg == False)
    if not calreg:
        pass
    else:
        callist = openCalRegion(calreg)
        
        keep_it = False
        indexlist = []

        for index, calstar_old in enumerate(combined_starlist):
            keep_it = False
            for calstar_new in callist:
                calstar_list = calstar_new.split(',')
                caldist = 206264.806247*(float(ephem.separation(((numpy.pi)*\
                    (1/180.)*(float(calstar_old[0])), (numpy.pi)*(1/180.)*\
                    (float(calstar_old[1]))),((numpy.pi)*(1/180.)*\
                    (float(calstar_list[0])), (numpy.pi)*(1/180.)*\
                    (float(calstar_list[1]))))))
                # only remove the calstar entry if it is not our target;
                # i.e. it is not identified as 'target' above
                # Remove the calstar entry if distance is small
                if caldist < 5 or calstar_old[7] == 'target':
                    keep_it = True
            if not keep_it:
                indexlist += [index]

        indexlist.sort()
        if verbose:
            print 'here is the indexlist of calibration stars which are not used'
            print indexlist
            print 'length of indexlist is'
            print len(indexlist)
        revindex = indexlist[::-1]
        for ind in revindex:
            if verbose:
                print "REMOVING %s" % str(combined_starlist[ind])
            del combined_starlist[ind]
            
    # Use the combined_starlist to calculate a zeropoint for the science image.
    zeropoint_list = []
    zeropoint_err_list = []
    
    '''In PAIRITEL, there are additional sources of uncertainty which
    have to be accounted for above the normal photometric statistical
    uncertainty.  One is an additive error which reduces by the sqrt
    of the number of dither positions used, and the second is a 
    multiplicitive factor due to the fact that we are typically 
    rebinning from 2"/pix to 1"/pix (should be a factor of two, but
    in testing with GRB 071025, Dan found it to be ~2.4).'''
    
    num_triplestacks = ditherdict['N_dither_pos']
    if num_triplestacks == 0:
        num_triplestacks = 1    

    for star in combined_starlist:
        ra = star[0]
        dec = star[1]
        tmass_mag = star[2]
        tmass_e_mag = star[3]
        ptel_mag = star[4]
        ptel_e_mag = star[5]
        ptel_flag = star[6]
        star_ra_rad = ra * 0.01745329252
        star_dec_rad = dec * 0.01745329252
        target_separation_arcsec = 206264.806247*(float(ephem.separation(
                (target_ra_rad, target_dec_rad), 
                (star_ra_rad, star_dec_rad))))
        # Don't add it to the zeropoint if it is our target
        if ((target_separation_arcsec > 5) and 
            ((ptel_flag == 0) or (ptel_flag == 2))):
            zeropoint_list.append(tmass_mag - ptel_mag)
            if star[7] != 'calib': # sanity check that different target matching works
                raise ValueError('Calib Star identification mismatch somewhere in the code; squash this bug')
            sub_zp_err = return_ptel_uncertainty(ptel_e_mag,tmass_e_mag,num_triplestacks)

            if caliblimit:
                calregpath = calreg
                calreglist = openCalRegion(calregpath)
                if sub_zp_err > 0.3:
                    errmessage = '''Star at RA: %f Dec: %f has a huge uncertainty.  
                    Investigate and possibly remove from calibration star list before continuing''' % (ra, dec)                  
                    # automatic culling of calibration stars
                    if autocull:
                        print errmessage
                        for index, calreg_stars in enumerate(calreglist):
                            calreg_ra = float(calreg_stars.split(',')[0])
                            calreg_dec = float(calreg_stars.split(',')[1])
                            if (abs(ra-calreg_ra) < 0.00001) & (abs(dec-calreg_dec) < 0.00001):
                                 print 'calreglist first'
                                 print calreglist
                                 del calreglist[index]
                                 print 'calibration star with calreg values RA: %f, Dec : %f has been removed' % (calreg_ra, calreg_dec)
                                 print 'calreglist second'
                                 print calreglist
                                 print 'calreglist second length is'
                                 print len(calreglist)
                    else:
                        raise Exception(errmessage)          
            
            # zeropoint_err_list.append(sqrt(tmass_e_mag**2 + (ptel_e_mag*2.4)**2)
            #             + (base_dither_error/sqrt(num_triplestacks))**2)
            zeropoint_err_list.append(sub_zp_err)
            
            #Sanity check on target separation
            if star[8] > 5:
                raise Exception("Investigate why the calib star/2mass separation is so large")
            
    # zeropoint = average(zeropoint_list) # average is not very robust
    zeropoint = numpy.median(zeropoint_list)
    if numpy.isnan(zeropoint):
        print 'ZEROPOINT IS NAN - something is wrong.  Here is zeropoint_list:'
        print zeropoint_list
        #raise(ValueError)
    
    # AMorgan edits to a different calculation of Zeropoint Uncertainty
    N_zp = len(zeropoint_list)
    zp_err_arr = numpy.array(zeropoint_err_list)
    zp_err_arr_sq = zp_err_arr**2
    sum_zp_err_arr_sq = numpy.sum(zp_err_arr_sq)
#    zeropoint_error = average(zeropoint_err_list)
    zeropoint_error = numpy.sqrt(sum_zp_err_arr_sq)/N_zp


    # Now apply the zeropoint to the instrumental magnitudes and create the 
    # final_starlist. Store the target photometry in target_mag and target_e_mag.
    final_starlist = []
    photdict.update({'calib_stars':{}})
    for star in combined_starlist:
        # If the star is our target . . .
        if (206264.806247*(float(ephem.separation(
            (target_ra_rad, target_dec_rad), 
            (star[0] * 0.01745329252, star[1] * 0.01745329252)))) < 5 and
            ((star[6] == 0) or (star[6] == 2) or (star[6] == 4) or (star[6] == 6) or 
            (star[6] == 7))):
            if star[7] != 'target': # sanity check that different target matching works
                raise ValueError('Target identification mismatch somewhere in the code; squash this bug')
            ra = star[0]
            dec = star[1] 
            tmass_mag = star[2]
            tmass_e_mag = star[3]
            ptel_mag = star[4]
            ptel_e_mag = star[5]
            ptel_flag = star[6]
            new_mag = ptel_mag + zeropoint
      
            '''In PAIRITEL, there are additional sources of uncertainty which
            have to be accounted for above the normal photometric statistical
            uncertainty.  One is an additive error which reduces by the sqrt
            of the number of dither positions used, and the second is a 
            multiplicitive factor due to the fact that we are typically 
            rebinning from 2"/pix to 1"/pix (should be a factor of two, but
            in testing with GRB 071025, Dan found it to be ~2.4).'''
            
            # preserving the target magnitude for later
            target_mag = new_mag
            
            # e_mag = sqrt((inst_err*2.4)**2 + zp_err**2 + (base_dither_error/sqrt(num_triplestacks))**2)

            # new_e_mag = float(sqrt(zeropoint_error**2 + (ptel_e_mag*2.4)**2 
            #                 + (base_dither_error/sqrt(num_triplestacks))**2))
            new_e_mag = return_ptel_uncertainty(ptel_e_mag,zeropoint_error,num_triplestacks)
            target_e_mag = new_e_mag
            target_instrumental_e_mag = ptel_e_mag

            target_flux = src_flux # Note does not take into account zp
            target_flux_err = src_flux_err # Note does not take into account zp
            final_starlist.append([ra, dec, tmass_mag, tmass_e_mag, 
                ptel_mag, ptel_e_mag, ptel_flag, new_mag, new_e_mag])
            if verbose:
                print ''
                print '******'
                print '!!!!YAY our target!!!!'
                print star
                print ''
            continue
        # If the star is just a field 2MASS star . . .
        elif star[6] == 0:
            ra = star[0]
            dec = star[1]
            tmass_mag = star[2]
            tmass_e_mag = star[3]
            ptel_mag = star[4]
            ptel_e_mag = star[5]
            ptel_flag = star[6]
            new_mag = ptel_mag + zeropoint
            delta_pos = star[8]
            delta_mag = new_mag-tmass_mag # difference between our magnitude and 2mass
            # new_e_mag = float(sqrt(zeropoint_error**2 + (ptel_e_mag*2.4)**2 
            #                 + (base_dither_error/sqrt(num_triplestacks))**2))
            new_e_mag = return_ptel_uncertainty(ptel_e_mag,zeropoint_error,num_triplestacks)
            final_starlist.append([ra, dec, tmass_mag, tmass_e_mag, 
                ptel_mag, ptel_e_mag, ptel_flag, new_mag, new_e_mag])
            photdict['calib_stars'].update({str(ra):{'ra':ra,'dec':dec,
                '2mass_mag':tmass_mag,'2mass_e_mag':tmass_e_mag,
                'inst_mag':ptel_mag,'inst_e_mag':ptel_e_mag,'new_mag':new_mag,
                'new_e_mag':new_e_mag, 'delta_mag':delta_mag, 'delta_pos':delta_pos}})
            if verbose:
                print ''
                print '******'
                print 'yay good star'
                print star
                print ''
        else:
            if verbose:
                print ''
                print '******'
                print 'BLARGH CANNOT USE THIS STAR'
                print star
                print ''
    # Calculate the midpoint heliocentric julian date of the exposure. We use a 
    # try/except clause in case something fails and use a placeholder hjd in that
    # instance.
    try:
        start_time = ephem.date(strt_cpu.replace("-", "/"))
        stop_time = ephem.date(stop_cpu.replace("-", "/"))
        hjd_start_time = heliocentric_julian_date(start_time, target_ra_rad, 
            target_dec_rad)
        hjd_stop_time = heliocentric_julian_date(stop_time, target_ra_rad, 
            target_dec_rad)
        hjd = (hjd_stop_time + hjd_start_time)/2
    except:
        hjd = 999
        hjd_start_time = 999
        hjd_stop_time = 999
    
    photdict.update({"HJD_mid":hjd})
    photdict.update({"HJD_start":hjd_start_time})
    photdict.update({"HJD_stop":hjd_stop_time})
    
    # Write out a final catalog of photometry data. This is useful for diagnostic
    # purposes to make sure that the newly observed magntidues are in line with the
    # archival 2MASS values. The target's archival data is replaced with the 
    # placeholder values (999).
    abs_2mass_deviation_list = []
    output_file = file(sexcat.replace("sex", "finalcat") + ".txt", "w")
    output_file.write(str(hjd) + "\nRA\t\t\tDEC\t\t\t\t2mass_mag\t\tnew_mag" + 
        "\t\t\tnew_e_mag\n")
    for star in final_starlist:
        output_file.write(("%f" + "\t" + "%f" + "\t\t" + 
            "%f" + "\t\t" + "%f" + "\t\t" + "%f" + "\n") % (star[0], 
            star[1], star[2], star[7], star[8]))
        if star[2] != 999:
            abs_2mass_deviation_list.append(abs(star[2]-star[7]))
    output_file.close()
    # Form the photometry string which will be printed out as the final result.
    photometry_string = "Source Extractor detects no source at target position."
    if str(target_mag) != "999":
        photometry_string = ("%s_mag: %.3f \terr: %f" % 
            (band_type, target_mag, target_e_mag))
        # Update the photometry string with updated values
        photdict.update({'targ_mag':(target_mag,target_e_mag)})
        photdict.update({'targ_flux':(target_flux,target_flux_err)})
        photdict.update({'targ_inst_err':target_instrumental_e_mag})
        photdict.update({'targ_s2n':src_s2n})
    # Print out photometry data.
    print progenitor_image_name, "HJD:", hjd
    print "Photometry Results:", photometry_string
    print ("2MASS Catalog comparison average absolute deviation: " + 
        str(average(abs_2mass_deviation_list)))
    print "Zeropoint:", zeropoint, "err", zeropoint_error
    
    photdict.update({'zp_list':zeropoint_list})
    photdict.update({'zp_err_list':zeropoint_err_list})
    photdict.update({'zp':(zeropoint,zeropoint_error)})
    photdict.update({'2mass_abs_avg_dev':str(average(abs_2mass_deviation_list))})
    
    
    # Sort the instrumental magnitudes to find the faintest detection.
    sex_inst_mag_list_sorted = sorted(sex_inst_mag_list, key=operator.itemgetter(1))
    faintest_mag = sex_inst_mag_list_sorted[-1][0]
    faintest_mag_err = sex_inst_mag_list_sorted[-1][1]
    faintest_ra = sex_inst_mag_list_sorted[-1][2]
    faintest_dec = sex_inst_mag_list_sorted[-1][3]
    sex_inst_flux_list_sorted = sorted(sex_inst_flux_list, key=operator.itemgetter(1))
    faintest_flux = sex_inst_flux_list_sorted[0][0]
    faintest_flux_err = sex_inst_flux_list_sorted[0][1]
    # should be redundant
    faintest_flux_ra = sex_inst_flux_list_sorted[0][2]
    faintest_flux_dec = sex_inst_flux_list_sorted[0][3]
    # but raise an exception just in case they're not
    # if faintest_flux_dec != faintest_dec or faintest_flux_ra != faintest_ra:
    #     raise Exception('The Faintest Flux != Faintest Mag.  WTF.')
        
    if do_upper:
        # Print the upper limit data.
        print "Number of sky aperatures: " + str(num_apertures)
        photdict.update({'Num_UL_apertures':num_apertures})
    
        greenpixavgweight = average(weight_flux_list_green)/(pi*aperture_size*aperture_size/4)
        green_ul_mag = -2.5*log(center_x_green + 3*x_width_green, 10) + zeropoint
        print ("Gaussian Fit Upper Limit green (avg weight pixel=" + 
            str(greenpixavgweight) + "): " + 
            str(green_ul_mag) + 
            " err " + str(zeropoint_error))
    
        photdict.update({'upper_green':green_ul_mag})
        photdict.update({'upper_green_avgpix':greenpixavgweight})
    
        if yellow:
            yellowpixavgweight = average(weight_flux_list_yellow)/(pi*aperture_size*aperture_size/4)
            yellow_ul_mag = -2.5*log(center_x_yellow + 3*x_width_yellow, 10) + zeropoint
            print ("Gaussian Fit Upper Limit yellow (avg weight pixel=" + 
                str(yellowpixavgweight) + "): " + 
                str(yellow_ul_mag) + 
                " err " + str(zeropoint_error))
            photdict.update({'upper_yellow':yellow_ul_mag})
            photdict.update({'upper_yellow_avgpix':yellowpixavgweight})
        
        if orange:
            orangepixavgweight = average(weight_flux_list_orange)/(pi*aperture_size*aperture_size/4)
            orange_ul_mag = -2.5*log(center_x_orange + 3*x_width_orange, 10) + zeropoint
            print ("Gaussian Fit Upper Limit orange (avg weight pixel=" + 
                str(orangepixavgweight) + "): " + 
                str(orange_ul_mag) + 
                " err " + str(zeropoint_error))
            photdict.update({'upper_orange':orange_ul_mag})
            photdict.update({'upper_orange_avgpix':orangepixavgweight})
        
        if red:
            redpixavgweight = average(weight_flux_list_red)/(pi*aperture_size*aperture_size/4)
            red_ul_mag = -2.5*log(center_x_red + 3*x_width_red, 10) + zeropoint
            print ("Gaussian Fit Upper Limit red (avg weight pixel=" + 
                str(redpixavgweight) + "): " + 
                str(red_ul_mag) + 
                " err " + str(zeropoint_error))
            photdict.update({'upper_red':red_ul_mag})
            photdict.update({'upper_red_avgpix':redpixavgweight})
    
    sex_faintest = faintest_mag + zeropoint
    sex_faintest_err = sqrt(faintest_mag_err**2 + zeropoint_error**2)
    print ("SExtractor faintest detection: " + str(sex_faintest) + 
        " err " + str(sex_faintest_err) + " at " + 
        str(faintest_ra) + ", " + str(faintest_dec))
#comments out s2n starts////////////////////////////
    faintest_s2n = faintest_flux/faintest_flux_err
    print ("SExtractor faintest flux: " + str(faintest_flux) + 
        " err " + str(faintest_flux_err) + " => S/N = " + str(faintest_s2n))
    
    photdict.update({'faintest_s2n':faintest_s2n})
#comments out s2n stops////////////////////////////////////////    
    photdict.update({'sex_faintest':(sex_faintest,sex_faintest_err)})
    photdict.update(ditherdict)

    # Clean up the photometry catalog files.
    
    if cleanup:
        system("rm viz_output.txt")
        #system("rm viz_output_cropped.txt")
        system("rm " + progenitor_image_name.replace(".fits", ".sex"))
        system("rm " + weight_image_name.replace(".fits", ".sex"))
    
     # Outputting the .reg file of the calibration stars (if calstar_reg_output==True).

    if autocull:
        # Write new calibration star calregion file
        tempname_unique = str(progenitor_image_name.split('_')[2]) # + str(progenitor_image_name.split('_')[4]).rstrip('.fits') 
        temppath = storepath + 'calstarregs/' + tempname_unique+'_autocull.reg'
        # Remove old autocull file if autocull is set to True
        if autocull:
            autocull_temp = './autoculltemp.reg'
            autocullpath = autocull_temp
            if os.path.isfile(autocullpath) == True:
                os.remove(autocullpath)
            
        tempreg = open(autocull_temp, 'w')
        tempreg.write('# Region file format: DS9 version 4.1\n')
        secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                                       '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                                       'delete=1 include=1 source=1\n'
        tempreg.write(secondstr)
        tempreg.write('fk5\n')
        # Add a few arcseconds to the first position to make sure we don't use it as a target
        #Why?
#        test_ra = float(calreg_stars.split(',')[0]) + 0.005
#        test_dec = float(calreg_stars.split(',')[1]) + 0.005
        print 'calreglist final is'
        print calreglist
        print 'length of calreglist is'
        print len(calreglist)
        for reg_lines in calreglist:
            reg_ra = float(reg_lines.split(',')[0])
            reg_dec = float(reg_lines.split(',')[1])
            tmp_reg_str = 'circle(%f,%f,4") # width=2 font="helvetica 16 normal"\n' % (reg_ra,reg_dec)                 
            tempreg.write(tmp_reg_str)       
        tempreg.close()
        if os.path.isfile(temppath) == True:
                os.remove(temppath)
        command = 'mv ./autoculltemp.reg' + ' '+temppath
        os.system(command)
        print 'Autoculled calibration star calregion created, please re-run photometry using this calibration file'

    if calstar_reg_output == True:
        uniquename = str(progenitor_image_name.split('_')[2]) # + str(progenitor_image_name.split('_')[4]).rstrip('.fits') 
        reg_name = storepath + 'calstarregs/'  + uniquename + '.reg'
        regfile=open(reg_name,'w')
        regfile.write('# Region file format: DS9 version 4.1\n')
        secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
        regfile.write(secondstr)
        regfile.write('fk5\n')
        for calstar in combined_starlist:
             #only add it to the region file if it is not our target
            if calstar[2] != 999:
                tmp_str = 'circle('+str(calstar[0])+','+str(calstar[1]) +','+str(aperture_size)+'")' + '\n'
                regfile.write(tmp_str)
        regfile.close
    print ("Photometry completed, " + 
        str(time() - t1) + " seconds required.")
        
    return photdict
    

def do_dir_phot(photdir='./',reg='PTEL.reg',ap=None, do_upper=False, auto_upper=True):
    '''
    if auto_upper = True, then rerun to find upper limit if no target mag found
    '''
    import glob
    mylist = []
    photdict = {}
    fulllist = glob.glob(photdir+'/*coadd*.fits')
    # Remove the weight images from the list
    for item in fulllist:
        if item.find('weight') == -1:
            mylist.append(item)
    for myfile in mylist:
        print "Now performing photometry for %s \n" % (myfile)
        photout = dophot(myfile,reg,ap)
        # If a target magnitude or upper limit isn't found, rerun the 
        # photometry code with do_upper = True to find an upper limit
        if 'targ_mag' not in photout and 'upper_green' not in photout and auto_upper:
            print '**Target Magnitude not found. Re-running to find UL**.'
            photout = dophot(myfile,reg,ap,do_upper=True)
        label = photout['FileName']
        photdict.update({label:photout})
    return photdict


def openCalRegion(reg_path, nonstore=False):
    if nonstore:
        reg_path = reg_path
    else:
        reg_path = storepath + reg_path
    regfile = open(reg_path,'r')
    reglist = regfile.readlines()
    callist = []
    for line in reglist:
        if 'circle' in line:
            callist += [(((line.rstrip('\n')).strip('circle')).rstrip(')')).strip('(')]
        else:
            pass
   
    return callist

def photreturn(GRBname, filename, clobber=False, reg=None, aper=None, \
    auto_upper=True, calregion = None, trigger_id = None, stardict = None, 
    caliblimit=True, verbose=False, autocull=False, find_fwhm=False,
    fallback_n_dither=0,reuse_old_vizcat=False):
    '''Returns the photometry results of a GRB that was stored in a pickle file. 
    If the pickle file does not exists, this function will create it. Use 
    clobber=True for overwriting existing pickle files. '''
    filepath = storepath + GRBname + 'ap' + str(aper)  
    if calregion:
        calibration_list = openCalRegion(calregion)
        n_calstars = len(calibration_list)
        filepath += '_WithCalReg' + str(n_calstars)
    if stardict:
        filepath += '_WithDeepStack'
    filepath += '.data'
    while clobber == False:
        if os.path.isfile(filepath) == True:
            data = qPickle.load(filepath)
            if filename in data:
                return data
            else:
                clobber = True
        else:
            clobber = True
    
    while clobber == True:
        if reg == None: 
            print 'Need to input reg file'
            break
        else:
            if os.path.isfile(filepath) == False:
                photdict = {}  
            else:
                #f = file(filepath)
                photdict = qPickle.load(filepath) # This line loads the pickle file, enabling photLoop to work                
            data = dophot(filename, reg, ap=aper, calreg = calregion, stardict=stardict, 
                caliblimit=caliblimit, verbose=verbose, autocull=autocull, find_fwhm = find_fwhm, 
                fallback_n_dither=fallback_n_dither, reuse_old_vizcat=reuse_old_vizcat)
            if 'targ_mag' not in data and 'upper_green' not in data and auto_upper:
                print '**Target Magnitude not found. Re-running to find UL**.'
                data = dophot(filename,reg,aper,do_upper=True, stardict=stardict, calreg=calregion, caliblimit=caliblimit, verbose=verbose, find_fwhm = find_fwhm, fallback_n_dither=fallback_n_dither)
            label = data['FileName']
            time = float(t_mid.t_mid(filename, trigger = trigger_id))
            terr = float(t_mid.t_mid(filename, delta = True, trigger = trigger_id))/2.
            timetuple = (time, terr)
            data.update({'t_mid':timetuple})
            photdict.update({label:data})
            qPickle.save(photdict, filepath, clobber = True)
            return photdict
            break
    
def photLoop(GRBname, regfile, ap=None, calregion = None, trigger_id = None, \
    stardict=None, clobber=False, auto_upper=True, caliblimit=True, verbose=False, 
    autocull=False, find_fwhm=False, global_fallback_n_dither=0, 
    reuse_vizcat_after_one_iteration=True,force_continue_if_error=True):
    '''Run photreturn on every file in a directory; return a dictionary
    with the keywords as each filename that was observed with photreturn
    '''
    import glob
    GRBlist = []
    GRBlistwweight = glob.glob('*.fits')
    # Remove the weight images from the list
    filepath = storepath + GRBname + '.data'
    filepath = storepath + GRBname + 'ap' + str(ap)
    if calregion:
        calibration_list = openCalRegion(calregion)
        n_calstars = len(calibration_list)
        filepath += '_WithCalReg' + str(n_calstars)
    if stardict:
        filepath += '_WithDeepStack'
    filepath += '.data'
    # Remove old autocull file if autocull is set to True
    if autocull: 
        autocullpath = storepath + 'calstarregs/' + GRBname+'_autocull.reg'
        if os.path.isfile(autocullpath) == True:
            os.remove(autocullpath)
    # Remove the stored pickle file if clobber is set to True
    if clobber == True:
        if os.path.isfile(filepath) == True:
            os.remove(filepath)
    for item in GRBlistwweight:
        if item.find('weight') == -1:
            GRBlist.append(item)
    if not GRBlist:
        raise ValueError('No files to perform photometry on. In right directory?')
    
    reuse=False
    count=0
    for mosaic in GRBlist:
        if count > 0 and reuse_vizcat_after_one_iteration: reuse = True
        print "=========================================================="
        print "Now performing photometry for %s \n" % (mosaic)
        try:
            photout = photreturn(GRBname, mosaic, clobber=clobber, reg=regfile, \
            aper=ap, calregion = calregion, trigger_id=trigger_id, stardict=stardict,
            auto_upper=auto_upper, caliblimit=caliblimit, verbose=verbose, autocull=autocull, 
            find_fwhm=find_fwhm, fallback_n_dither=global_fallback_n_dither,
            reuse_old_vizcat=reuse)
        except:
            if force_continue_if_error:
                print '****Error was encountered, but continuing to next image****'
                pass
            else:
                raise Exception
        count += 1
    return photout
    
def plotzp(photdict,ylim=None):
    '''Plots a graph of zp from the pickle output of the photreturn function'''
    import matplotlib
    import glob
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    vallist = []
    errlist = []

    for mosaics in photdict:
        print 'now doing ' + str(mosaics)
        time = photdict[mosaics]['t_mid'][0]
        terr = photdict[mosaics]['t_mid'][1]
       
        valu = float(photdict[mosaics]['zp'][0])
        verr = float(photdict[mosaics]['zp'][1])

#there's probably a prettier way to do this, the second if statements are 
#there so that only 1 label per filter is on the legend

        if 'h_' in mosaics:
            if h == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = '.', linestyle ='None', mfc = 'green', mec = 'green',\
                    ecolor = 'green')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = '.', linestyle ='None', mfc = 'green', mec = 'green',\
                    ecolor = 'green', label = 'h')
                h = True

        elif 'j_' in mosaics:            
            if j == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = '.', linestyle ='None', mfc = 'blue', mec = 'blue', \
                    ecolor = 'blue')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = '.', linestyle ='None', mfc = 'blue', mec = 'blue', \
                    ecolor = 'blue', label = 'j')
                j = True

        elif 'k_' in mosaics:
            if k == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = '.', linestyle ='None', mfc = 'red', \
                    mec = 'red', ecolor = 'red')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = '.', linestyle ='None', mfc = 'red', \
                    mec = 'red', ecolor = 'red', label = 'k')
                k = True

    ax = matplotlib.pyplot.gca()
    if ylim:
        ax.set_ylim(ylim)
    ax.set_ylim(ax.get_ylim()[::-1])
    
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('zero point')
    matplotlib.pyplot.legend()
    uniquename = photdict.keys()[0].split('_')[2]
    savepath = storepath + uniquename + '_zeropoint.png'
    print 'zeropoint plot saved to ' + savepath
    savefig(savepath)           
    matplotlib.pyplot.close()


def plots2n(photdict):
    '''Plots a graph of s/n versus duration from the pickle output of the 
    photreturn function. This function also takes the output from photLoop'''
    import matplotlib
    import glob
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    vallist = []
    errlist = []

    for mosaics in photdict:
        print 'now doing ' + str(mosaics)
        #time = float(photdict[mosaics]['EXPTIME'])
        time = photdict[mosaics]['t_mid'][0]
       
        valu = float(photdict[mosaics]['targ_s2n'])

#there's probably a prettier way to do this, the second if statements are 
#there so that only 1 label per filter is on the legend

        if 'h_' in mosaics:
            if h == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'red', mec = 'green', \
                    ecolor = 'red')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'red', mec = 'green', \
                    ecolor = 'red', label = 'h')
                h = True

        elif 'j_' in mosaics:            
            if j == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'blue', mec = 'green', \
                    ecolor = 'blue')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'blue', mec = 'green', \
                    ecolor = 'blue', label = 'j')
                j = True

        elif 'k_' in mosaics:
            if k == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'yellow', mec = 'green',\
                    ecolor = 'yellow')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = 0, xerr= 0, \
                    marker = 's', linestyle ='None', mfc = 'yellow', mec = 'green',\
                    ecolor = 'yellow', label = 'k')
                k = True
    
    matplotlib.pyplot.xlabel('duration(s)')
    matplotlib.pyplot.ylabel('s/n')
    matplotlib.pyplot.legend()
    #matplotlib.pyplot.semilogx()
    uniquename = photdict.keys()[0].split('_')[2]
    #savepath = storepath + uniquename + '_s2nvsdurationlogx.png'
    savepath = storepath + uniquename + '_s2n.png'
    print 's/n plot saved to ' + savepath
    savefig(savepath)           
    matplotlib.pyplot.close()

def plotproperty(photdict, prop):
    '''Plots a graph of a property from the pickle output of the photreturn function
    CURRENTLY NOT IMPLEMENTED?
    '''
    import matplotlib
    import glob

    matplotlib.pyplot.clf()
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    vallist = []
    errlist = []

    for mosaics in photdict:
        print 'now doing ' + str(mosaics)
        time = photdict[mosaics]['t_mid'][0]
        terr = photdict[mosaics]['t_mid'][1]

        if prop == 'FWHM':
            valu = float(photdict[mosaics]['FWHM'][0])
            verr = float(photdict[mosaics]['FWHM'][1])

#there's probably a prettier way to do this, the second if statements are 
#there so that only 1 label per filter is on the legend

        if 'h_' in mosaics:
            if h == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = 'o', linestyle ='None', mfc = 'green', mec = 'green',\
                    ecolor = 'green')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = 'o', linestyle ='None', mfc = 'green', mec = 'green',\
                    ecolor = 'green', label = 'h')
                h = True

        elif 'j_' in mosaics:            
            if j == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                    ecolor = 'blue')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, \
                    marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                    ecolor = 'blue', label = 'j')
                j = True

        elif 'k_' in mosaics:
            if k == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = 'o', linestyle ='None', mfc = 'red', \
                    mec = 'green', ecolor = 'red')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr,\
                    marker = 'o', linestyle ='None', mfc = 'red', \
                    mec = 'green', ecolor = 'red', label = 'k')
                k = True

    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    if prop == 'FWHM':
        axis_label = 'FWHM (arcseconds)'
    matplotlib.pyplot.ylabel(axis_label)
    matplotlib.pyplot.legend()
    uniquename = photdict.keys()[0].split('_')[2]
    propname = '_'+prop+'.png'
    savepath = storepath + uniquename + propname
    print '%s plot saved to ' % (prop) + savepath
    savefig(savepath)           
    matplotlib.pyplot.close()

def oldtextoutput(dict,filt=None):
    '''outputs a text file from photdict.  If filt specified, only output that
    particular filter.
    '''
    filts = ['j','h','k']
    if filt:
        if not filt in filts:
            raise ValueError('Unacceptable value for filt. Needs to be j, h, or k')
            
    uniquename = dict.keys()[0].split('_')[2]
    if filt:
        uniquename = uniquename + '_' + filt
    savepath = storepath + uniquename + '_lightcurve.txt'
    text = file(savepath, "w")
    namelist = ['FileName', 'Aperture', 'STRT_CPU', 'STOP_CPU', 't_mid' , 'EXPTIME',\
        'targ_mag', 'targ_mag error']
    text.write('\t'.join(namelist))
    text.write('\n')
    
    timelist = []
    for key in dict:
        timelist.append((dict[key]['t_mid'][0], key))
    timelist.sort()
    
    for tup in timelist:
        values = dict[tup[1]]
        if filt:
            if filt != values['filter']:
                continue
        if 'targ_mag' not in values:
            try:
                mag = str(values['upper_green'])
            except:
                print 'Neither Targ Mag nor Upper limit found!'
                print tup                
            magerr = '3sig'
        else:
            mag = str(values['targ_mag'][0]) 
            magerr = str(values['targ_mag'][1])
        time = str(values['t_mid'][0])
        datalist = [str(values['FileName']), str(values['Aperture']), str(values['STRT_CPU']), \
            str(values['STOP_CPU']), time, str(values['EXPTIME']), mag, magerr]
        text.write('\t'.join(datalist))
        text.write('\n')
    text.close()


def tmp_phot_plot(photdict):
    '''quick write-up of getting the photometry because I am lazy '''
    import pylab
    timearr = []
    timeerrarr = []
    photarr = []
    photerrarr = []
    for imgresult in photdict.itervalues():
        timearr.append(imgresult['HJD_mid'])
        duration = imgresult['HJD_start'] - imgresult['HJD_stop']
        timeerrarr.append(duration)
        photarr.append(imgresult['targ_mag'][0])
        photerrarr.append(imgresult['targ_mag'][1])
        
    pylab.errorbar(y=photarr,x=timearr,xerr=timeerrarr,yerr=photerrarr,fmt='ro')
    pylab.ylim(pylab.ylim()[1],pylab.ylim()[0]) #transpose axes
    pylab.xlabel('Time since Burst (s)')
    pylab.ylabel('J Mag')
    pylab.semilogx()


def photplot(photdict,ylim=None,xlim=None, plotname=False, overplot=False, no_ap=False, 
            ylab="Mag", xlab="Time since Burst (s)", divx=1.0, divy=1.0, logx=True,
            title=False, round_x_offset=False):
    '''Plots a graph from photdict
    
    title: override the default title with one of your own
    ylab: Y axis label
    xlab: X axis label
    divx: factor to divide the X values by; useful if e.g. want to plot days instead of seconds.
    divy: Factor to divide the Y values by
    logx: take the log of the X values?
    round_x_offset: make the offset of the tickmarks a more resonable round value (useful for when using mjd)
    '''
    import matplotlib
    import glob
    import datetime
    
    #filepath = path + GRBname + '.data'
    #photdict = qPickle.load(filepath)

    matplotlib.pyplot.clf()
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    vallist = []
    errlist = []

    ap = None
    

    for mosaics in photdict:
        print 'now doing ' + str(mosaics)
        if ap:
            if photdict[mosaics]['Aperture'] != ap:
                if no_ap:
                    pass
                else:
                    raise ValueError ('The aperture used in this mosaic does not match the last one.')
        if no_ap:
            pass
        else:
            ap = photdict[mosaics]['Aperture']
        time = photdict[mosaics]['t_mid'][0]/divx
        terr = photdict[mosaics]['t_mid'][1]/divx
        
        if 'targ_mag' in photdict[mosaics]: 
            valu = float(photdict[mosaics]['targ_mag'][0])/divy
            verr = float(photdict[mosaics]['targ_mag'][1])/divy

            #there's probably a prettier way to do this, the second if statements 
            #are there so that only 1 label per filter is on the legend

            if 'h_' in mosaics:
                if h == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green', label = 'h')
                    h = True

            elif 'j_' in mosaics:            
                if j == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue', label = 'j')
                    j = True

            elif 'k_' in mosaics:
                if k == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'k')
                    k = True
        elif 'upper_green' in photdict[mosaics]: 
                valu = float(photdict[mosaics]['upper_green'])
                if photdict[mosaics]['filter'] == 'h':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                if photdict[mosaics]['filter'] == 'j':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                if photdict[mosaics]['filter'] == 'k':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
    if round_x_offset:
        locs,labels = matplotlib.pyplot.xticks()
        matplotlib.pyplot.xticks(locs, map(lambda x: "%g" % x, locs-min(locs)))
        matplotlib.pyplot.text(0.92, -0.07, "+%g" % min(locs), fontsize=10, transform = matplotlib.pyplot.gca().transAxes)
    
    
    matplotlib.pyplot.xlabel(xlab)
    matplotlib.pyplot.ylabel(ylab)
    #matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    if logx:
        ax.set_xscale('log')
    matplotlib.pyplot.legend()
    uniquename = photdict.keys()[0].split('_')[2]
    if title:
        # override the unique name with the custom title.
        uniquename = title
    if no_ap:
        matplotlib.pyplot.title(uniquename)
    else:
        matplotlib.pyplot.title(uniquename+' Lightcurve: ap = ' + str(ap))
    savepath = storepath + uniquename + '_' + 'ap' + str(ap) + '_lightcurve.eps'
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if overplot:
        pass
    else:
        F = pylab.gcf()
        DefaultSize = F.get_size_inches()
        DPI = F.get_dpi()
        # F.set_size_inches( (DefaultSize[0]*2.5, DefaultSize[1]*2.5) )
        # was getting incresingly larger with multiple runs
        F.set_size_inches((20, 15))

    if plotname:
        savepath = plotname

    if overplot:
        pass
    else:
        print 'lightcurve saved to ' + savepath
        savefig(savepath)    
        matplotlib.pyplot.close()

def findOptimalAperture(GRBname, regfile, calregion=None, trigger_id = None, plot=True, 
        caliblimit=True, global_fallback_n_dither=0):
    ''' Due to sampling and dithering issues, finding the optimal aperture in
    PAIRITEL is not as simple as simply measuring the FWHM of an image.  Here,
    we loop around images in a directory using PhotLoop and measure the 
    photometry of the calibration stars.
    '''
    import pylab
    pylab.clf()
    j_delta_med_list = []
    h_delta_med_list = []
    k_delta_med_list = []
    aplist = [2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0]
    for ap in aplist:
        
        # Can set auto_upper = False since dont care if source is detected or not
        data = photLoop(GRBname,regfile,calregion=calregion, ap=ap, \
            clobber=True, trigger_id=trigger_id, auto_upper=False, caliblimit=caliblimit,
            global_fallback_n_dither=global_fallback_n_dither) 
    
        k_delta_list = []
        h_delta_list = []
        j_delta_list = []

        for key, val in data.iteritems():
            for calkey, calval in val['calib_stars'].iteritems():
                delta = abs(calval['delta_mag'])
                if val['filter']=='j':
                    j_delta_list.append(delta)
                elif val['filter']=='h':
                    h_delta_list.append(delta)
                elif val['filter']=='k':
                    k_delta_list.append(delta)
                else:
                    raise ValueError
    
        j_delta_med_list.append(numpy.median(j_delta_list))
        h_delta_med_list.append(numpy.median(h_delta_list))
        k_delta_med_list.append(numpy.median(k_delta_list))
    
    if plot:
        pylab.scatter(aplist,j_delta_med_list, color='b')
        pylab.scatter(aplist,h_delta_med_list, color='g')
        pylab.scatter(aplist,k_delta_med_list, color='r')
        
        pylab.legend(['j','h','k'])
        pylab.xlabel('Aperture Size (")')
        pylab.ylabel('Median Calibration star Uncertainty (mag)')
        pylab.title('Aperture Optimization for ' + GRBname)
        
        savepath = storepath + GRBname + 'ApFind.png'
        savefig(savepath)           
        matplotlib.pyplot.close()
        print 's/n plot saved to ' + savepath

    return (j_delta_med_list,h_delta_med_list,k_delta_med_list)

def colorevo(photdict, photdict_lcurve, JorK, ylim=None,xlim=None, big=False):
    '''Plots a color evolution J-H, H-K from photdict. There are two photdicts, one for the color evolution, and another for the lightcurve. This is because they might have different time resolution. Of course, they can be the same.'''
    import matplotlib
    import glob
    import datetime
    from operator import itemgetter
    #filepath = path + GRBname + '.data'
    #photdict = qPickle.load(filepath)

    matplotlib.pyplot.clf()
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    h_list = []
    h_errlist = []
    j_list = []
    j_errlist = []
    k_list = []
    k_errlist = []
    mosaiclist=[]
    for mosaics in photdict:
        mosaiclist += [photdict[mosaics]]
    #sorting w.r.t time
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)
    for mosaics in mosaiclist:

        if 'targ_mag' in mosaics: 
            valu = float(mosaics['targ_mag'][0])
            verr = float(mosaics['targ_mag'][1])
            
            if mosaics['filter'] == 'h':                
                h_list += [valu]
                h_errlist += [verr]
                timlist += [mosaics['t_mid'][0]]
                terlist += [mosaics['t_mid'][1]]

            elif mosaics['filter'] == 'j':
                j_list += [valu]
                j_errlist += [verr]
            elif mosaics['filter'] == 'k':
                k_list += [valu]
                k_errlist += [verr]
            
        elif 'upper_green' in mosaics: 
                valu = float(mosaics['upper_green'])
                if mosaics['filter'] == 'h':
                     h_list += [valu]
                     h_errlist += [0]
                     timlist += [mosaics['t_mid'][0]]
                     terlist += [mosaics['t_mid'][1]]
                if mosaics['filter'] == 'j':
                     j_list += [valu]
                     j_errlist += [0]
                if mosaics['filter'] == 'k':
                     k_list += [valu]
                     k_errlist += [0]
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)

    timlist = numpy.array(timlist)
    terlist = numpy.array(terlist)
    h_list = numpy.array(h_list)
    j_list = numpy.array(j_list)
    k_list = numpy.array(k_list)
    h_errlist = numpy.array(h_errlist)
    j_errlist = numpy.array(j_errlist)
    k_errlist = numpy.array(k_errlist)

    jminh = j_list - h_list
    hmink = h_list - k_list

    fig = pylab.figure()

    ax_lcurve = fig.add_axes([0.1,0.55,0.8,0.4])
    pylab.subplots_adjust(hspace=0.001)    
    #Lightcurve
    
    h = False
    j = False
    k = False
    ap = None
    
    for mosaics in photdict_lcurve:
        print 'now doing ' + str(mosaics)
        if ap:
            if photdict_lcurve[mosaics]['Aperture'] != ap:
                raise ValueError ('The aperture used in this mosaic does not match the last one.')
        ap = photdict_lcurve[mosaics]['Aperture']
        time = photdict_lcurve[mosaics]['t_mid'][0]
        terr = photdict_lcurve[mosaics]['t_mid'][1]
        
        if 'targ_mag' in photdict_lcurve[mosaics]: 
            valu = float(photdict_lcurve[mosaics]['targ_mag'][0])
            verr = float(photdict_lcurve[mosaics]['targ_mag'][1])

            #there's probably a prettier way to do this, the second if statements 
            #are there so that only 1 label per filter is on the legend

            if 'h_' in mosaics:
                if h == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green', label = 'h')
                    h = True

            elif 'j_' in mosaics:            
                if j == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue', label = 'j')
                    j = True

            elif 'k_' in mosaics:
                if k == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'k')
                    k = True
        elif 'upper_green' in photdict_lcurve[mosaics]: 
                valu = float(photdict_lcurve[mosaics]['upper_green'])
                if photdict_lcurve[mosaics]['filter'] == 'h':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                if photdict_lcurve[mosaics]['filter'] == 'j':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                if photdict_lcurve[mosaics]['filter'] == 'k':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
#    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    #matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    matplotlib.pyplot.legend()
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    uniquename = photdict_lcurve.keys()[0].split('_')[2]
    matplotlib.pyplot.title(uniquename+' Lightcurve: ap = ' + str(ap))

    # Color evolution plot    

    ax = fig.add_axes([0.1,0.1,0.8,0.3])
    
    if JorK == 'j':
        error_array = sqrt(j_errlist**2 + h_errlist**2)
        ax.errorbar(timlist, jminh, yerr = error_array, xerr=terlist,\
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'j-h')
    elif JorK =='k':
        error_array = sqrt(h_errlist**2 + k_errlist**2)
        ax.errorbar(timlist, hmink, yerr = error_array, xerr=terlist,\
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue', label = 'h-k')
    elif JorK == 'both':
        jh_error_array = sqrt(j_errlist**2 + h_errlist**2)
        hk_error_array = sqrt(h_errlist**2 + k_errlist**2)
        ax.errorbar(timlist, jminh, yerr=jh_error_array, xerr=terlist,\
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'j-h')
        ax.errorbar(timlist, hmink, yerr=hk_error_array,xerr=terlist,\
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue', label = 'h-k')
        
    #ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    #ax.set_xscale('log')
    matplotlib.pyplot.legend()
    uniquename = photdict.keys()[0].split('_')[2]
    if JorK == 'j':
        cname = 'J-H'
    elif JorK == 'k':
        cname = 'H-K'
    elif JorK == 'both':
        cname = 'both'
    #matplotlib.pyplot.title(uniquename+' Color Evolution '+'('+cname+')')
    savepath = storepath + uniquename + '_' + cname + '_colorevo.png'
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
        
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    DPI = F.get_dpi()
    if big:
        F.set_size_inches((20, 15))
    
    print 'colorevolution saved to ' + savepath
    savefig(savepath)    
    matplotlib.pyplot.close()


def textoutput(photdict,utburst='00:00:00',name=None,filt=None, day=False):
    '''outputs a text file from photdict.  If filt specified, only output that
    particular filter. Requires a utburst string in the form of hh:mm:ss. 
    Days are integers (e.g. second day observation => day=2. WARNING: If the 
    observations consists of multiple days, care needs to be exercised! 
    Current fix: do it manually!
    '''

    import datetime
    from operator import itemgetter
    
    filts = ['j','h','k']
    if filt:
        if not filt in filts:
            raise ValueError('Unacceptable value for filt. Needs to be j, h, or k')
    if name:
        uniquename = name
    else:  
        uniquename = photdict.keys()[0].split('_')[2]
    if filt:
        uniquename = uniquename + '_' + filt
    savepath = storepath + uniquename + '_data.dat'
    text = file(savepath, "w")

    text.write('@inunit=days')
    text.write('\n')
    text.write('@expunit=sec')
    text.write('\n')
    utburst_text = '@utburst=' + str(utburst)
    text.write(utburst_text)
    text.write('\n')
    text.write('@source=PAIRITEL')
    text.write('\n')
    namelist = ['%', 't_mid', 'ut_start_date', 'ut_start_time', 'ut_end_date', 'ut_end_time', 'exp', 'filt', '=', 'mag', 'emag', 'lim']
    text.write(' '.join(namelist))
    text.write('\n')

    mosaiclist=[]
    for mosaics in photdict:
        mosaiclist += [photdict[mosaics]]
    #sorting w.r.t time
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)

    h_list = []
    her_list = []
    j_list = []
    jer_list = []
    k_list = []
    ker_list = []
    tstart_list = []
    tstop_list = []
    exp_list = []

    burst_time_sec = float(utburst.split(':')[0])*3600 + float(utburst.split(':')[1])*60 + float(utburst.split(':')[2]) 
    
    for mosaics in mosaiclist:
        
        strt_cpu = mosaics['STRT_CPU']
        stop_cpu = mosaics['STOP_CPU']
        t_mid_days = mosaics['t_mid'][0]/86400.
        
        timestart_str = mosaics['STRT_CPU'].split(' ')[1]
        timestart_sec = float(timestart_str.split(':')[0])*3600. + float(timestart_str.split(':')[1])*60. + float(timestart_str.split(':')[2])
        if day:
            start_after_burst_sec = timestart_sec - burst_time_sec + (day-1)*86400
        else:
            start_after_burst_sec = timestart_sec - burst_time_sec
        timestop_str = mosaics['STOP_CPU'].split(' ')[1]
        timestop_sec = float(timestop_str.split(':')[0])*3600. + float(timestop_str.split(':')[1])*60. + float(timestop_str.split(':')[2])
        if day:
            stop_after_burst_sec = timestop_sec - burst_time_sec + (day-1)*86400
        else:
            stop_after_burst_sec = timestop_sec - burst_time_sec
        exp = mosaics['EXPTIME']
        
        if 'targ_mag' in mosaics: 
            valu = float(mosaics['targ_mag'][0])
            verr = float(mosaics['targ_mag'][1])
            
            # IS IT NECESSARY TO REPEAT THESE EACH TIME? JUST GRAB FILT DIRECTLY
            if mosaics['filter'] == 'h':                
                mag = valu
                magerr = verr
                filt = 'H'
            elif mosaics['filter'] == 'j':
                mag = valu
                magerr = verr
                filt = 'J'
            elif mosaics['filter'] == 'k':
                mag = valu
                magerr = verr
                filt = 'K'
            
        elif 'upper_green' in mosaics: 
                valu = float(mosaics['upper_green'])
                if mosaics['filter'] == 'h':
                     mag = valu
                     magerr = 0
                     filt = 'H'
                if mosaics['filter'] == 'j':
                     mag = valu
                     magerr = 0
                     filt = 'J'
                if mosaics['filter'] == 'k':
                     mag = valu
                     magerr = 0
                     filt = 'K'
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)    
 
        if 'upper_green' in mosaics:
            datalist = [str(t_mid_days), str(strt_cpu), str(stop_cpu), str(exp), filt, '=', str(mag), str(magerr), 'yes']
        else:   
            datalist = [str(t_mid_days), str(strt_cpu), str(stop_cpu), str(exp), filt, '=', str(mag), str(magerr)]
        # old    
        # if 'upper_green' in mosaics:
        #     datalist = [str(start_after_burst_sec), str(stop_after_burst_sec), str(exp), filt, '=', str(mag), str(magerr), 'yes']
        # else:   
        #     datalist = [str(start_after_burst_sec), str(stop_after_burst_sec), str(exp), filt, '=', str(mag), str(magerr)]
        # 

        text.write(' '.join(datalist))
        text.write('\n')
    text.close()

def plotindex(photdict, photdict_lcurve, ylim=None,xlim=None, big=False, dif_ap=False, just_two=False):
    '''Plots spectral indexes vs. time from photdict. If there is no detection in the k band, set just_two = jh'''
    import matplotlib
    import glob
    import datetime
    import numpy
    from operator import itemgetter
    from scipy import optimize
    from scipy import log10
    from MiscBin import q
    from Phot import pofit

    #filepath = path + GRBname + '.data'
    #photdict = qPickle.load(filepath)

    matplotlib.pyplot.clf()
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    h_list = []
    h_errlist = []
    j_list = []
    j_errlist = []
    k_list = []
    k_errlist = []
    mosaiclist=[]
    for mosaics in photdict:
        mosaiclist += [photdict[mosaics]]
    #sorting w.r.t time
    get = itemgetter('t_mid')
    mosaiclist.sort(key=get)
    for mosaics in mosaiclist:

        if 'targ_mag' in mosaics: 
            valu = float(mosaics['targ_mag'][0])
            verr = float(mosaics['targ_mag'][1])
            
            if mosaics['filter'] == 'h':                
                h_list += [valu]
                h_errlist += [verr]
                timlist += [mosaics['t_mid'][0]]
                terlist += [mosaics['t_mid'][1]]

            elif mosaics['filter'] == 'j':
                j_list += [valu]
                j_errlist += [verr]
            elif mosaics['filter'] == 'k':
                k_list += [valu]
                k_errlist += [verr]
            
        elif 'upper_green' in mosaics: 
                valu = float(mosaics['upper_green'])
                if mosaics['filter'] == 'h':
                     h_list += [valu]
                     h_errlist += [0] # this error = 0 causes optiphot to fail at fitting SED
                     timlist += [mosaics['t_mid'][0]]
                     terlist += [mosaics['t_mid'][1]]
                if mosaics['filter'] == 'j':
                     j_list += [valu]
                     j_errlist += [0]
                if mosaics['filter'] == 'k':
                     k_list += [valu]
                     k_errlist += [0]
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)

    timlist = numpy.array(timlist)
    terlist = numpy.array(terlist)
    h_list = numpy.array(h_list)
    j_list = numpy.array(j_list)
    k_list = numpy.array(k_list)
    h_errlist = numpy.array(h_errlist)
    j_errlist = numpy.array(j_errlist)
    k_errlist = numpy.array(k_errlist)

    # We have to fit these one by one.
    # Note that the length of h, j, k must be the same.

    betalist = []
    notime_list = []
    for index, mosaic in enumerate(j_list):
        if just_two == 'jh':
            h_mag = h_list[index]
            j_mag = j_list[index] # sort by ascending frequencies
        elif just_two == 'kh':
            k_mag = k_list[index]
            h_mag = h_list[index] # sort by ascending frequencies
        else:
            k_mag = k_list[index]
            h_mag = h_list[index]
            j_mag = j_list[index] # sort by ascending frequencies

        # These aren't fluxes.. they are converted to AB magnitude
        #k_flux = -2.5*numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) - 48.60 # zeropoints are in erg/s/cm^2/Hz, convert to AB mag
        #h_flux = -2.5*numpy.log10(q.mag2flux(h_mag,0,1.024e-20)) - 48.60
        #j_flux = -2.5*numpy.log10(q.mag2flux(j_mag,0,1.594e-20)) -48.60
        if just_two == 'jh':
            h_flux = numpy.log10(q.mag2flux(h_mag,0,1.024e-20))
            j_flux = numpy.log10(q.mag2flux(j_mag,0,1.594e-20))
        elif just_two == 'kh':
            k_flux = numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) 
            h_flux = numpy.log10(q.mag2flux(h_mag,0,1.024e-20))
        else:
            k_flux = numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) # zeropoints are in erg/s/cm^2/Hz, convert to AB mag
            h_flux = numpy.log10(q.mag2flux(h_mag,0,1.024e-20))
            j_flux = numpy.log10(q.mag2flux(j_mag,0,1.594e-20)) 
#        k_flux = q.mag2flux(k_mag,0,666.7e6) # zeropoints are in microjansky
#        h_flux = q.mag2flux(h_mag,0,1024e6)
#        j_flux = q.mag2flux(j_mag,0,1594e6) 

        if just_two == 'jh':
            specmag = numpy.array([h_flux, j_flux])
            specmag_error = numpy.array([h_errlist[index], j_errlist[index]])
            p_freqs = numpy.array([1.816923988e+14,2.398339664e+14]) 
            log_freqs = log10(p_freqs) 
            log_specmag_err = numpy.array(specmag_error)/numpy.array(specmag)

        elif just_two == 'kh':
            specmag = numpy.array([k_flux, h_flux])
            specmag_error = numpy.array([k_errlist[index], h_errlist[index]])
            p_freqs = numpy.array([1.394383526e+14,1.816923988e+14]) 
            log_freqs = log10(p_freqs) 
            log_specmag_err = numpy.array(specmag_error)/numpy.array(specmag)

        else:
            specmag = numpy.array([k_flux, h_flux, j_flux])
            specmag_error = numpy.array([k_errlist[index], h_errlist[index], j_errlist[index]])
            p_freqs = numpy.array([1.394383526e+14,1.816923988e+14,2.398339664e+14]) 
            log_freqs = log10(p_freqs) 
            log_specmag_err = numpy.array(specmag_error)/numpy.array(specmag)

#        fitfunc = lambda p, x: p[0] + p[1] * x # the linear function in log-log space
#        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err # the error function, this is just the data minus the line of best fit divided by the error
#        pinit = [0., -1.0] # initial guess

#        out = optimize.leastsq(errfunc, pinit, args=(log_freqs, specmag, log_specmag_err), full_output=1)
#        pfinal = out[0] # Getting the coefficients
#        covar = out[1] # The covariance        

        print '----'
        print 'index number %s' % index
        if just_two == 'jh':
            print h_flux
            print j_flux
        elif just_two == 'kh':
            print k_flux
            print h_flux
        else:
            print k_flux
            print h_flux
            print j_flux
        print specmag
        print '----'   

        print 'specmag _error is'
        print specmag_error
        try:
            the_fit = pofit.fit(p_freqs, specmag, specmag_error, 'beta', name='spectral fit')
        except:
            notime_list += [index]
            print 'cannot fit!'
            print 'specmag_error is:'
            print specmag_error
            print 'Go investigate!'
            break
        results = the_fit[0]
        betalist += [results[1]] # Getting the spectral index

    fig = pylab.figure()
    ax_lcurve = fig.add_axes([0.1,0.55,0.8,0.4])
    pylab.subplots_adjust(hspace=0.001)
    
    #Lightcurve
    
    h = False
    j = False
    k = False
    ap = None
    
    for mosaics in photdict_lcurve:
        print 'now doing ' + str(mosaics)
        if ap:
            if not dif_ap:
                if photdict_lcurve[mosaics]['Aperture'] != ap:
                    raise ValueError ('The aperture used in this mosaic does not match the last one.')
                else:
                    print 'The aperture used in this mosaic does not match the last one, continuing...'
        ap = photdict_lcurve[mosaics]['Aperture']
        time = photdict_lcurve[mosaics]['t_mid'][0]
        terr = photdict_lcurve[mosaics]['t_mid'][1]
        
        if 'targ_mag' in photdict_lcurve[mosaics]: 
            valu = float(photdict_lcurve[mosaics]['targ_mag'][0])
            verr = float(photdict_lcurve[mosaics]['targ_mag'][1])

            #there's probably a prettier way to do this, the second if statements 
            #are there so that only 1 label per filter is on the legend

            if 'j_' in mosaics:            
                if j == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue', label = 'J')
                    j = True

            elif 'h_' in mosaics:
                if h == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green', label = 'H')
                    h = True


            elif 'k_' in mosaics:
                if k == True: 
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
                else:
                    matplotlib.pyplot.errorbar(time, valu, yerr=verr, xerr=terr, \
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red', label = 'Ks')
                    k = True
        elif 'upper_green' in photdict_lcurve[mosaics]: 
                valu = float(photdict_lcurve[mosaics]['upper_green'])
                if photdict_lcurve[mosaics]['filter'] == 'h':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'green', mec = 'green', \
                        ecolor = 'green')
                if photdict_lcurve[mosaics]['filter'] == 'j':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')
                if photdict_lcurve[mosaics]['filter'] == 'k':
                    matplotlib.pyplot.errorbar(time, valu, xerr=terr, \
                        marker = 'v', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')
        else:
            print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (mosaics)
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
#    matplotlib.pyplot.xlabel('Time since Trigger (s)')
    matplotlib.pyplot.ylabel('Mag', fontsize=26)
    #matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    matplotlib.pyplot.legend()
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    uniquename = photdict_lcurve.keys()[0].split('_')[2]
#    matplotlib.pyplot.title(uniquename+' Lightcurve: ap = ' + str(ap))

    # spectral index plot
#    beta_err = the_fit[1][1] #why take square-root again? This looks wrong
    beta_err = numpy.std(betalist)
    ax = fig.add_axes([0.1,0.1,0.8,0.3])
    print 'betalist is'
    print betalist
    timlist = list(timlist)
    terlist = list(terlist)
    for time in notime_list:
        del timlist[time]
        del terlist[time]
    timlist = numpy.array(timlist)
    terlist = numpy.array(terlist)
    print 'timlist is'
    print timlist
    ax.errorbar(timlist, betalist, yerr = beta_err, xerr=terlist,\
                        marker = 'o', linestyle ='None', mfc = 'red', mec = 'green', \
                        ecolor = 'red')

    #ax = matplotlib.pyplot.gca()
    #ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    
    matplotlib.pyplot.xlabel('Time since Trigger (s)', fontsize=19)
    #matplotlib.pyplot.ylabel('$\beta$', fontsize=26)
    ylabel(r'$\beta$', fontsize=24)
    matplotlib.pyplot.semilogx()
    ax = matplotlib.pyplot.gca()
    #ax.set_xscale('log')
    matplotlib.pyplot.legend()
    uniquename = photdict.keys()[0].split('_')[2]
    cname = 'beta'
    #matplotlib.pyplot.title(uniquename+' Color Evolution '+'('+cname+')')
    savepath = storepath + uniquename + '_' + cname + '_colorevo.eps'
    print 'savepath'
    print savepath
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
        
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    DPI = F.get_dpi()
    if big:
        F.set_size_inches((20, 15))
    
    print 'colorevolution saved to ' + savepath
    savefig(savepath)    
    matplotlib.pyplot.close()

def test_fit():
    from MiscBin import q
    k_mag = 13.3697
    j_mag = 15.1144
    k_flux = -2.5*numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) - 48.60 # zeropoints are in erg/s/cm^2/Hz, convert to AB mag
    j_flux = -2.5*numpy.log10(q.mag2flux(j_mag,0,1.594e-20)) - 48.60

    k_length = -2.5*numpy.log10(2.398*10**14) - 48.60
    j_length = -2.5*numpy.log10(1.394*10**14) - 48.60

    result = (k_flux - j_flux)/(k_length - j_length)
    
    return result

def update_z(photdict, z, z_err):
    '''Updates the redshifts of photdicts'''
    z_dict = {'z': (z,z_err)}
    photdict.update(z_dict)

def update_all_z(photloop_result, z, z_err):
    '''Updates all redshifts in a photloop dictionary result'''
    for epoch in photloop_result:
        update_z(photloop_result[epoch] , z, z_err)
