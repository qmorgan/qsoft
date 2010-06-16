"""
This program is the bee's knees. I'll add more explanatory text and usage 
examples to this space as the program gets closer to completion.

It is super important that you have the files
    pairitel_photo.sex
    pairitel_photo.param
    default.nnw
    star_find.param
in the same directory as this program file.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Author: Christopher Klein
Contact: cklein@berkeley.edu
Date: began 2010/02/11
"""


# --------------------------    LIBRARY and MODULE IMPORTS  --------------------

import sys
import ephem
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
from MiscBin import t_mid
from MiscBin import qPickle
from pylab import close
from time import time

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
sextractor_bin = "/opt/local/bin/sex"

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
    fwhm = mean(fwhm_list)
    fwhm_stdev = std(fwhm_list)
    return fwhm, fwhm_stdev
    

# --------------------------    BEGIN PROGRAM   --------------------------------

def dophot(progenitor_image_name,region_file, ap=None, find_fwhm = False, do_upper = False):
    # Begin timing
    t1 = time()
     
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
        "-CHECKIMAGE_NAME check.fits " + \
        "-PHOT_APERTURES 5.8"
    
    system(syscmd)


    objects_image = "check.fits"
    hdulist = pyfits.open(objects_image)
    objects_data = hdulist[0].data
    hdulist.close()
    system("rm check.fits")

    sex_file = file(storepath + "star_cat.txt", "r")
    sat_locations = []
    for line in sex_file:
        if (line[0] == "#"):
            continue
        else:
            if float(line.split()[2]) == 0:
                x_index = int(round(float(line.split()[0]) - 1))
                y_index = int(round(float(line.split()[1]) - 1))
                sat_locations.append([y_index, x_index])
    
    if find_fwhm or not ap:
        if not ap:
            print 'Aperture is not specified; fitting FWHM to determine'
        fwhm, fwhm_stdev = fit_fwhm(sat_locations, objects_data, 2.8, 0.8)
        print "Initial fwhm fit complete."
        print ("FWHM mean: " + str(fwhm) + " stdev: " + str(fwhm_stdev) + ".")
        fwhm, fwhm_stdev = fit_fwhm(sat_locations, objects_data, fwhm, fwhm_stdev)
        print "Sigma clipping of fwhm complete."
        print ("FWHM mean: " + str(fwhm) + " stdev: " + str(fwhm_stdev) + ".")
        photdict.update({'FHWM':(fwhm,fwhm_stdev)})        
    
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

    # Read in the image and convert data to array. Also store STRT_CPU, STOP_CPU, 
    # and FILTER header keyword values for future use. If the start and stop times
    # are not in the header, use placeholder values.
    imagefile_hdulist = pyfits.open(progenitor_image_name)
    imagefile_header = imagefile_hdulist[0].header
    image_data = imagefile_hdulist[0].data
    try:
        strt_cpu = str(imagefile_header["STRT_CPU"])
        stop_cpu = str(imagefile_header["STOP_CPU"])
    except:
        strt_cpu = 999
        stop_cpu = 999
    band = str(imagefile_header["FILTER"])
    imagefile_hdulist.close()
    
    photdict.update({'STRT_CPU':strt_cpu})
    photdict.update({'STOP_CPU':stop_cpu})
    photdict.update({'filter':band})
    
    if do_upper:
        # Store the image's dimensions.
        height, width = image_data.shape[0], image_data.shape[1]
        # We will drop apertures of 3 pixel radius everywhere on the image where there
        # exists science data. Since apertures are circular there will be gaps between
        # adjacent apertures. Science data is determined as any pixel with a weightmap
        # value > 0. To make sure we drop all possible apertures we need to create an 
        # expanded weight mask. This is a mask array with value 1 where aperatures can 
        # be centered and 0 elsewhere. It is formed by shrinking slightly the science
        # data area from the original image (by about 4 pixels).
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
        # To run SourceExtractor in two-image mode we need to first generate the fake
        # stars image which dictates where the apertures will be dropped on the science
        # image. We make fake stars at evenly spaced (6 pixels) intervals, blur them a 
        # bit to make them look like stars, and then use the expanded weight mask to 
        # enforce the boundaries.
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
                        fake_star_data[int((j+1-1./(order*2))*spacing+aperture_size*h)][int((i+1-1./(order*2))*spacing+aperture_size*g)] = peak
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

        print len(aperture_data), len(weight_aperture_data)
    

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
        # Create a histogram from the flux_array. Convert the flux_histogram_data into a
        # list for use with the gaussian fitting. Convert the flux_histogram_bins into a
        # list of length one less with values replaced with the average value of each 
        # bin. This is necessary because the output bins array is actually the bin
        # boundaries.
        # Repeat this code block 4 times, once for each color/quality subset.
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
    # Define the obj_string and band_type from the progenitor_image_name. If the 
    # progenitor_image_name is not in the format used by PARITIEL reduction
    # pipeline 3, these will not work properly.
    obj_string = progenitor_image_name.split("_")[2]
    band_type = (progenitor_image_name.split("_")[0] + "_" + 
        progenitor_image_name.split("_")[1])
    # Run source extractor on the science image.
    sexcat = make_sex_cat(progenitor_image_name, weight_image_name, aperture_size)
    # Store placehold values for the target photometry.
    target_mag = 999
    target_e_mag = 999
    # Create the catalog of 2MASS stars from the vizier data which corrresponds to
    # the filter/band used for the science image.
    vizcat_file = file(storepath + "viz_output_cropped.txt", "r")
    vizcat_starlist = []
    for line in vizcat_file:
        data_list = line.rstrip().lstrip().split(";")
        if len(data_list) == 13: # Check formatting 
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
    vizcat_file.close()
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
    # Compare the entries in sexcat_starlist and vizcat_starlist to create a 
    # combined_starlist which has entries for sources with both archival and new
    # instrumental data. The target need not be included in the archive.
    combined_starlist = []
    for sexcat_star in sexcat_starlist:
        sexcat_star_ra_rad = sexcat_star[0] * 0.01745329252
        sexcat_star_dec_rad = sexcat_star[1] * 0.01745329252
        target_separation_arcsec = 206264.806247*(float(ephem.separation(
            (target_ra_rad, target_dec_rad), 
            (sexcat_star_ra_rad, sexcat_star_dec_rad))))
        if target_separation_arcsec < 5:
            combined_starlist.append([sexcat_star[0], sexcat_star[1], 
                999, 999, 
                sexcat_star[2], sexcat_star[3], 
                sexcat_star[4]])
            continue
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
                    sexcat_star[4]])
    # Use the combined_starlist to calculate a zeropoint for the science image.
    zeropoint_list = []
    zeropoint_err_list = []
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
        if ((target_separation_arcsec > 5) and 
            ((ptel_flag == 0) or (ptel_flag == 2))):
            zeropoint_list.append(tmass_mag - ptel_mag)
            zeropoint_err_list.append(sqrt(tmass_e_mag*tmass_e_mag + 
                ptel_e_mag*ptel_e_mag))
    zeropoint = average(zeropoint_list)
    zeropoint_error = average(zeropoint_err_list)
    # Now apply the zeropoint to the instrumental magnitudes and create the 
    # final_starlist. Store the target photometry in target_mag and target_e_mag.
    final_starlist = []
    for star in combined_starlist:
        # If the star is our target . . .
        if (206264.806247*(float(ephem.separation(
            (target_ra_rad, target_dec_rad), 
            (star[0] * 0.01745329252, star[1] * 0.01745329252)))) < 5 and
            ((star[6] == 0) or (star[6] == 2) or (star[6] == 4) or (star[6] == 6) or 
            (star[6] == 7))):
            ra = star[0]
            dec = star[1]
            tmass_mag = star[2]
            tmass_e_mag = star[3]
            ptel_mag = star[4]
            ptel_e_mag = star[5]
            ptel_flag = star[6]
            new_mag = ptel_mag + zeropoint
            new_e_mag = sqrt(zeropoint_error*zeropoint_error + 
                ptel_e_mag*ptel_e_mag)
            target_mag = new_mag
            target_e_mag = new_e_mag
            final_starlist.append([ra, dec, tmass_mag, tmass_e_mag, 
                ptel_mag, ptel_e_mag, ptel_flag, new_mag, new_e_mag])
            continue
        # If the star is just a field 2MASS star . . .
        if star[6] == 0:
            ra = star[0]
            dec = star[1]
            tmass_mag = star[2]
            tmass_e_mag = star[3]
            ptel_mag = star[4]
            ptel_e_mag = star[5]
            ptel_flag = star[6]
            new_mag = ptel_mag + zeropoint
            new_e_mag = sqrt(zeropoint_error*zeropoint_error + 
                ptel_e_mag*ptel_e_mag)
            final_starlist.append([ra, dec, tmass_mag, tmass_e_mag, 
                ptel_mag, ptel_e_mag, ptel_flag, new_mag, new_e_mag])
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
        photdict.update({'targ_mag':(target_mag,target_e_mag)})
    # Print out photometry data. 
    print progenitor_image_name, "HJD:", hjd
    print "Photometry Results:", photometry_string
    print ("2MASS Catalog comparison average absolute deviation: " + 
        str(average(abs_2mass_deviation_list)))
    print "Zeropoint:", zeropoint, "err", zeropoint_error
    
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
    
    faintest_s2n = faintest_flux/faintest_flux_err
    print ("SExtractor faintest flux: " + str(faintest_flux) + 
        " err " + str(faintest_flux_err) + " => S/N = " + str(faintest_s2n))
    
    photdict.update({'faintest_s2n':faintest_s2n})
    
    photdict.update({'sex_faintest':(sex_faintest,sex_faintest_err)})
    
    # Clean up the photometry catalog files.
    system("rm viz_output.txt")
    system("rm viz_output_cropped.txt")
    system("rm " + progenitor_image_name.replace(".fits", ".sex"))
    system("rm " + weight_image_name.replace(".fits", ".sex"))
    
    print ("Photometry completed, " + 
        str(time() - t1) + " seconds required.")
        
    return photdict
    
def do_dir_phot(photdir='./',reg='PTEL.reg',ap=None):
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
        if 'targ_mag' not in photout and 'upper_green' not in photout:
            print '**Target Magnitude not found. Re-running to find UL**.'
            photout = dophot(myfile,reg,ap,do_upper=True)
        label = photout['FileName']
        photdict.update({label:photout})
    return photdict

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

def photloop(filename, reg, aper=None):
    '''Loops through all the mosaics on h/j/k_mosaics.txt from CoaddWrap, running q_super_photometry on each. Output is stored in a text file with the format t_mid|t_miderror|magnitude|magnitudeerror'''
    close()
    f = file(filename)
    textname = filename[0:1] + '_photometry_results.txt'
    r = open(textname, 'w')
    #vallist = []
    #errlist = []
    timlist = []

    for lameline in f.readlines():
        line = lameline.rstrip('\n')
        print line
        magnitude = dophot(line,reg,ap=aper)
        if magnitude['targ_mag'][0] == 'nan':
            pass
        else:
            vallist = magnitude['targ_mag'][0]
            errlist = magnitude['targ_mag'][1]
        time = t_mid.t_mid(line)
        time_err = t_mid.t_mid(line, delta = True) 

        result = str(time) + '|' + str(time_err) + '|' + str(vallist) + '|' + str(errlist) + '\n'
        r.write(result)
        
    f.close()
    r.close()
    
def photplot(GRBname):
    '''Plots a graph from the pickle output of the photreturn function'''
    import matplotlib
    import glob

    filepath = GRBname + '.data'
    photdict = qPickle.load(filepath)
    
    h = False
    j = False
    k = False

    timlist = []
    terlist = []
    vallist = []
    errlist = []

    for mosaics in photdict:
        time = float(t_mid.t_mid(mosaics))
        terr = float(t_mid.t_mid(mosaics, delta = True))/2.
        valu = float(photdict[mosaics]['targ_mag'][0])
        verr = float(photdict[mosaics]['targ_mag'][1])

#there's probably a prettier way to do this, the second if statements are there so that only 1 label per filter is on the legend

        if 'h_' in mosaics:
            if h == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'red', mec = 'green', ecolor = 'red')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'red', mec = 'green', ecolor = 'red', label = 'h')
                h = True

        elif 'j_' in mosaics:            
            if j == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'blue', mec = 'green', ecolor = 'blue')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'blue', mec = 'green', ecolor = 'blue', label = 'j')
                j = True

        elif 'k_' in mosaics:
            if k == True: 
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'yellow', mec = 'green', ecolor = 'yellow')
            else:
                matplotlib.pyplot.errorbar(time, valu, yerr = verr, xerr= terr, marker = 's', linestyle ='None', mfc = 'yellow', mec = 'green', ecolor = 'yellow', label = 'k')
                k = True

    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    matplotlib.pyplot.semilogx()
    matplotlib.pyplot.legend()
    
    savefig('lightcurve')    
    matplotlib.pyplot.close()


def photreturn(GRBname, filename, Clobber=False, reg=None, aper=None):
    
    filepath = GRBname + '.data'
    while Clobber == False:
        if os.path.isfile(filepath) == True:
            data = qPickle.load(filepath)
            if filename in data:
                return data[filename]
            else:
                Clobber = True
        else:
            Clobber = True

    while Clobber == True:
        if reg == None: 
            print 'Need to input reg file'
            break
        else:
            if os.path.isfile(filepath) == False:
                photdict = {}  
            else:
                #f = file(filepath)
                photdict = qPickle.load(filepath)
            data = dophot(filename, reg, aper)
            label = data['FileName']
            photdict.update({label:data})
            #f = file(filepath, 'w')
            qPickle.save(photdict, filepath, clobber = True)
            #.close()
            return data
            break


