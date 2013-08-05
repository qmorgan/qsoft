"""
Migration.py
Author: Adam Morgan
Created: July 24 2013
	
"""
import os
import glob
import numpy as np
from os import system

sextractor_bin = "sex"
loadpath = './'

def read_xml_configuration(configfile):
    import xml.etree.ElementTree as ET
    tree = ET.parse(configfile)
    root = tree.getroot()
    
    for item in root.iter('PVScan'):
        scan_date = item.attrib['date']
    
    #format is:
    # <Key key="positionCurrent_ZAxis" permissions="WriteSave" value="-518.6" />
    z_values = []
    x_values = []
    microns_per_pixel_list = []
    y_values = []
    z_diff_values = []
    # parse xml file and build list of values
    for key in root.iter('Key'):
        if key.attrib['key'] == "positionCurrent_ZAxis":
            z_values.append(float(key.attrib['value']))
        elif key.attrib['key'] == "positionCurrent_XAxis":
            x_values.append(float(key.attrib['value']))
        elif key.attrib['key'] == "positionCurrent_YAxis":
            y_values.append(float(key.attrib['value']))
        elif key.attrib['key'] == "micronsPerPixel_XAxis":
            microns_per_pixel_list.append(float(key.attrib['value']))
        elif key.attrib['key'] == "micronsPerPixel_YAxis":
            microns_per_pixel_list.append(float(key.attrib['value']))
    # determine difference between each z position
    for i in np.arange(len(z_values)-1):
        z_diff_values.append(z_values[i]-z_values[i+1])
    z_diff_mean = np.mean(z_diff_values)
    microns_per_pixel = np.mean(microns_per_pixel_list)
    # Error checking
    if np.std(z_diff_values) > 1E-9:
        raise ValueError("Something is wrong; Variance in z_diff detected")
    if np.std(x_values) > 1E-9:
        raise ValueError("Something is wrong; Variance in x_values detected")
    if np.std(y_values) > 1E-9:
        raise ValueError("Something is wrong; Variance in y_values detected")    
    if np.std(microns_per_pixel_list) > 1E-9:
        raise ValueError("Something is wrong; Variance in microns per pixel detected")
        
    config = {'date':scan_date,'z_diff':z_diff_mean,'pixel_scale':microns_per_pixel}
    return config
    
def tiff2fits(filelist,truncate_filename=True):
    '''Convert all tiff files in the current directory to a fits file
    truncate_filename renames the file as the last characters after the last 
    '_' in the input filename.
    '''
    working_directory = os.getcwd()
    outlist = []
    for filepath in filelist:
        assert os.path.exists(filepath)
        filedir = os.path.dirname(filepath)
        filebase = os.path.basename(filepath)
        basename, baseextension = filebase.split('.')
        if truncate_filename:
            basename = basename.split('_')[-1]
        if baseextension[0:3] != 'tif':
            raise ValueError('Expected tiff file; got ' + baseextension)
        outpath = working_directory + '/' + basename + '.fits'
        
        magickcommand = "convert %s %s" % (filepath, outpath)
        try:
            os.system(magickcommand)
            outlist += outpath
        except:
            print "Do not have Imagemagick, cannot convert svg to png"
            
    return outlist


def make_sex_cat(image_name, checkimages=False):
    '''Create the Source EXtraction CATalog'''
    basename = image_name.split('.')[0]
    sexcat_file = 'cat'+basename+'.txt'
    sex_command = sextractor_bin + " " + image_name + " -c " + loadpath + "cell.sex " + \
        " -CATALOG_NAME " + sexcat_file
    if checkimages:
        sex_command += ' -CHECKIMAGE_TYPE APERTURES,SEGMENTATION,BACKGROUND '
        sex_command += ' -CHECKIMAGE_NAME ap_%s.fits,seg_%s.fits,back_%s.fits ' % (basename,basename,basename)
    print sex_command
    system(sex_command)
    return sexcat_file

def make_ass_cat(image_name, assoc_image_name):
    '''Create the ASSociation CATalaog'''
    basename = image_name.split('.')[0]
    assoc_basename = assoc_image_name.split('.')[0]
    assoc_catname = assoc_basename + '.txt'
    if not os.path.exists(assoc_catname):
        # errmsg = 'SEx catalog for %s does not exist; cannot run ASSOC on it' % (assoc_catname)
        # raise ValueError(errmsg)
        make_sex_cat(assoc_image_name)
    sexcat_file = basename+'_'+assoc_basename+'_assoc.txt'
    sex_command = sextractor_bin + " " + image_name + " -c " + loadpath + "cell_assoc.sex " + \
        " -CATALOG_NAME " + sexcat_file + " ASSOC_NAME " + assoc_catname 
    system(sex_command)
    return sexcat_file


def sex_loop(checkimages=False):

    # give previous_image as None if first in the list; give following_image as None if last
    fits_list = glob.glob('0*.fits')
    for image in fits_list:
        make_sex_cat(image,checkimages=checkimages)
    
        