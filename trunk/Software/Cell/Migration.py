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
    rel_time_values = []
    abs_time_values = []
    indices = []
    # parse xml file and build list of values
    for frame in root.iter('Frame'):
        rel_time_values.append(float(frame.attrib['relativeTime']))
        abs_time_values.append(float(frame.attrib['absoluteTime']))
        indices.append(int(frame.attrib['index']))
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
        
    data_dict = {'z_index':indices,'z_pos':z_values,'x_pos':x_values,'y_pos':y_values,'t_rel':rel_time_values,'t_abs':abs_time_values}
    config = {'date':scan_date,'z_diff':z_diff_mean,'pixel_scale':microns_per_pixel,'data_dict':data_dict}

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
            print "Do not have Imagemagick, cannot convert tif to fits"
            
    return outlist


def _convert_verification_pngs():
    '''This takes several minutes, recommended only if you need to verify that
    your detections are being populated correctly 
    '''
    ap_list = glob.glob('ap_0*.fits')
    for ap_img in ap_list:
        png_name = ap_img.replace('fits','png')
        im_cmd = 'convert %s -fx "log(abs((u.r+u.g+u.b)*80))" %s' % (ap_img,png_name)
        os.system(im_cmd)
        # cp_cmd = 'cp %s %s' % (png_name,png_name.replace('ap_','detections_'))
        # os.system(cp_cmd)
def _overlay_verification_text(obj_str,xx,yy,theta,double_size=True):    
    
    if theta <= 0: 
        ang = abs(theta)
    else: 
        ang = 360. - theta
    
    if not double_size:
        xpos = xx - 256.  #subtracting 256 since changing gravity setting not only 
        # refers to the position on the background image, but also the position 
        # in the overlay image that is to be drawn.
        ypos = 512.- yy -256.
    else:
        xpos = xx*2 - 512.
        ypos = 1024. - yy*2 - 512
    
    if xpos < 0:
        xsign="-"
    else:
        xsign="+"
    if ypos < 0:
        ysign="-"
    else:
        ysign="+"
    
    im_cmd_append = ' -fill red -annotate %fx%fXSIGN%fYSIGN%f "%s" ' % (ang,ang,abs(xpos),abs(ypos),obj_str)
    im_cmd_append = im_cmd_append.replace('XSIGN',xsign)
    im_cmd_append = im_cmd_append.replace('YSIGN',ysign)
    return im_cmd_append
    
    #convert testlog.png -fill red -annotate %fx%f%s%f%s%f "%i" testannotate.png
    # % (ang,ang,xsign,xpos,ysign,ypos,det_index)
    
def _create_verification_images(detection_table,object_table,link_table):
     
    #convert ap_000017.fits -fx "log(abs((u.r+u.g+u.b)*80))" testlog.png
    
    # get the detection table from the 5th image slice
    current_img_detections = detection_table[detection_table.z_index==5]
    # join it with the detection_object link table
    current_link = link_table.join(current_img_detections.z_index,on='det_id',how='right')
        
    # join it with the object stats table
    current_objects=current_link.join(object_table,on='obj_id',how='left')
    # loop through the objects and make a label for each
    im_cmd = 'convert ap_000005.png -scale 200% -gravity Center'
    ### here we're centering on each DETECTION and printing a list of each OBJECT that detection is associated with
    for index,row in current_img_detections.iterrows():
        matching_objects = current_objects[current_objects['det_id'] == index]
        obj_list_str = str(list(matching_objects.obj_id)).rstrip(']').lstrip('[')
        im_cmd+=_overlay_verification_text(obj_list_str,row['x'],row['y'],row['theta'])
    ### here we're centering on the average position of each OBJECT and printing it at this average position
    # for index, row in current_objects.iterrows():
    #     im_cmd+=_overlay_verification_text(str(row['obj_id']),row['x_mean'],row['y_mean'],0)
    im_cmd += 'detections_000005.png'
    os.system(im_cmd)
    
    #convert testlog.png -fill red -annotate %fx%f%s%f%s%f "10" testannotate.png
    # % (ang,ang,xsign,xpos,ysign,ypos)
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
    if assoc_image_name == None:
        return None
    basename = image_name.split('.')[0]
    assoc_basename = assoc_image_name.split('.')[0]
    assoc_catname = assoc_basename + '.txt'
    if not os.path.exists(assoc_catname):
        # errmsg = 'SEx catalog for %s does not exist; cannot run ASSOC on it' % (assoc_catname)
        # raise ValueError(errmsg)
        make_sex_cat(assoc_image_name)
    sexcat_file = 'assoc_'+basename+'_'+assoc_basename+'.txt'
    sex_command = sextractor_bin + " " + image_name + " -c " + loadpath + "cell_assoc.sex " + \
        " -CATALOG_NAME " + sexcat_file + " -ASSOC_NAME " + assoc_catname 
    print sex_command
    system(sex_command)
    return sexcat_file


def sex_loop(checkimages=False):
    # give previous_image as None if first in the list; give following_image as None if last
    fits_list = glob.glob('0*.fits')
    for image in fits_list:
        make_sex_cat(image,checkimages=checkimages)
    assoc_list = [None] + fits_list + [None]
    # for index in np.arange(len(fits_list))+1:
    #     make_ass_cat(assoc_list[index], assoc_list[index-1])
    #     make_ass_cat(assoc_list[index], assoc_list[index+1])
    # 
    # assuming file format of cat000119.txt
    # catlist = glob.glob('cat0*.txt')
    # if len(catlist) == 0:
    #     raise IOError("List of catalogs is empty. Check path and verify presence of cat0xxxxx.txt files")
    # for cat in catlist:
    #     catid = int(item[3:].lstrip('0').rstrip('.txt'))
    #     

    # print catnumlist

def find_all_objects():
    import pandas as pd
    # df=pd.read_csv('cat000010.txt',delim_whitespace=True)
    configfile = 'ZSeries-07192013-1148-028.xml'
    print 'Reading config File %s' % configfile
    config = read_xml_configuration(configfile)
    image_table = pd.DataFrame(config['data_dict'])
    
    detection_names = ['z_index','x','y','A','B','ellipticity','theta','r_kron',
        'r_flux','mu_max','mag_auto','magerr_auto','mag_iso','magerr_iso',
        'fwhm','flags']
    
    count = 0
    cat_list = glob.glob('cat0*.txt')
    for cat_text_file in cat_list:
        catid = int(cat_text_file[3:].lstrip('0').rstrip('.txt'))
        tmp_table=pd.read_csv(cat_text_file,header=None,delim_whitespace=True,names=detection_names)
        tmp_table.z_index = catid
        if count == 0:
            detection_table = tmp_table # do this for the first item, then append on to it
        else:
            detection_table = detection_table.append(tmp_table,ignore_index=True)
        count += 1
    print str(count) + ' tables read in'
    return detection_table

def get_dat_ass(det,rad=3.0,slice_rad=10.0):
    import pandas as pd
    
    ass_rad = rad # association radius in number of pixels
    # slice_rad = 10 # association radius in slices
    
    detections_accounted_for_set = set()
    pot_obj_list = []
    
    z_mean_list = []
    z_std_list = []
    x_mean_list = []
    x_std_list = []
    y_mean_list = []
    y_std_list = []
    det_list_list = []
    
    # get database of associations
    for index,row in det.iterrows():
        current_x = row['x']
        current_y = row['y']
        current_z = row['z_index']

        # det[np.sqrt((current_x-det.x)**2+(current_y-det.y)**2) < ass_rad]
        # Another common operation is the use of boolean vectors to filter the data. 
        # The operators are: | for or, & for and, and ~ for not. These must be 
        # grouped by using parentheses.
        boolean_array = (np.sqrt((current_x-det.x)**2+(current_y-det.y)**2) \
            < ass_rad) & (np.abs(current_z - det.z_index) < slice_rad) \
            & (np.abs(current_z - det.z_index) != 0) # added this to avoid double-counting of nearby cells in same frame 
        
        # get integer index of label index [just in case the lengths are different]
        current_int_index = np.where(det.index == index)[0][0]
        
        # set boolean value to True for current index position so it is included
        boolean_array[current_int_index] = True
        
        potential_object_frame = det[boolean_array]  
        
        
        if not len(potential_object_frame.index) >= 3: # only include if there are at least 3 frames within the search radius at this position
            print "too few detections for index %s" % (str(index))
            print potential_object_frame.index
        # check if all associated items have already been accoutned for by seeing if all detections
        # in the object frame are already in the detections_accounted_for_set
        if not set(potential_object_frame.index).issubset(detections_accounted_for_set):
            detections_accounted_for_set.update(set(potential_object_frame.index))
            pot_obj_list.append(potential_object_frame)
            
        

        # get data frame at specific location (using labels)
        #.loc is strictly label based, will raise KeyError when the items are not found
        # det.loc[[1091,1092,1093]] # returns values based on this index
        # det data frame at specific location (using integer position)
        # .iloc is strictly integer position based (from 0 to length-1 of the axis), 
        # will raise IndexError when the requested indicies are out of bounds.
        # det.loc[[1091,1092,1093]] # returns values based on this index
        
    # loop again through object list and populate the object_table and detection_object_link_table?
    # may be a more efficient way to do this.. 
    # 
    object_names = ['z_mean','z_std','x_mean','x_std','y_mean','y_std']
    
    count = 0
    for obj_df in pot_obj_list:
        means = obj_df.mean()
        stds = obj_df.std()
        ind_obj_list = [means['z_index'],stds['z_index'],means['x'],stds['x'],means['y'],stds['y']]
        tmp_table = pd.DataFrame([ind_obj_list],columns=object_names)
        obj_det_link_dict = {'obj_id':count,'det_id':obj_df.index}
        tmp_table2 = pd.DataFrame(obj_det_link_dict)
        if count == 0:
            object_table = tmp_table
            detection_object_link_table = tmp_table2
        else:
            object_table = object_table.append(tmp_table,ignore_index=True)
            detection_object_link_table = detection_object_link_table.append(tmp_table2,ignore_index=True)
        count += 1
    
    #     catid = int(cat_text_file[3:].lstrip('0').rstrip('.txt'))
    #     tmp_table=pd.read_csv(cat_text_file,header=None,delim_whitespace=True,names=detection_names)
    #     tmp_table.z_index = catid
    #     if count == 0:
    #         detection_table = tmp_table # do this for the first item, then append on to it
    #     else:
    #         detection_table = detection_table.append(tmp_table,ignore_index=True)
    #     count += 1
    # print str(count) + ' tables read in'
    
    
    return object_table,detection_object_link_table 
        
def test_pandas():
    import pandas as pd
    # df=pd.read_csv('cat000010.txt',delim_whitespace=True)
    configfile = 'ZSeries-07192013-1148-028.xml'
    print 'Reading config File %s' % configfile
    config = read_xml_configuration()
    image_table = pd.DataFrame(config['data_dict'])
    
    detection_names = ['z_index','x','y','A','B','ellipticity','theta','r_kron',
        'r_flux','mu_max','mag_auto','magerr_auto','mag_iso','magerr_iso',
        'fwhm','flags']
    
    count = 0
    cat_list = glob.glob('cat0*.txt')
    for cat_text_file in cat_list:
        catid = int(cat_text_file[3:].lstrip('0').rstrip('.txt'))
        tmp_table=pd.read_csv(cat_text_file,header=None,delim_whitespace=True,names=detection_names)
        tmp_table.z_index = catid
        if count == 0:
            detection_table = tmp_table # do this for the first item, then append on to it
        else:
            detection_table = detection_table.append(tmp_table,ignore_index=True)
        count += 1
    print str(count) + ' tables read in'
    return detection_table
        
    ass_rad = 3.0 # association radius in number of pixels
    slice_rad = 10 # association radius in slices 
    
    current_x = 111.0
    current_y = 369.0
    current_z = 70
    
    det[np.sqrt((current_x-det.x)**2+(current_y-det.y)**2) < ass_rad]
    # Another common operation is the use of boolean vectors to filter the data. 
    # The operators are: | for or, & for and, and ~ for not. These must be 
    # grouped by using parentheses.
    det[(np.sqrt((current_x-det.x)**2+(current_y-det.y)**2) < ass_rad) & (np.abs(current_z - det.z_index) < slice_rad)]
    
    # get data frame at specific location (using labels)
    #.loc is strictly label based, will raise KeyError when the items are not found
    det.loc[[1091,1092,1093]] # returns values based on this index
    # det data frame at specific location (using integer position)
    # .iloc is strictly integer position based (from 0 to length-1 of the axis), 
    # will raise IndexError when the requested indicies are out of bounds.
    # det.loc[[1091,1092,1093]] # returns values based on this index
    