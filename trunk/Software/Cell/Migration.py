"""
Migration.py
Author: Adam Morgan
Created: July 24 2013
	
"""
import os
import glob
import numpy as np
from os import system
import pandas as pd


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
def _overlay_verification_text(obj_str,xx,yy,theta,color='red',double_size=True):    
    
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
    
    im_cmd_append = ' -fill %s -pointsize 20 -annotate %fx%fXSIGN%fYSIGN%f "%s" ' % (color,ang,ang,abs(xpos),abs(ypos),obj_str)
    im_cmd_append = im_cmd_append.replace('XSIGN',xsign)
    im_cmd_append = im_cmd_append.replace('YSIGN',ysign)
    return im_cmd_append
    
    # get all the detections associated with object 60
    # det.loc[link.loc[link.obj_id==60].det_id]
    
    
    #convert testlog.png -fill red -annotate %fx%f%s%f%s%f "%i" testannotate.png
    # % (ang,ang,xsign,xpos,ysign,ypos,det_index)
    
def _create_verification_images(detection_table,object_table,link_table):
     
    #convert ap_000017.fits -fx "log(abs((u.r+u.g+u.b)*80))" testlog.png
    ap_list = glob.glob('ap_0*.png')
    
    for ap_img in ap_list:
        num = int(ap_img.replace('ap_','').lstrip('0').replace('.png',''))
        out_img = ap_img.replace('ap_','detections_')
        
        # get the detection table from the 5th image slice
        current_img_detections = detection_table[detection_table.z_index==num]
        # join it with the detection_object link table
        current_link = link_table.join(current_img_detections.z_index,on='det_id',how='right')
        
        # join it with the object stats table
        current_objects=current_link.join(object_table,on='obj_id',how='left')
        # loop through the objects and make a label for each
        im_cmd = 'convert ' + ap_img + ' -scale 200% -gravity Center ' 
        ### here we're centering on each DETECTION and printing a list of each OBJECT that detection is associated with
        for index,row in current_img_detections.iterrows():
            matching_objects = current_objects[current_objects['det_id'] == index]
            obj_list_str = str(list(matching_objects.obj_id)).rstrip(']').lstrip('[').replace('.0','')
            if obj_list_str == 'nan': 
                color = 'red'
            else:
                color = 'green'
            im_cmd+=_overlay_verification_text(obj_list_str,row['x'],row['y'],row['theta'],color=color)
        ### here we're centering on the average position of each OBJECT and printing it at this average position
        # for index, row in current_objects.iterrows():
        #     im_cmd+=_overlay_verification_text(str(row['obj_id']),row['x_mean'],row['y_mean'],0)
        im_cmd += out_img 
        os.system(im_cmd)
    print im_cmd
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

def build_image_table(configfile):
    '''
    Build up the image table from the configuration file.
    '''
    # df=pd.read_csv('cat000010.txt',delim_whitespace=True)
    configfile = 'ZSeries-07192013-1148-028.xml'
    print 'Reading config File %s' % configfile
    config = read_xml_configuration(configfile)
    image_table = pd.DataFrame(config['data_dict'])
    return image_table

def find_all_objects():
    '''    
    Build up the full detection table by reading in all the ascii outputs from 
    source extractor.  
    
    currently must run within the directory containing the cat0*.txt files output 
    by sextractor
    
    '''

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

def ignore_this_code():
        pass 
        #### ALGORITHM 2
        ## Create a database of unique objects
        # Is there an object at detected position in the potential object database ?
        #   NO: Create one
        #   YES: is there a detection for this object in the same frame as this one?
        #       YES: Create a new object for this detection
        #       NO: is the detection more than 20(?) frames away from the object already in the database?
        #           YES: Create a new object for this detection
        #           NO: it is probably associated with an object in the database already
        
        # loop through the unique objects and find possible associations of detections
        
        # if there are multiple associations with a detection in another slice, choose the nearest one?
    
        # ## ALGORITHM 1
        # testtab =det.loc[link.loc[link.obj_id==107].det_id]
        # qq = testtab
        # tt = testtab
        # qq['self_index']=qq.index # assigning as a new column
        # tt['self_index']=tt.index
        # qq2=qq.loc[:,['self_index','z_index','x','y']]
        # tt2=tt.loc[:,['self_index','z_index','x','y']]
        # new=pd.merge(tt2,qq2,on='z_index',suffixes=('_1','_2'))
        # duplicates = new[new.self_index_1 < new.self_index_2]
        # dup_ind = duplicates.z_index
        # 
        # # list of the ones where there were only SINGLE detections - use a logical not
        # no_neighbors = qq[np.logical_not(qq.z_index.isin(list(dup_ind)))]
        # # list of only the ones where there were multiple detections in the same frame:
        # neighbors = qq[qq.z_index.isin(list(dup_ind))]
        # # so far i can only deal with double neighbors, so im going to put this here: [could possibly loop thru this process until no neighbors left]
        # assert len(list(neighbors.z_index))==len(set(neighbors.z_index))
        # # the following however, should ALWAYS be true, unless i screwed up
        # assert len(list(no_neighbors.z_index))==len(set(no_neighbors.z_index))
        # 
        # # for all the doubles found, find the pairs that have the *nearest* positions
        # # 
        # # In [667]: duplicates
        # # Out[667]:
        # #     self_index_1  z_index  self_index_2
        # # 4            896       46           897
        # # 8            917       47           918
        # # 12           940       48           941
        # # 16           959       49           960
        # # 20           976       50           977
        # # 24           994       51           995
        # # 28          1010       52          1011
        # # 32          1024       53          1025
        # # 36          1037       54          1038
        # # take the two lists of duplicates self_index 1 and self_index 2
        # # take the first detection and compare it slice by slice to see which
        # # of the two possible sources it is nearest.  take the nearest one and store
        # # it in object 1, and the other should go to object 2. Repeat for the second duplicate
        # # and verify you get the opposite answer. 
        
        
        # IMPLEMENT THIS    
        #must have three unique NONBLENDED detections.. but still count the blended ones as once
        
def get_dat_ass_old(full_det_table,xy_rad=5.0,slice_rad=30.0):
    '''
    DEPRECIATED
    '''
            
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
    for index,row in full_det_table.iterrows():
        current_x = row['x']
        current_y = row['y']
        current_z = row['z_index']
        
        
        #######
        # STEP 1: Gather table of nearby detections within distance and slice threshold
        # full_det_table[np.sqrt((current_x-full_det_table.x)**2+(current_y-full_det_table.y)**2) < xy_rad]
        # Another common operation is the use of boolean vectors to filter the data. 
        # The operators are: | for or, & for and, and ~ for not. These must be 
        # grouped by using parentheses.
        
        # Here is the logic for determining table of nearby detections
        nearby_detections_boolean = (np.sqrt((current_x-full_det_table.x)**2+(current_y-full_det_table.y)**2) \
            < xy_rad) & (np.abs(current_z - full_det_table.z_index) < slice_rad) \
            & (current_z != full_det_table.z_index) # added this to avoid double-counting of nearby cells in same frame 
        
        # get integer index of label index [just in case the lengths are different]
        current_int_index = np.where(full_det_table.index == index)[0][0]
        
        # set boolean value to True for current index position so it is included
        nearby_detections_boolean[current_int_index] = True
        
        # extract the subset of the full detections 
        nearby_detections_table = full_det_table[nearby_detections_boolean]  
        
        # Don't add as an object or continue down the algorithm if there are
        # fewer than 3 detections. Probably not a real object.  
        if not len(nearby_detections_table.index) >= 3: # only include if there are at least 3 frames within the search radius at this position
            print "Too few detections for index %s" % (str(index))
            print nearby_detections_table.index
            continue
        
        
        # check if all associated items have already been accounted for by seeing if all detections
        # in the object frame are already in the detections_accounted_for_set
        if not set(nearby_detections_table.index).issubset(detections_accounted_for_set):
            detections_accounted_for_set.update(set(nearby_detections_table.index))
            
            # build up list of potential objects. CHANGE THIS.
            pot_obj_list.append(nearby_detections_table)
        
        # IMPLEMENT THIS    
        #must have three unique NONBLENDED detections.. but still count the blended ones as once
        
        
        ## ALGORITHM 3
        from collections import Counter
        from itertools import groupby
        
        # create a list of the modes of slices (if there are multiple) by counting the frequency of each slice in the detection list
        slice_freq = groupby(Counter())
        

        # get data frame at specific location (using labels)
        #.loc is strictly label based, will raise KeyError when the items are not found
        # full_det_table.loc[[1091,1092,1093]] # returns values based on this index
        # full_det_table data frame at specific location (using integer position)
        # .iloc is strictly integer position based (from 0 to length-1 of the axis), 
        # will raise IndexError when the requested indicies are out of bounds.
        # full_det_table.loc[[1091,1092,1093]] # returns values based on this index
        
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
        
        # attempting to make a multiindex table
        # obj_id_arr = np.zeros(len(obj_df.index))
        # det_id_arr = np.array(obj_df.index)
        # arrays = []
        # nevermind for now
        
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
    
    print str(len(object_table)) + ' potential objects found'
    
    return object_table,detection_object_link_table 
 
 
def get_dat_ass(full_det_table,xy_rad=5.0,slice_rad=30.0):
    '''
    Get database of associations
    
    Parameters
    ----------
    full_det_table : pandas detection table of all detection across all slices
        generated by find_all_objects()
    xy_rad : distance search radius in the x-y plane in number of pixels. This
        should be relatively large and able to capture nearby sources in the 
        same z-slice as well as the same source across z-slices
    slice_rad : slice search radius in the z-plane. Single sources can persist
        across 10+ slices, but should not be many slices away.
        
        TODO: add a check at the object assignment stage to make sure there
        are no gaps of few+ slices between any two detections
    
    '''
            
    detections_accounted_for_set = set()
    set_of_detections_accounted_for = set()
    
    pot_obj_list = []
    
    z_mean_list = []
    z_std_list = []
    x_mean_list = []
    x_std_list = []
    y_mean_list = []
    y_std_list = []
    det_list_list = []
    
    object_table = None
    object_names = ['z_mean','z_std','x_mean','x_std','y_mean','y_std']
    object_count = 0
    
    # get database of associations
    for index,row in full_det_table.iterrows():
        current_x = row['x']
        current_y = row['y']
        current_z = row['z_index']
        
        
        #######
        # STEP 1: Gather table of nearby detections within distance and slice threshold
        # full_det_table[np.sqrt((current_x-full_det_table.x)**2+(current_y-full_det_table.y)**2) < xy_rad]
        # Another common operation is the use of boolean vectors to filter the data. 
        # The operators are: | for or, & for and, and ~ for not. These must be 
        # grouped by using parentheses.
        
        # Here is the logic for determining table of nearby detections
        nearby_detections_boolean = (np.sqrt((current_x-full_det_table.x)**2+(current_y-full_det_table.y)**2) \
            < xy_rad) & (np.abs(current_z - full_det_table.z_index) < slice_rad) \
            & (current_z != full_det_table.z_index) # added this to avoid double-counting of nearby cells in same frame 
        
        # get integer index of label index [just in case the lengths are different]
        current_int_index = np.where(full_det_table.index == index)[0][0]
        
        # set boolean value to True for current index position so it is included
        nearby_detections_boolean[current_int_index] = True
        
        # extract the subset of the full detections 
        nearby_detections_table = full_det_table[nearby_detections_boolean]  
        
        # Don't add as an object or continue down the algorithm if there are
        # fewer than 3 detections. Probably not a real object.  
        if not len(nearby_detections_table.index) >= 3: # only include if there are at least 3 frames within the search radius at this position
            print "Too few detections for index %s" % (str(index))
            print nearby_detections_table.index
            continue
        
        
        ## STEP 2: Find the z-positions in the nearby detection table that have 
        # the most detections. 
        from collections import Counter
        from itertools import groupby
        
        
        ## STEP 3: Iterate until we find a number of slices with a common
        # number of detections. This fights against single outliers that 
        # only appear in a single image
        real_detection_threshold = 2
        most_common_slices = []
        while len(most_common_slices) < real_detection_threshold:
            
            # create a list of the modes of slices (if there are multiple) by 
            # counting the frequency of each slice in the detection list
            slice_counter = Counter(list(nearby_detections_table.z_index))
            slice_freq_iterator = groupby(slice_counter.most_common(), lambda x:x[1])
            most_common_slices = [freqval for freqval,freqcount in slice_freq_iterator.next()[1]]

            # get the subset table of detections have the most number of nearby 
            # detections in a single slice 
            most_det_in_ind_slice_table = nearby_detections_table[nearby_detections_table.z_index.isin(most_common_slices)]
            
            if len(most_det_in_ind_slice_table) == 1: 
                print "Down to the lowest level."
                break
            if len(most_common_slices) == 1:
                print "Just one instance of %s object detections in slice %s; moving to next level" % (str(len(most_det_in_ind_slice_table)),str(most_common_slices))
                nearby_detections_table = nearby_detections_table[np.logical_not(nearby_detections_table.z_index.isin(most_common_slices))]

   
            
        # if all the detections in each of these slices is accounted for already, 
        # then the associated object with them should exist
        # if SOME detections are accounted for already but not all, 
        # this is an indication that something is wrong and the threshold is too
        # low. 
        
        # need both of these?
        cont = False
        for inddd in most_det_in_ind_slice_table.index:
            if inddd in set_of_detections_accounted_for:
                cont = True
        if cont == True:
            continue 
            
        if not set(most_det_in_ind_slice_table.index).issubset(set_of_detections_accounted_for): 
            set_of_detections_accounted_for.update(set(most_det_in_ind_slice_table.index))
        else:
            for inddd in most_det_in_ind_slice_table.index:
                assert inddd not in set_of_detections_accounted_for
            continue
        
    

        
        # take the middle slice from this list as the reference point for 
        # object assignment. The middle slice is most likely to have the 
        # clearest detections of the objects and best positions
        slice_ref_index = len(most_common_slices)/2
        reference_slice = most_common_slices[slice_ref_index]
        
        reference_table = most_det_in_ind_slice_table[most_det_in_ind_slice_table.z_index == reference_slice]
        
        
        # gather the close by (in xy) detections by looping through the most_det_in_ind_slice_table 
        # and adding the closest detections to each object in the reference table
        for ref_ind, ref_row in reference_table.iterrows():
            # each reference detection will get assigned to an object.
            # create a table with just the reference info so far, and then add detections
            # have to do nearby_detections_table.z_index.loc[ref_ind] instead of say ref_row['z_index']
            # in order to preserve the data type
            check_table = pd.DataFrame([[0,nearby_detections_table.z_index.loc[ref_ind],ref_ind]],columns=['distance','z_index','det_id'])
            for sub_ind, sub_row in most_det_in_ind_slice_table.iterrows():
                dist = np.sqrt((ref_row['x']-sub_row['x'])**2+(ref_row['y']-sub_row['y'])**2)
                tmp_chk_table = pd.DataFrame([[dist,nearby_detections_table.z_index.loc[sub_ind],sub_ind]],columns=['distance','z_index','det_id'])
                # if no detection from that slice has been recorded yet, append it
                # if index == 41 and sub_ind == 86: raise Exception
                
                # if not sub_row['z_index'] in check_table.z_index: #THIS FAILS for some reason sometimes
                # completely boggles my mind. 
                # ipdb> check_table.z_index
                # 0    5
                # 1    4
                # Name: z_index, dtype: int64
                # ipdb> 4 in check_table.z_index
                # False
                # doing this workaround using boolean indexing instead
                
                if len(check_table[check_table['z_index'] == int(sub_row['z_index'])]) == 0:
                    check_table = check_table.append(tmp_chk_table, ignore_index=True)
                else:
                    #grab the index of the detection for the current slice
                    prev_index = check_table[check_table.z_index==sub_row['z_index']].index
                    # if the current distance is shorter, replace that row
                    if dist < check_table.loc[prev_index]['distance']:
                        check_table.distance[prev_index] = dist
                        check_table.det_id[prev_index] = sub_ind
            
            # if there's only one instance of this many objects in a slice, it's 
            # probably a spurious detection. skip and move down to the next most 
            # common detection level.
            if len(check_table.det_id) < 3:
                print "only %i instance(s) of %s object detections in slice %s at x=%f y=%f" % (len(check_table.det_id),str(len(reference_table)),str(int(ref_row['z_index'])),ref_row['x'],ref_row['y'])
                # could instead flag it as a spurious detection
                # or have a column in the object table of how many slices went into it
                
            
            # now that the detection ids of the closest objects have been found,
            # add these all to a table  
            
            obj_df = most_det_in_ind_slice_table.loc[check_table.det_id]
            means = obj_df.mean()
            stds = obj_df.std()
            object_names = ['z_mean','z_std','x_mean','x_std','y_mean','y_std','detections','associations']     
            ind_obj_list = [means['z_index'],stds['z_index'],means['x'],stds['x'],means['y'],stds['y'],len(check_table.det_id),len(check_table.det_id)]
            tmp_table = pd.DataFrame([ind_obj_list],columns=object_names)
            obj_det_link_dict = {'obj_id':object_count,'det_id':obj_df.index}
            tmp_table2 = pd.DataFrame(obj_det_link_dict)
            if object_count == 0:
                object_table = tmp_table
                detection_object_link_table = tmp_table2
                association_object_link_table = tmp_table2
            else:
                object_table = object_table.append(tmp_table,ignore_index=True)
                detection_object_link_table = detection_object_link_table.append(tmp_table2,ignore_index=True)
                association_object_link_table = association_object_link_table.append(tmp_table2,ignore_index=True)
            object_count += 1
            
            
            # if index > 1051: raise Exception
        
        # assign each detection in the reference table as a potential object

        
        # detection id, distance, slice
        
        
        
        # check if all associated items have already been accounted for by seeing if all detections
        # in the object frame are already in the detections_accounted_for_set
        if not set(nearby_detections_table.index).issubset(detections_accounted_for_set): 
            detections_accounted_for_set.update(set(nearby_detections_table.index))
            
            # build up list of potential objects. CHANGE THIS.
            pot_obj_list.append(nearby_detections_table)
        
        # IMPLEMENT THIS    
        #must have three unique NONBLENDED detections.. but still count the blended ones as once
        
            
    # # loop again through object list and populate the object_table and detection_object_link_table?
    # # may be a more efficient way to do this.. 
    # # 
    # object_names = ['z_mean','z_std','x_mean','x_std','y_mean','y_std']
    # 
    # count = 0
    # for obj_df in pot_obj_list:
    #     means = obj_df.mean()
    #     stds = obj_df.std()
    #     ind_obj_list = [means['z_index'],stds['z_index'],means['x'],stds['x'],means['y'],stds['y']]
    #     tmp_table = pd.DataFrame([ind_obj_list],columns=object_names)
    #     obj_det_link_dict = {'obj_id':count,'det_id':obj_df.index}
    #     tmp_table2 = pd.DataFrame(obj_det_link_dict)
    #     if count == 0:
    #         object_table = tmp_table
    #         detection_object_link_table = tmp_table2
    #     else:
    #         object_table = object_table.append(tmp_table,ignore_index=True)
    #         detection_object_link_table = detection_object_link_table.append(tmp_table2,ignore_index=True)
    #     
    #     # attempting to make a multiindex table
    #     # obj_id_arr = np.zeros(len(obj_df.index))
    #     # det_id_arr = np.array(obj_df.index)
    #     # arrays = []
    #     # nevermind for now
    #     
    #     count += 1
    #     
    #     catid = int(cat_text_file[3:].lstrip('0').rstrip('.txt'))
    #     tmp_table=pd.read_csv(cat_text_file,header=None,delim_whitespace=True,names=detection_names)
    #     tmp_table.z_index = catid
    #     if count == 0:
    #         detection_table = tmp_table # do this for the first item, then append on to it
    #     else:
    #         detection_table = detection_table.append(tmp_table,ignore_index=True)
    #     count += 1
    # print str(count) + ' tables read in'
    
    print str(len(object_table)) + ' potential objects found'
    
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
        
    xy_rad = 3.0 # association radius in number of pixels
    slice_rad = 10 # association radius in slices 
    
    current_x = 111.0
    current_y = 369.0
    current_z = 70
    
    detection_table[np.sqrt((current_x-detection_table.x)**2+(current_y-detection_table.y)**2) < xy_rad]
    # Another common operation is the use of boolean vectors to filter the data. 
    # The operators are: | for or, & for and, and ~ for not. These must be 
    # grouped by using parentheses.
    detection_table[(np.sqrt((current_x-detection_table.x)**2+(current_y-detection_table.y)**2) < xy_rad) & (np.abs(current_z - detection_table.z_index) < slice_rad)]
    
    # get data frame at specific location (using labels)
    #.loc is strictly label based, will raise KeyError when the items are not found
    detection_table.loc[[1091,1092,1093]] # returns values based on this index
    # detection_table data frame at specific location (using integer position)
    # .iloc is strictly integer position based (from 0 to length-1 of the axis), 
    # will raise IndexError when the requested indicies are out of bounds.
    # detection_table.loc[[1091,1092,1093]] # returns values based on this index
    
    
    