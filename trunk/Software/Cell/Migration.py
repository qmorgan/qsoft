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
import copy
import matplotlib.pyplot as plt
from matplotlib import rc
import cPickle as pickle

rc('font', family='Times New Roman', size=14)

sextractor_bin = "sex"


class ImageStack:
    def __init__(self, image_directory='./',output_directory='./',config_directory='./'):
        
        if image_directory[0:2]=='./' or image_directory[0:3]=='../':
            # assume the user is trying to run in current working directory
            os.chdir(image_directory)
            image_directory = os.getcwd() + '/'
        
        if output_directory[0:2]=='./' or output_directory[0:3]=='../':
            # assume the user is trying to run in current working directory
            os.chdir(output_directory)
            output_directory = os.getcwd() + '/'
        
        if config_directory[0:2]=='./' or config_directory[0:3]=='../':
            # assume the user is trying to run in current working directory
            os.chdir(config_directory)
            config_directory = os.getcwd() + '/'
            
        if ' ' in image_directory:
            raise ValueError("Spaces not allowed in image_directory: " + image_directory)
        if ' ' in output_directory:
            raise ValueError("Spaces not allowed in output_directory: " + output_directory)
        if ' ' in config_directory:
            raise ValueError("Spaces not allowed in config_directory: " + config_directory)

        self.config_directory = config_directory #path for configuration files cell.sex, cell.param, gauss_3.0_7x7.conv
        self.image_directory = image_directory # where the .tif files are located
        self.output_directory = output_directory # where all the output files should go
        
        xml_path = self.image_directory + 'ZSeries*.xml'
        xml_list = glob.glob(xml_path)
        if len(xml_list) > 1:
            raise IOError("Too many .xml configuration files in " + xml_path)
        elif len(xml_list) == 0:
            try:
                print "Couldn't find .xml configuration file; checking in parent directory"
                xml_path = os.path.dirname(self.image_directory.rstrip('/'))+'/ZSeries*.xml'
                xml_list = glob.glob(xml_path)
                assert len(xml_list) == 1
                self.xml_config = xml_list[0]
            except:
                raise IOError("Cannot find .xml configuration in path " + xml_path)
        else:
            self.xml_config = xml_list[0]
        self.name = os.path.basename(self.xml_config).replace('.xml','')
        self.read_xml_configuration()
    
    def PrepareImages(self):
        self.tiflist = glob.glob(self.image_directory+'ZSeries*.tif')
        if len(self.tiflist) == 0:
            print "Cannot find any ZSeries*.tif images to convert in " + self.image_directory
        print "Converting tiff files to fits files..."
        self.tiff2fits(truncate_filename=True)
        if len(self.fitslist) == 0:
            raise IOError("Couldn't convert tif files to fits files...")
        print "Done. Ready to FindDetections()"
    
    
    def FindDetections(self,checkimages=False):
        # give previous_image as None if first in the list; give following_image as None if last
        
        if not hasattr(self,"fitslist"):
            print "ImageStack instance does not have fitslist attribute."
            print "Either convert .tif data to fits using PrepareImages(), "
            print " or grab list of already generated .fits files in the "
            print " output_directory using grab_fits_list()."
            return
        self.catlist = []
        fits_list = self.fitslist
        
        for image in fits_list:
            self.catlist.append(self.make_sex_cat(image,checkimages=checkimages))
        print "Source Extraction complete. Wrote " + str(len(self.catlist)) + " catalogs."
        print "Ready to ReadCatalogs()"
        
    def Save(self,clobber=False):
        savePickle(self,self.output_directory+self.name+'.pkl',clobber=clobber)
    
    def ExportToExcel(self):
        excelpath = self.output_directory+self.name+'.xlsx'
        writer = pd.ExcelWriter(excelpath)
        self.object_table.to_excel(writer,sheet_name="Objects")
        self.detection_table.to_excel(writer,sheet_name="Detections")
        self.detection_object_link_table.to_excel(writer,sheet_name="ObjectLinks")
        self.association_object_link_table.to_excel(writer,sheet_name="AssociationLinks")
        writer.save()
        
    def read_xml_configuration(self):
        import xml.etree.ElementTree as ET
        tree = ET.parse(self.xml_config)
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
        config = {'filename':self.xml_config,'date':scan_date,'z_diff':z_diff_mean,'pixel_scale':microns_per_pixel,'data_dict':data_dict}
        
        self.config = config
        self.image_table = pd.DataFrame(config['data_dict'])
        
    def grab_fits_list(self,make_new=False):
        self.fitslist=glob.glob(self.output_directory+'0*.fits')
        if len(self.fitslist) == 0 and make_new == True:
            self.tiff2fits()
        
        
    def tiff2fits(self,truncate_filename=True):
        '''Convert all tiff files in the current directory to a fits file
        truncate_filename renames the file as the last characters after the last 
        '_' in the input filename.
        '''
        print 'Converting tiff files to fits files'
        working_directory = self.output_directory
        outlist = []
        for filepath in self.tiflist:
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
                outlist.append(outpath)
            except:
                print "Do not have Imagemagick, cannot convert tif to fits"
    
        self.fitslist = outlist
    
    def make_sex_cat(self,image_name, checkimages=False):
        '''Create the Source EXtraction CATalog'''
        basename = os.path.basename(image_name).replace('.fits','')
        sexcat_file = self.output_directory+'cat'+basename+'.txt'
        sex_command = sextractor_bin + " " + image_name + " -c " + self.config_directory + "cell.sex " + \
            " -CATALOG_NAME " + sexcat_file
        sex_command += ' -PARAMETERS_NAME  %scell.param ' % (self.config_directory)
        sex_command += ' -FILTER_NAME %sgauss_3.0_7x7.conv ' % (self.config_directory)
        if checkimages:
            sex_command += ' -CHECKIMAGE_TYPE APERTURES,SEGMENTATION,BACKGROUND '
            sex_command += ' -CHECKIMAGE_NAME %sap_%s.fits,%sseg_%s.fits,%sback_%s.fits ' % (self.output_directory,basename,self.output_directory,basename,self.output_directory,basename)
        else:
            sex_command += ' -CHECKIMAGE_TYPE NONE '
        print sex_command
        system(sex_command)
        return sexcat_file
    

            
    def ReadCatalogs(self):
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
        if not hasattr(self,"catlist"):
            print "catlist attribute doesn't exist; looking for already generated"
            print " catalogs in the output_directory."
            self.catlist = glob.glob(self.output_directory+'cat??????.txt')
        if len(self.catlist) == 0:
            print 'No *.cat files found in %s! Run FindDetections() first.' % (self.output_directory)
            return
        for cat_text_path in self.catlist:
            cat_text_file = os.path.basename(cat_text_path)
            catid = int(cat_text_file[3:].lstrip('0').rstrip('.txt'))
            tmp_table=pd.read_csv(cat_text_path,header=None,delim_whitespace=True,names=detection_names)
            tmp_table.z_index = catid
            if count == 0:
                detection_table = tmp_table # do this for the first item, then append on to it
            else:
                detection_table = detection_table.append(tmp_table,ignore_index=True)
            count += 1
        print str(count) + ' tables read in'
        self.detection_table = detection_table
        print "Ready for Object/Association assignment with FindObjects()"
    
    
    def FindObjects(self,xy_rad=20.0,slice_rad=20.0):
        '''
        Get database of associations

        Parameters
        ----------
        self.detection_table : pandas detection table of all detection across all slices
            generated by ReadCatalogs()
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
        for index,row in self.detection_table.iterrows():
            current_x = row['x']
            current_y = row['y']
            current_z = row['z_index']


            #######
            # STEP 1: Gather table of nearby detections within distance and slice threshold
            # self.detection_table[np.sqrt((current_x-self.detection_table.x)**2+(current_y-self.detection_table.y)**2) < xy_rad]
            # Another common operation is the use of boolean vectors to filter the data. 
            # The operators are: | for or, & for and, and ~ for not. These must be 
            # grouped by using parentheses.

            # Here is the logic for determining table of nearby detections
            nearby_detections_boolean = (np.sqrt((current_x-self.detection_table.x)**2+(current_y-self.detection_table.y)**2) \
                < xy_rad) & (np.abs(current_z - self.detection_table.z_index) < slice_rad) \
                & (current_z != self.detection_table.z_index) # added this to avoid double-counting of nearby cells in same frame 

            # get integer index of label index [just in case the lengths are different]
            current_int_index = np.where(self.detection_table.index == index)[0][0]

            # set boolean value to True for current index position so it is included
            nearby_detections_boolean[current_int_index] = True

            # extract the subset of the full detections 
            nearby_detections_table = self.detection_table[nearby_detections_boolean]  

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


            ## STEP 2A: Iterate until we find a set number of slices with a common
            # number of nearby detections. This fights against single outliers that 
            # only appear in a single image. Stop the iteration if we're at a 
            # single detection across slices.  
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

            ## STEP 2B: Only continue with the object assignment process if the object
            # is not assigned already.   
            # if all the detections in each of these slices is accounted for already, 
            # then the associated object with them should exist
            # if SOME detections are accounted for already but not all, 
            # this is an indication that something is wrong and the threshold is too
            # low.  

            unaccounted_objects = []
            for inddd in most_det_in_ind_slice_table.index:
                if inddd not in set_of_detections_accounted_for:
                    unaccounted_objects.append(inddd)

            most_det_in_ind_slice_table = most_det_in_ind_slice_table.loc[unaccounted_objects]

            # if not set(most_det_in_ind_slice_table.index).issubset(set_of_detections_accounted_for): 

            # else:
            #     for inddd in most_det_in_ind_slice_table.index:
            #         assert inddd not in set_of_detections_accounted_for
            #     continue


            ## STEP 3: Take one slice as a reference slice and compare each of the
            # detections in the other slices with the same number of nearby 
            # detections and assign the grouping of them to a single potential object.
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
                # common detection level. EDIT this is now taken care of higher up
                if len(check_table.det_id) < 3:
                    print "only %i instance(s) of %s object detections in slice %s at x=%f y=%f" % (len(check_table.det_id),str(len(reference_table)),str(int(ref_row['z_index'])),ref_row['x'],ref_row['y'])
                    # could instead flag it as a spurious detection
                    # or have a column in the object table of how many slices went into it

                # now that the detection ids of the closest objects have been found,
                # add these all to a table. The object table will have the means
                # and standard deviations of the z, x, and y positions of the 
                # DETECTIONS, which is the slices in which the objects were separated
                # from each other. Not to be confused with the ASSOCIATIONS, which 
                # is ones that overlap with a detection.
                obj_df = most_det_in_ind_slice_table.loc[check_table.det_id]
                means = obj_df.mean()
                stds = obj_df.std()
                # the slice at which mu_max [indicator of max surface brightness] 
                # appears to be a good indicator of when the cell is most in focus
                # and thus near the center of the current slice
                z_min_mu_max = obj_df.z_index.loc[obj_df.mu_max.idxmin()]
                object_names = ['z_mean','z_std','x_mean','x_std','y_mean','y_std','z_min_mu_max','detections','associations']     
                ind_obj_list = [means['z_index'],stds['z_index'],means['x'],stds['x'],means['y'],stds['y'],z_min_mu_max,len(check_table.det_id),len(check_table.det_id)]
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
                set_of_detections_accounted_for.update(set(obj_df.index))
                # if index == 1091: raise Exception
                # if index > 1051: raise Exception


            # # check if all associated items have already been accounted for by seeing if all detections
            # # in the object frame are already in the detections_accounted_for_set
            # if not set(nearby_detections_table.index).issubset(detections_accounted_for_set): 
            #     detections_accounted_for_set.update(set(nearby_detections_table.index))
            #     
            #     # build up list of potential objects. CHANGE THIS.
            #     pot_obj_list.append(nearby_detections_table)
            # 
            # # IMPLEMENT THIS    
            # #must have three unique NONBLENDED detections.. but still count the blended ones as once
            # 

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

        ## Consistency checks!
        # if all the detections of an object are in adjacent slices to each other, 
        # the standard deviation of the z_index 
        # should be equal to std(np.arange(N)) where N is the number of slices that 
        # the object appears.  

        print str(len(object_table)) + ' potential objects found'
        print str(len(set_of_detections_accounted_for))+ '/' + str(len(self.detection_table)) + ' detections accounted for'
        print 'Attempting to assign associations to the remaining detections'

        ## STEP 5: iterate through the remaining detections and associate each
        # object within its elliptical area to it. This covers the large detections
        # that may be a blurred combination of several objects.
        # copy the list of detections accounted for
        set_of_associations_accounted_for = set(copy.copy(list(set_of_detections_accounted_for)))
        for index,row in self.detection_table.iterrows():
            if index in set_of_associations_accounted_for:
                # only do these association calculations for unassigned detections
                continue
            else:
                pass
            current_x = row['x']
            current_y = row['y']
            current_z = self.detection_table.z_index[index]
            current_A = row['A'] # semimajor axis
            current_B = row['B']
            current_kron = row['r_kron'] 
            current_theta = row['theta']

            # A*r_kron is the semimajor axis shown on the verification images.
            # 25% of this radius should be more than enough for an association
            # could make this user-specifyable for more tweaking ability. 
            # making it too big will lead to an overload of associations.
            association_radius_thresh = 1
            kron_A = current_kron*current_A*association_radius_thresh
            kron_B = current_kron*current_B*association_radius_thresh

            # initial truncation of possible assotiations; to be within an ellipse,
            # the object must be within the box centered at the origin of the 
            # ellipse with sides = 2*semimajor axis
            # [actually more constraining is that it is within the circle of 
            # radius r=semimajor_axis, but this requires an additional calculation]
            possible_associations_table = object_table.loc[(abs(object_table.x_mean-current_x) < kron_A) & \
                (abs(object_table.y_mean-current_y) < kron_A)]

            verified_indices = []
            maximum_association_distance = 7
            for poss_association_ind, poss_association_row in possible_associations_table.iterrows():
                potential_association_z_list = list(self.detection_table.loc[detection_object_link_table.loc[detection_object_link_table.obj_id == poss_association_ind].det_id].z_index)
                if current_z in potential_association_z_list:
                    print "In slice %i, not associating detection %i with object %i as it exists already in that slice." % (current_z,index,poss_association_ind)
                    continue # don't make an association with an object that already exists in the current slice
                z_distances = abs(np.array(potential_association_z_list) - current_z)
                if min(z_distances) > maximum_association_distance:
                    print "In slice %i, detection %i is close to object %i in x,y space but is %i slices away." % (current_z,index,poss_association_ind,min(z_distances))
                    continue
                ellipse_of_detection = Ellipse(0,0,kron_A,kron_B)
                xyang=np.arctan((poss_association_row['y_mean']-current_y)/(poss_association_row['x_mean']-current_x))
                px,py = ellipse_of_detection.pointFromAngle(xyang,theta=current_theta)
                # if the distance from the center of the ellipse to px,py is 
                # further than the distance to x_mean,y_mean of the object, this
                # means that the object is within the elliptical threshhold.
                px = px[0]
                py = py[0]
                dist_to_ellipse_edge = np.sqrt((px)**2 + (py)**2)
                dist_to_potential_object = np.sqrt((poss_association_row['x_mean']-current_x)**2 + (poss_association_row['y_mean']-current_y)**2)

                if dist_to_ellipse_edge < dist_to_potential_object:
                    print "In slice %i, object %i is nearby detection %i but is outside its elliptical threshold" % (current_z,poss_association_ind,index)
                    continue
                object_table.associations.loc[poss_association_ind] += 1 # add that we have an association

                verified_indices.append(poss_association_ind)

            verified_associations_table = possible_associations_table.loc[verified_indices]
            #print possible_associations_table
            ass_det_link_dict = {'obj_id':verified_associations_table.index,'det_id':index}
            tmp_table3 = pd.DataFrame(ass_det_link_dict)
            association_object_link_table = association_object_link_table.append(tmp_table3,ignore_index=True)
            if len(verified_associations_table.index) > 0:
                set_of_associations_accounted_for.add(index)
            else:
                print 'Could not find an association for detection %i in slice %i' % (index, current_z)
        print str(len(set_of_associations_accounted_for)) + '/' + str(len(self.detection_table)) + ' detections associated'
        
        self.total_detections = len(self.detection_table)
        self.total_objects = len(object_table)
        self.total_detections_assigned_to_objects = len(set_of_detections_accounted_for)
        self.total_detections_associated_with_objects = len(set_of_associations_accounted_for)
        self.object_table = object_table
        self.detection_object_link_table = detection_object_link_table
        self.association_object_link_table = association_object_link_table
        print "Object/association table generation complete."
        print "Can create verification images [MakeVerificationImages()]"
        print "Or generate histograms [PlotCountHist()]"

    def PrepareVerificationImages(self):
        '''This takes several minutes, recommended only if you need to verify that
        your detections are being populated correctly 
        '''
        print "Converting aperture images to scaled pngs. This will take a long time."
        ap_list = glob.glob(self.output_directory+'ap_0*.fits')
        if len(ap_list) == 0:
            print "Can't find any ap_0*.fits files. Rerun FindDetections() with checkimages=True?"
            return
        for ap_img in ap_list:
            png_name = ap_img.replace('fits','png')
            im_cmd = 'convert %s -fx "log(abs((u.r+u.g+u.b)*80))" %s' % (ap_img,png_name)
            os.system(im_cmd)
            # cp_cmd = 'cp %s %s' % (png_name,png_name.replace('ap_','detections_'))
            # os.system(cp_cmd)
        print "Images converted. Ready for MakeVerificationImages()"
    
    def _create_verification_images(self,link_table,textcolor='gold',overlay=False):
        detection_table = self.detection_table
        object_table = self.object_table
        
        
        #convert ap_000017.fits -fx "log(abs((u.r+u.g+u.b)*80))" testlog.png
        if not overlay:
            ap_list = glob.glob(self.output_directory+'ap_0*.png')
        else:
            ap_list = glob.glob(self.output_directory+'detections_0*.png')
        
        if len(ap_list) == 0:
            print "Need to convert the sextractor verification fits files to pngs first using PrepareVerificationImages()"
            return
            
        for ap_img in ap_list:
            if not overlay:
                base = os.path.basename(ap_img)
                num = int(base.replace('ap_','').lstrip('0').replace('.png',''))
                out_img = ap_img.replace('ap_','detections_')
            else:
                base = os.path.basename(ap_img)
                num = int(base.replace('detections_','').lstrip('0').replace('.png',''))
                out_img = ap_img

            # get the detection table from the 5th image slice
            current_img_detections = detection_table[detection_table.z_index==num]
            # join it with the detection_object link table
            current_link = link_table.join(current_img_detections.z_index,on='det_id',how='right')

            # join it with the object stats table
            current_objects=current_link.join(object_table,on='obj_id',how='left')
            # loop through the objects and make a label for each
            im_cmd = 'convert ' + ap_img + ' -scale 200% -gravity Center ' 
            if overlay:
                im_cmd = 'convert ' + ap_img + ' -gravity Center ' #dont double the size

            ### here we're centering on each DETECTION and printing a list of each OBJECT that detection is associated with
            for index,row in current_img_detections.iterrows():
                matching_objects = current_objects[current_objects['det_id'] == index]
                obj_list_str = str(list(matching_objects.obj_id)).rstrip(']').lstrip('[').replace('.0','')
                if obj_list_str == 'nan': 
                    color = 'red'
                    if overlay:
                        continue
                else:
                    color = textcolor
                im_cmd+=self._overlay_verification_text(obj_list_str,row['x'],row['y'],row['theta'],color=color)
            ### here we're centering on the average position of each OBJECT and printing it at this average position
            # for index, row in current_objects.iterrows():
            #     im_cmd+=_overlay_verification_text(str(row['obj_id']),row['x_mean'],row['y_mean'],0)
            im_cmd += out_img 
            os.system(im_cmd)
        print im_cmd
        #convert testlog.png -fill red -annotate %fx%f%s%f%s%f "10" testannotate.png
        # % (ang,ang,xsign,xpos,ysign,ypos)

    
    def _overlay_verification_text(self,obj_str,xx,yy,theta,color='red',double_size=True):    

        if theta <= 0: 
            ang = abs(theta)
        else: 
            ang = 360. - theta

        if not double_size:
            xpos = xx - 256.  #subtracting 256 since changing gravity setting not only 
            # refers to the position on the background image, but also the position 
            # in the overlay image that is to be drawn.
            ypos = 512.- yy -256.
            ptsize = 10
        else:
            xpos = xx*2 - 512.
            ypos = 1024. - yy*2 - 512
            ptsize = 20

        if xpos < 0:
            xsign="-"
        else:
            xsign="+"
        if ypos < 0:
            ysign="-"
        else:
            ysign="+"


        im_cmd_append = ' -fill %s -stroke black -pointsize 30 -annotate %fx%fXSIGN%fYSIGN%f "%s" ' % (color,ang,ang,abs(xpos),abs(ypos),obj_str)
        im_cmd_append = im_cmd_append.replace('XSIGN',xsign)
        im_cmd_append = im_cmd_append.replace('YSIGN',ysign)
        return im_cmd_append

        # get all the detections associated with object 60
        # det.loc[link.loc[link.obj_id==60].det_id]


        #convert testlog.png -fill red -annotate %fx%f%s%f%s%f "%i" testannotate.png
        # % (ang,ang,xsign,xpos,ysign,ypos,det_index)
    
    
    def MakeVerificationImages(self):
        if not hasattr(self,"object_table"):
            print "Need to generate object/association tables with FindObjects() first."
            return
        self._create_verification_images(self.association_object_link_table,textcolor='gold')
        self._create_verification_images(self.detection_object_link_table,textcolor='cyan',overlay=True)
    
    
    def PlotCountHist(self):
        if not hasattr(self,"object_table"):
            print "Need to generate object/association tables with FindObjects first."
            return
        object_table = self.object_table
        savename=self.name+'_hist.png'
        plt.bar(np.round(object_table.z_mean).value_counts().index,np.round(object_table.z_mean).value_counts(),color='blue',alpha=0.5)
        plt.bar(object_table.z_min_mu_max.value_counts().index,object_table.z_min_mu_max.value_counts(),color='red',alpha=0.5)
        plt.legend(['mean slice #','minimum $\mu_{max}$ slice #'],'upper right')
        plt.ylabel('Counts')
        plt.xlabel('Slice #')
        plt.title(self.name)
        plt.savefig(self.output_directory+savename)
    
    def ret_detections_from_obj_id(self, obj_id):
        if not hasattr(self,"detection_object_link_table"):
            print "Need to generate object/association tables with FindObjects first."
            return
        full_det_table = self.detection_table
        detection_object_link_table = self.detection_object_link_table
        return full_det_table.loc[detection_object_link_table.loc[detection_object_link_table.obj_id == obj_id].det_id]
    
    

def make_ass_cat(image_name, assoc_image_name):
    '''Create the ASSociation CATalaog. Depreciated??'''
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
    sex_command = sextractor_bin + " " + image_name + " -c " + config_directory + "cell_assoc.sex " + \
        " -CATALOG_NAME " + sexcat_file + " -ASSOC_NAME " + assoc_catname 
    print sex_command
    system(sex_command)
    return sexcat_file


# def FindDetections(checkimages=False):
#     # give previous_image as None if first in the list; give following_image as None if last
#     fits_list = glob.glob('0*.fits')
#     for image in fits_list:
#         make_sex_cat(image,checkimages=checkimages)
#     # assoc_list = [None] + fits_list + [None]
#     # for index in np.arange(len(fits_list))+1:
#     #     make_ass_cat(assoc_list[index], assoc_list[index-1])
#     #     make_ass_cat(assoc_list[index], assoc_list[index+1])
#     # 
#     # assuming file format of cat000119.txt
#     # catlist = glob.glob('cat0*.txt')
#     # if len(catlist) == 0:
#     #     raise IOError("List of catalogs is empty. Check path and verify presence of cat0xxxxx.txt files")
#     # for cat in catlist:
#     #     catid = int(item[3:].lstrip('0').rstrip('.txt'))
#     #     
# 
#     # print catnumlist





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
 



from math import copysign, cos, sin, sqrt
class Ellipse:

    def __init__(self, mx, my, rh, rv):
        self.mx = mx # x center
        self.my = my # y center
        self.rh = rh # semimajor axis
        self.rv = rv # semiminor axis
       
    def pointFromAngle(self, xyangle, theta=0):
        '''
        Determine the point on the ellipse that intersects the ray from the
        origin of the ellipse to another point.  
        
        xyangle is the angle between the x axis and the point of interest
            = arctan((y-y_0)/(x-x_0))
            where (x_0,y_0) is the origin of the ellipse and (x,y) is the 
            point of interest
        
        theta is the rotation angle of the ellipse from the x axis.  
        
        see http://i.imgur.com/2m7ymjC.png for a graphical explination 
        '''
        psi = xyangle - theta
        c = np.cos(psi)
        s = np.sin(psi)
        ta = s / c  ## tan(psi)
        tt = ta * self.rh / self.rv  ## tan(t)
        d = 1. / np.sqrt(1. + tt * tt)
        x = self.mx + np.copysign(self.rh * d, c)
        y = self.my + np.copysign(self.rv * tt * d, s)
        try:
            len(x)
        except:
            x=np.array([x])
            y=np.array([y])
        
        if theta != 0:
            #rotate
            xlist=[]
            ylist=[]
            xyarr = np.array(zip(x,y))
        
            for xy in xyarr:
                rot_matrix = np.array([[np.cos(theta),-1*np.sin(theta)],[np.sin(theta),np.cos(theta)]])
                xy_prime = rot_matrix.dot(xy.T)
                xlist.append(xy_prime[0])
                ylist.append(xy_prime[1])
            x=np.array(xlist)
            y=np.array(ylist)
        return x, y
        
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
    



def loadPickle(pklpath):
    '''Given an input object and an output path, load a pickle file.'''
    if os.path.exists(pklpath):
        storefile=open(pklpath)
        loadedpkl = pickle.load(storefile)
        storefile.close()
        print "Loaded pickle file for this object."
        return loadedpkl
    else:
        print "Pickle file %s does not exist" % pklpath
        return None

def savePickle(input,outpath,clobber=False):
    path_existed = os.path.exists(outpath)
    if path_existed and not clobber:
        print '%s already exists and clobber == False; not overwriting.' % (outpath)
        return
    else:
        storefile = open(outpath,'w')
        pickle.dump(input,storefile)
        storefile.close
        if not path_existed:
            print "No Pickle file existed, so one was created:"
            print outpath
        else:
            print "Overwrote Pickle file:"
            print outpath    
    