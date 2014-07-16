Automated image analysis pipeline to quantify biological cell migration
====================

* Scientific lead: *Mariana B. Garcia* (UC Berkeley Vision Science)
* Software development: *Adam N. Morgan* (UC Berkeley Astrophysics)

Background
-----------

![Assay Preparation](http://i.imgur.com/8Lq6iJr.png)

Our pilot study consisted of four replicates of each condition, resulting in 256 Z-stacks (4 imaged areas/gel, 4 gels/condition, 4 conditions/day, 4 days). Each hydrogel area was imaged to a depth of about 2mm, requiring roughly 100-300 10um slices.  Each of these 3D image stacks contains anywhere from a few to over a hundred cells, each appearing over about 3-10 of the 10 um image slices.  This results in a total of over **50,000** images containing more than **8000** cells (with over **60,000** detections of those cells over all image slices). 

Automated Image Processing 
--------------------------

Every image [A] in each stack was processed with our custom Python code wrapped around [Source Extractor][1], a source detection algorithm developed to identify and analyze galaxies in astronomical images. The processing steps include:

1. Background measurement and subtraction [B]
2. Smoothing the image with a gaussian filter
3. Object identification through thresholding 
4. Deblending/segmentation of nearby objects [C]
5. Measurement of object shapes and positions [D]
6. Brightness measurement and catalog output

![Image Processing Steps](http://i.imgur.com/CydObU5.png)

The resultant output database for each image stack was analyzed with a custom cell association algorithm written for this study.  The code iteratively searches through the hundreds to thousands of source detections to determine which comprise the same cell, and then measures the image location and migration depth of each identified cell (shown in blue numbers in the figure to the right).  The algorithm accounts for **crowded sources**, **possible associations** (one or more cells that likely become better resolved further down the stack; shown in yellow), and **spurious detections** (shown in red).

![Cell Count Identification](http://i.imgur.com/nBSwbbk.png)

This allows for the rapid and automated analysis of cell count as a function of depth in the hydrogel, allowing for quick measurement of cell migration.  In addition to eliminating the speed bottleneck, this also improves the quality and consistency of the image analysis, freeing it from potential human error.

For example, the animation below illustrates how this method was used to show how a coverslip covered with many cells isn't perfectly level in the x-y plane, manifesting itself as a line of object detections passing through a number of z-stacks. 

![Cell Count Animation](http://i.imgur.com/85MQl0o.gif)


Running the Code
------------

### Dependencies
* This code: [https://github.com/qmorgan/qsoft/tree/master/Software/Cell](https://github.com/qmorgan/qsoft/tree/master/Software/Cell)
* ImageMagick: [http://www.imagemagick.org/](http://www.imagemagick.org/)
* Source Extractor: [http://www.astromatic.net/software/SExtractor](http://www.astromatic.net/software/SExtractor)
* PyLSM: [https://pypi.python.org/pypi/pyLSM](https://pypi.python.org/pypi/pyLSM)

### Image Preparation 
* Initialize ImageStack object

    `imgstack = Migration.ImageStack(image_directory=inpath,output_directory=outpath,config_directory=configpath)`
    
A stack of images at a particular x-y location (a "z-stack") defines an ImageStack object upon which further operations are performed.  This is initialized by defining the paths of the original image locations, the desired output locations, and the location of the SExtractor configuration files. 

The configuration file output by the microscope imaging program ZSeries*.xml is loaded and read with the `read_xml_configuration()` method.
    

* Convert the images to .fits files:

    `imgstack.PrepareImages()`

The images for a z-stack are originally .tif files located in a
single directory. This function loads all ZSeries*.tif files within
self.image_directory and converts them to .fits files using 
[ImageMagick](http://www.imagemagick.org/)
    

### Source Detection
* Run SExtractor on the images
    `imgstack.FindDetections(checkimages=True)`

Loops through each image in the stack and and run SExtractor on 
each one. 

The following SExtractor configuration files are loaded:

 * `cell.sex` - SExtractor configuration
 * `cell.param` - parameters to include in SExtractor output catalog 
 * `gauss_3.0_7x7.conv` - image background filter

Set `checkimages=True` for additional calibration output images: 

 *  `back_%s.fits` - measured background, subtracted from final images
 * `seg_%s.fits` - segmentation image, showing detection separation
 * `ap_%s.fits` - shows elliptical aperatures calculated for each object
    
* Read in the SExtractor output into a table

    `imgstack.ReadCatalogs()`
    
Builds up the full detection table by reading in all the ascii outputs from 
Source Extractor.  

### Object/Association Assignment 
* Run the object/association assignment algorithm

    `imgstack.FindObjects()`

  The steps undergone are as follows:

  1. Gather table of nearby detections within distance and slice threshold. Must: 
     * be under the xy-distance search radius: `(np.sqrt((current_x-self.detection_table.x)**2+(current_y-self.detection_table.y)**2) < xy_rad)`
     * be under the z-distance search radius: `(np.abs(current_z - self.detection_table.z_index) < slice_rad)`
     * not double-counting nearby cells in the same frame: `current_z != self.detection_table.z_index`
     * find more than 3 detections over all stacks for a particular object to avoid spurious detections

  2. Find the z-positions in the nearby detection table that have the most detections.
     * Iterate until we find a set number of slices with a common number of nearby detections. This fights against single outliers that only appear in a single image. Stop the iteration if we're at a single detection across slices. In each iteration, we create a list of the modes of slices (if there are multiple) by counting the frequency of each slice in the detection list. Then we get the subset table of detections have the most number of nearby detections in a single slice.
     * Only continue with the object assignment process if the object is not assigned already. If all the detections in each of these slices is accounted for already, then the associated object with them should exist. If *some* detections are accounted for already but not all, this is an indication that something is wrong and the threshold is too low.  
   
  3. Take one slice as a reference slice and compare each of the
  detections in the other slices with the same number of nearby 
  detections and assign the grouping of them to a single potential object.
  Take the middle slice from this list as the reference point for 
  object assignment. The middle slice is most likely to have the 
  clearest detections of the objects and best positions.
    * Each reference detection will get assigned to an object.
  We create a table with just the reference info so far, and then add detections. If no detection from that slice has been found yet, append it to the list of objects. 
    * now that the detection ids of the closest objects have been found,
      add these all to a table. The object table will have the means
      and standard deviations of the z, x, and y positions of the 
      **detections**, which is the slices in which the objects were separated
      from each other. Not to be confused with the **associations**, which 
      is ones that overlap with a detection.
    * The z-location for a particular object can be defined by either the slice at which `mu_max` [indicator of max surface brightness] is a minimum, or the mean slice number of all slice locations for a particular object. Both appear to be a good indicator of when the cell is most in focus and thus the actual average z-position of the cell.
   
  4. Now that all objects have been assigned, we iterate through the remaining detections and assign *associations* for each
     object within its elliptical area to it. This covers the large detections
     that may be a blurred combination of several objects. 
    * The semimajor/semiminor axes of the ellipses as measured from the SExtractor output is multiplied by its Kron radius ([Kron 1980][2]) and an additional scaling factor to determine the association threshold.  The Kron radius is a brightness-weighted distance which defines the "first moment" of the object, defined by ![Kron](http://i.imgur.com/2cCJbaC.gif) 


### Analysis

* Plot count histogram for the stack: 

    `imgstack.PlotCountHist()`
    
Plots the number of objects as a function of depth for an image stack. Plots
for two definitions for the "center" z-position for a given object (the mean
slice number for a given object, and the slice number containing the minimum
value of `mu_max` [indicator of max surface brightness]. Both give generally
consistent values, but we recommend using the `mu_max` definition to further
defend against outliers that may lead to skewed values of the mean slice
number.
    
![Single Stack Histogram](http://i.imgur.com/CP6LeYX.png)
    
* Plot combined count histograms for many stacks: 
    
    `Migration.CombinedHist(objectlist,outname,title)`
    
Plots the results from all images for each gel. In the plot below, each of the
4 gels (which can be thought of as a "trial" for a given set of conditions) is
plotted with a different base color (red, yellow, blue, purple), and each of
the 4 different imaged areas for that gel is plotted as a different shade of
that base color. If the number of gels and imaged areas for each gel is large
enough, this will give an approximation of the true distribution of cell
migration for a given set of conditions.

![Multi-Stack Histogram](http://i.imgur.com/0O3YjUg.png)

### Saving/Exporting 
   
Options to save/load the Image Stack Object are provided for future analysis.
In addition, you can export the object detection/association tables to an
Excel file if desired.
   
* Save the ImageStack object

    `imgstack.Save()`

* Load the ImageStack object

    `imgstack = Migration.loadPickle("PathToSavedFile.pkl")`

* Optional: export the ImageStack tables to an Excel file, if desired.

    `imgstack.ExportToExcel()`

The exported tables include: 

  * `self.object_table`
  * `self.detection_table`
  * `self.detection_object_link_table`
  * `self.association_object_link_table`

### Image Verification

The source detection, object assignment, and association parameters have been optimized for this specific study, and may need to be changed for different experiments. To check that the algorithm is performing as expected, verification images can be made for each image slice within a stack. In the verification image shown [above](http://i.imgur.com/chjZ7Sx.png), ellipses are drawn over each detection in that slice (dotted lines indicate a nearby source). If this detection corresponds to a single object, its object ID is overplotted in blue text. If multiple objects may be associated with this detection, their object IDs are plotted in yellow. If no object is associated with a detection, "nan" is plotted in red (likely outliers). 

This process takes several minutes, and is recommended only if you need to verify that
your detections are being populated correctly. The steps are as follows:

* Prepare verification images (convert .fits files to scaled .png files)

    `imgstack.PrepareVerificationImages()`

* Create verification images (overlay the object IDs)

    `imgstack.MakeVerificationImages()`
    
    
Acknowledgements
------------

[1]: http://www.astromatic.net/software/SExtractor     "Source Extractor"
[2]: http://adsabs.harvard.edu/cgi-bin/bib_query?1980ApJS...43..305K "Kron (1980)"

Financial support for this study was provided by National Institutes of Health Grant R01 EY012392
