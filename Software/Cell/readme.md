Quantifying 3D cell migration from intact monolayers into hydrogels
====================

Scientific lead: *Mariana B. Garcia* (UC Berkeley Vision Science)
Software development: *Adam N. Morgan* (UC Berkeley Astrophysics)

### Background

![Assay Preparation](http://i.imgur.com/8Lq6iJr.png)

Our pilot study consisted of four replicates of each condition, resulting in 256 Z-stacks (4 imaged areas/gel, 4 gels/condition, 4 conditions/day, 4 days). Each hydrogel area was imaged to a depth of about 2mm, requiring roughly 100-300 10um slices.  Each of these 3D image stacks contains anywhere from a few to over a hundred cells, each appearing over about 3-10 of the 10 um image slices.  This results in a total of over *50,000* images containing more than *8000* cells (with over *60,000* detections of those cells over all image slices). 

### Automated Image Processing 

![Image Processing Steps](http://i.imgur.com/CydObU5.png)

Every image [A] in each stack was processed with our custom Python code wrapped around Source Extractor (*http://www.astromatic.net/software/sextractor), a source detection algorithm developed to identify and analyze galaxies in astronomical images. The processing steps include:
1. Background measurement and subtraction [B]
2. Smoothing the image with a gaussian filter
3. Object identification through thresholding 
4. Deblending/segmentation of nearby objects [C]
5. Measurement of object shapes and positions [D]
6. Brightness measurement and catalog output

The resultant output database for each image stack was analyzed with a custom cell association algorithm written for this study.  The code iteratively searches through the hundreds to thousands of source detections to determine which comprise the same cell, and then measures the image location and migration depth of each identified cell (shown in blue numbers in the figure to the right).  The algorithm accounts for *crowded sources*, *possible associations* (one or more cells that likely become better resolved further down the stack; shown in yellow), and *spurious detections* (shown in red).

![Cell Count Identification](http://i.imgur.com/T66iJ4U.png)

![Cell Count Animation](http://i.imgur.com/85MQl0o.gif)