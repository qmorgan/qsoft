import sys
import time
import os
try:
    from PIL import Image
    from PIL import ImageOps
except:
    sys.exit('You do not have PIL.  Download it.')
import base64, Image, string
import glob
from MiscBin.q import dec2sex
from MiscBin.q import sex2dec
from MiscBin import qErr

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have WCSTOOLS installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

class qImage:
    '''Instance of a qImage class'''
    def __init__(self,image_path=None):
        self.image_path = image_path
        self.survey_str = 'Unknown'
    
    def dss_grab(self,ra,dec,size,survey="dss2red"):
        
        if survey=="dss2red":
            self.survey = 'poss2ukstu_red'  # 'poss2ukstu_ir', 'poss2ukstu_blue', 'poss1_red', 
        url_str =  'http://stdatu.stsci.edu/cgi-bin/dss_search?v=' + self.survey + '&r=' + str(ra) + '&d=' + str(dec) +'&e=J2000&h=' + str(size) +'&w=' + str(size) +'&f=gif&c=none&fov=NONE&v3='
        self.img_ra = ra
        self.img_dec= dec
        self.wcs_size=size #in arcmin
        
        self.image_path=storepath+'dss_grab.gif'
        
        downloadImage(url_str,self.image_path)
    
    def sdss_grab(self,ra,dec,scale,width):
        url_str = 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra=' + str(ra) + '&dec=' + str(dec) + '&scale=' + str(scale) + '&width=' + str(width) + '&height=' + str(width) + '&opt=L'
        self.img_ra = ra
        self.img_dec = dec
        self.wcs_size= scale * width / 60.0 # arcminutes 
        self.survey = 'sdss'

        image_path=storepath+'sdss_grab.jpeg'
        downloadImage(url_str,image_path)      
        
        # PIL sometimes can't handle jpeg.  conver to png.
        newpath = storepath + 'sdss_grab.png'
        magickcommand = "convert %s %s" % (image_path, newpath)
        try:
            os.system(magickcommand)
            self.image_path = newpath
        except:
            print "Do not have Imagemagick, cannot convert svg to png"
    
    def invert_and_resize(self,newsize=600):
        '''
        Given a position and uncertainty, create a finding chart from the DSS
        (Move this to another module)
    
        Requires PIL
        http://pythonmac.org/packages/py25-fat/index.html
        http://www.p16blog.com/p16/2008/05/appengine-installing-pil-on-os-x-1053.html
    
        '''
        # SDSS - later http://casjobs.sdss.org/ImgCutoutDR5/getjpeg.aspx?ra=264.191&dec=-25.212&scale=0.8&width=750&height=750
    
        self.img_size = newsize # in pixels

        try:
            im1 = Image.open(self.image_path)
        except:
            errstring = 'Could not open %s; image may be malformed.' % (self.image_path)
            raise IOError(errstring)

        # Invert Colors
        im1 = ImageOps.invert(im1)
        # Convert to color mode
        im1 = im1.convert('RGB')
        # Resize with a cubic spline - it looks nice.
        # Other options are NEAREST, BILINEAR, ANTIALIAS
        im1 = im1.resize((newsize,newsize),Image.BICUBIC)
        self.save_str = storepath+'FCBase.png'
        im1.save(self.save_str)

            
    def overlay_finding_chart(self,ra,dec,uncertainty,src_name='unknown source',pos_label='UVOT',cont_str='',\
        uncertainty_shape='combo',suppress_scale=False):
        
        # cont_str='Contact: Adam N. Morgan (qmorgan@gmail.com, 510-229-7683)'
        
        self.loc_ra = ra
        self.loc_dec = dec
        self.loc_uncertainty=uncertainty
        self.src_name = src_name
        self.pos_label = pos_label
        self.cont_str = cont_str
        
        self.pixel_scale = self.wcs_size * 60.0 / self.img_size # arcsec/pixel
        
        self.str_loc_pos = dec2sex((self.loc_ra,self.loc_dec))
        
        update_time = time.strftime('%Y %h %d')
        
        head_buffer = 50
        uncertainty_circle_size = int(self.loc_uncertainty/self.pixel_scale)
        pos_pos_x = 300 #position of the error circle on the image , hard coded for now
        pos_pos_y = 350
        
        # Define scale line physical length based on how large the image is
        # Divide by 5 so that line is never more than ~20% of image
        if not self.wcs_size < 5.0:
            scale_line_physical_length = int(self.wcs_size/5) #arcmin
        elif not self.wcs_size < 2.0:
            scale_line_physical_length = 0.5
        else:
            scale_line_physical_length = 0.25
        scale_line_length = int(60*scale_line_physical_length/self.pixel_scale) 
        scale_line_start = head_buffer + self.img_size/2 - scale_line_length/2
        scale_line_stop = head_buffer + self.img_size/2 + scale_line_length/2
        scale_line_label = str(scale_line_physical_length) + ' arcmin'
        # Get rid of scale bar if desired
        if suppress_scale:
            scale_line_stop = scale_line_start
            scale_line_label = ''
        scale_line_label_pos_y = scale_line_start - 15
        
        pos_label_pos_y = pos_pos_y - uncertainty_circle_size - 15 
        
        epoch = 'J2000'
        pos_label_full = self.pos_label + ' Pos. (%s)' % epoch        
               
        if self.survey == 'poss2ukstu_red': 
            self.survey_str = 'DSS2 (red)'
        elif self.survey == 'poss1_red':
            self.survey_str = 'DSS1 (red)'
        elif self.survey == 'sdss':
            self.survey_str = 'SDSS'
        else: 
            print 'Unknown survey: %s' % survey
            sys.exit(1)

        ra_str = 'RA = ' + self.str_loc_pos[0]
        dec_str = 'Dec = ' + self.str_loc_pos[1]
        uncertainty_str = 'Uncertainty: ' + str(self.loc_uncertainty) + '"'
        
        if uncertainty_shape.lower() == 'circle':
            uncertainty_shape_string = '''
            <ellipse cx="%d" cy="%d" rx="%d" ry="%d"
            style="fill:none;
            stroke:rgb(0,0,255);stroke-width:2"/>
            ''' % (pos_pos_x, pos_pos_y, uncertainty_circle_size,uncertainty_circle_size)
        elif uncertainty_shape.lower() == 'cross': 
            cross_left_pos = pos_pos_x-uncertainty_circle_size-15
            cross_right_pos = pos_pos_x-uncertainty_circle_size
            cross_down_pos = pos_pos_y+uncertainty_circle_size
            cross_up_pos = pos_pos_y+uncertainty_circle_size+15
            uncertainty_shape_string = '''
            <line x1="%d" y1="%d" x2="%d" y2="%d"
            style="stroke:rgb(0,0,255);stroke-width:2"/>
            <line x1="%d" y1="%d" x2="%d" y2="%d"
            style="stroke:rgb(0,0,255);stroke-width:2"/>
            '''% (pos_pos_x, cross_down_pos, pos_pos_x, cross_up_pos,\
            cross_right_pos, pos_pos_y, cross_left_pos, pos_pos_y)
        else:
            #default to cross + circle
            cross_left_pos = pos_pos_x-uncertainty_circle_size-15
            cross_right_pos = pos_pos_x-uncertainty_circle_size
            cross_down_pos = pos_pos_y+uncertainty_circle_size
            cross_up_pos = pos_pos_y+uncertainty_circle_size+15
            uncertainty_shape_string = '''
            <ellipse cx="%d" cy="%d" rx="%d" ry="%d"
            style="fill:none;
            stroke:rgb(0,0,255);stroke-width:2"/>
            <line x1="%d" y1="%d" x2="%d" y2="%d"
            style="stroke:rgb(0,0,255);stroke-width:2"/>
            <line x1="%d" y1="%d" x2="%d" y2="%d"
            style="stroke:rgb(0,0,255);stroke-width:2"/>
            '''% (pos_pos_x, pos_pos_y, uncertainty_circle_size,uncertainty_circle_size,pos_pos_x, cross_down_pos, pos_pos_x, cross_up_pos,\
            cross_right_pos, pos_pos_y, cross_left_pos, pos_pos_y)
        
        im = Image.open(self.save_str)
        assert isinstance(im, Image.Image)
        im.fp.seek(0)
        imstr = base64.b64encode(im.fp.read())
        imstr = [ imstr[i:i+76] for i in range(0, len(imstr), 76) ]
        imstr = string.join(imstr, '\n')
        
        fc_width = im.size[0] + 200.0
        fc_height = im.size[1] + 100.0
        
        return '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
           width="%d" height="%d">
             <image y="%d" x="0"   width="%d" height="%d" xlink:href="data:image/%s;base64,%s" />
        

         <text x = "15" y = "10" fill = "black" font-size = "28">
         <tspan x="15" dy="1em">%s: Finding Chart</tspan>
         </text>
         <text x = "15" y = "660" fill = "black" font-size = "16">
         <tspan x="15" dy="1em">%s</tspan>
         </text>
         <text x = "15" y = "630" fill = "#222222" font-size = "11">
         <tspan x="15" dy="1em">Created from http://fc.qmorgan.com on %s.  Contact: Adam N. Morgan (qmorgan@gmail.com).</tspan>
         </text>

         <text x ="288" y="%d" fill="blue" font-size="12">
         <tspan x="288" dy="1em">%s</tspan>
         </text>
         
         %s

         <text x="615" y ="55" fill = "black" font-size = "18">
         <tspan x="615" dy="1em">%s</tspan>
         <tspan x="615" dy="1em">RA = %s</tspan>
         <tspan x="615" dy="1em">Dec = %s</tspan>
         <tspan x="615" dy="1em">Uncertainty = %s</tspan>
         <tspan x="615" dy="1em">(arcseconds)</tspan>
         </text>

         <text x="640" y="325" fill="black" font-size="18">
         <tspan x="640" dy="1em">Image from:</tspan>
         <tspan x="640" dy="1em">%s</tspan>
         </text>

         <text x="640" y="500" fill="black" font-size="18">
         <tspan x="640" dy="1em">Image is %s</tspan>
         <tspan x="640" dy="1em">arcmin/side</tspan>
         </text>

         <text x="625" y="640" fill="black" font-size ="18">
         <tspan x="625" dy="1em">E</tspan>
         </text>
         <line x1="640" y1="650" x2="700" y2="650" style="stroke:Black;stroke-width:2"/>
         <line x1="640" y1="650" x2="650" y2="660" style="stroke:Black;stroke-width:2"/>
         <line x1="640" y1="650" x2="650" y2="640" style="stroke:Black;stroke-width:2"/>

         <text x="694" y="570" fill="black" font-size ="18">
         <tspan x="694" dy="1em">N</tspan>
         </text>
         <line x1="700" y1="650" x2="700" y2="590" style="stroke:Black;stroke-width:2"/>
         <line x1="690" y1="600" x2="700" y2="590" style="stroke:Black;stroke-width:2"/>
         <line x1="710" y1="600" x2="700" y2="590" style="stroke:Black;stroke-width:2"/>


         <line x1="130" y1="%d" x2="130" y2="%d"
         style="stroke:rgb(255,0,0);stroke-width:2"/>
         <text x ="107" y ="%d" fill ="red" font-size = "12">
         <tspan x="107" dy="1em">%s</tspan>
         </text>
        
        </svg>''' % ( fc_width, fc_height, head_buffer, im.size[0], im.size[1], im.format.lower(), imstr,\
                self.src_name,  self.cont_str, update_time, pos_label_pos_y, self.pos_label, uncertainty_shape_string,\
                 pos_label_full, self.str_loc_pos[0], self.str_loc_pos[1],\
                 str(self.loc_uncertainty),self.survey_str, str(self.wcs_size), scale_line_start,scale_line_stop,scale_line_label_pos_y, scale_line_label)
    



def stealStuff(file_name,file_mode,base_url, timeout=20, verbose=False):
    from urllib2 import Request, urlopen, URLError, HTTPError
    
    #create the url and the request
    url = base_url + file_name
    req = Request(url)
    successful_download = False
    count = 0
    
    # Open the url:
    while not successful_download and count < 6:
        count += 1
        trys_left = 5-count
        try:
            f = urlopen(req, timeout=timeout)
            if verbose: 
                print "downloading " + url

            # Open our local file for writing
            local_file = open(file_name, "w" + file_mode)
            #Write to our local file
            local_file.write(f.read())
            local_file.close()
            successful_download = True

        #handle errors
        except HTTPError, e:
            print "HTTP Error:",e.code , url
            print "Trying again: %i attempts remaining" % (trys_left+1)
            if trys_left <= -1: qErr.qErr()
        except URLError, e:
            print "URL Error:",e.reason , url
            print "Trying again: %i attempts remaining" % (trys_left+1)
            if trys_left <= -1: qErr.qErr()
        except:
            print "Couldn't Download!"
            print "Trying again: %i attempts remaining" % (trys_left+1)
            if trys_left <= -1: qErr.qErr()

def downloadImage(img_url,out_name=None, timeout=20):
    
    #create file name based on known pattern
    # Now download the image. If these were text files,
    # or other ascii types, just pass an empty string
    # for the second param ala stealStuff(file_name,'',base_url)
    if not out_name:
        try:
            out_name = img_url.split('/')[-1]
        except:
            out_name = 'qImage.jpg'
    stealStuff(out_name,"b",img_url, timeout=timeout)

def MakeFindingChart(ra=198.40130,dec=8.09730,uncertainty=1.8,src_name='GRB090313',pos_label='XRT',survey='dss2red',cont_str='Contact: Test', size=3.0,err_shape='cross',incl_scale=True,return_svg=False):
    '''if return_svg, actually return the svg text rather than saving it to file.  Used for the online finding chart generator.'''
    fc = qImage()
    # define pixel scale from side size
    try:
        uncertainty = float(uncertainty)
        ra = float(ra)
        dec = float(dec)
    except:
        try:
            float_ra_dec = sex2dec(ra,dec)
            uncertainty = float(uncertainty)
            ra=float(float_ra_dec[0])
            dec=float(float_ra_dec[1])
        except:
            raise ValueError('RA/Dec or uncertainty misformatted.')
    try:
       size = float(size)
    except:
       pass  #assume we wanted auto

    if size:
        try:
            if size.upper() == 'AUTO' and uncertainty > 0:
                # Auto sets error circle size to ~10 pixels
                if not uncertainty < 1.0: 
                    side_size_arcmin = round(uncertainty/10.0/60.0*600.0)
                else:
                    side_size_arcmin = 1.0
            else:
                print 'Invalid string for size; defaulting to 3 arcmin.'
                side_size_arcmin = 3.0
        except(AttributeError):
            try:
                side_size_arcmin = float(size)
            except:
                print 'Invalid entry for size; defaulting to 3 arcmin'
                side_size_arcmin = 3.0
    # There's a limit to how big an image you can request; set it here to 20'
    if side_size_arcmin > 20.0:
        side_size_arcmin = 20.0
    side_size = side_size_arcmin * 60.0 # arcseconds
    img_size = 600 #pixels
    pixel_scale = side_size/img_size #arcsec/pixel
    
    if survey != 'sdss':
        fc.dss_grab(ra,dec,side_size_arcmin,survey)
        fc.invert_and_resize(img_size)
    else:
        fc.sdss_grab(ra,dec,pixel_scale,img_size)
        fc.invert_and_resize(img_size)
    
    if str(incl_scale).lower()=='false' or str(incl_scale).lower()=='no':
        suppress_scale=True
    else:
        suppress_scale=False

    if not return_svg:
        svgout = fc.overlay_finding_chart(ra,dec,uncertainty,src_name,pos_label,cont_str,err_shape,suppress_scale)
        outname = storepath+src_name+'_fc'
        outnamesvg = outname+'.svg'
        outnamepng = outname+'.png'
        f=open(outnamesvg,'w')
        f.write(svgout)
        f.close()
        magickcommand = "convert %s %s" % (outnamesvg, outnamepng)
        try:
            os.system(magickcommand)
        except:
            print "Do not have Imagemagick, cannot convert svg to png"
        print storepath+outname+'*'
        outlist = glob.glob(outname+'*')
        return outlist
    
    if return_svg:
        return fc.overlay_finding_chart(ra,dec,uncertainty,src_name,pos_label,cont_str,err_shape,suppress_scale)
