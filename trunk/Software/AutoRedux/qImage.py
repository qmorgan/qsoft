import sys
import os
try:
    from PIL import Image
    from PIL import ImageOps
except:
    sys.exit('You do not have PIL.  Download it.')
import base64, Image, string
import glob
from MiscBin import q

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
        self.wcs_size=size
        
        self.image_path=storepath+'dss_grab.gif'
        
        downloadImage(url_str,self.image_path)
        
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

        im1 = Image.open(self.image_path)
        # Invert Colors
        im1 = ImageOps.invert(im1)
        # Convert to color mode
        im1 = im1.convert('RGB')
        # Resize with a cubic spline - it looks nice.
        # Other options are NEAREST, BILINEAR, ANTIALIAS
        im1 = im1.resize((newsize,newsize),Image.BICUBIC)
        self.save_str = storepath+'FCBase.jpg'
        im1.save(self.save_str)

            
    def overlay_finding_chart(self,ra,dec,uncertainty,src_name='unknown source',pos_label='UVOT',cont_str='Adam N. Morgan (qmorgan@gmail.com, 510-229-7683)'):
        
        self.loc_ra = ra
        self.loc_dec = dec
        self.loc_uncertainty=uncertainty
        self.src_name = src_name
        self.pos_label = pos_label
        self.cont_str = cont_str
        
        self.pixel_scale = self.wcs_size * 60.0 / self.img_size # arcsec/pixel
        
        self.str_loc_pos = q.dec2sex((self.loc_ra,self.loc_dec))
        
        head_buffer = 50
        uncertainty_circle_size = int(self.loc_uncertainty/self.pixel_scale)
        pos_pos_x = 300 #position of the error circle on the image , hard coded for now
        pos_pos_y = 350
        
        scale_line_length = int(60/self.pixel_scale) #1 arcmin
        scale_line_start = head_buffer + self.img_size/2 - scale_line_length/2
        scale_line_stop = head_buffer + self.img_size/2 + scale_line_length/2
        scale_line_label_pos_y = scale_line_start - 15
        
        pos_label_pos_y = pos_pos_y - uncertainty_circle_size - 15 
        
        epoch = 'J2000'
        pos_label_full = self.pos_label + ' Pos. (%s)' % epoch        
               
        if self.survey == 'poss2ukstu_red': 
            self.survey_str = 'DSS2 (red)'
        elif survey == 'poss1_red':
            self.survey_str = 'DSS1 (red)'
        else: 
            print 'Unknown survey: %s' % survey
            sys.exit(1)

        ra_str = 'RA = ' + self.str_loc_pos[0]
        dec_str = 'Dec = ' + self.str_loc_pos[1]
        uncertainty_str = 'Uncertainty: ' + str(self.loc_uncertainty) + '"'
        
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
         <tspan x="15" dy="1em">Contact: %s</tspan>
         </text>

         <text x ="288" y="%d" fill="blue" font-size="12">
         <tspan x="288" dy="1em">%s</tspan>
         </text>
         <ellipse cx="%d" cy="%d" rx="%d" ry="%d"
         style="fill:none;
         stroke:rgb(0,0,255);stroke-width:2"/>

         <text x="615" y ="55" fill = "black" font-size = "18">
         <tspan x="615" dy="1em">%s</tspan>
         <tspan x="615" dy="1em">RA = %s</tspan>
         <tspan x="615" dy="1em">Dec = %s</tspan>
         <tspan x="615" dy="1em">Uncertainty = %s"</tspan>
         </text>

         <text x="640" y="325" fill="black" font-size="18">
         <tspan x="640" dy="1em">Image from:</tspan>
         <tspan x="640" dy="1em">%s</tspan>
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


         <line x1="200" y1="%d" x2="200" y2="%d"
         style="stroke:rgb(255,0,0);stroke-width:2"/>
         <text x ="177" y ="%d" fill ="red" font-size = "12">
         <tspan x="177" dy="1em">1 arcmin</tspan>
         </text>
        
        </svg>''' % ( fc_width, fc_height, head_buffer, im.size[0], im.size[1], im.format.lower(), imstr,\
                self.src_name,  self.cont_str, pos_label_pos_y, self.pos_label, pos_pos_x, pos_pos_y, uncertainty_circle_size,uncertainty_circle_size,\
                 pos_label_full, self.str_loc_pos[0], self.str_loc_pos[1],\
                 str(self.loc_uncertainty),self.survey_str, scale_line_start,scale_line_stop,scale_line_label_pos_y)
    


def stealStuff(file_name,file_mode,base_url):
	from urllib2 import Request, urlopen, URLError, HTTPError
    
	#create the url and the request
	url = base_url + file_name
	req = Request(url)

	# Open the url
	try:
		f = urlopen(req)
		print "downloading " + url

		# Open our local file for writing
		local_file = open(file_name, "w" + file_mode)
		#Write to our local file
		local_file.write(f.read())
		local_file.close()

	#handle errors
	except HTTPError, e:
		print "HTTP Error:",e.code , url
	except URLError, e:
		print "URL Error:",e.reason , url


def downloadImage(img_url,out_name='qImage.jpg'):
    
	#create file name based on known pattern
	# Now download the image. If these were text files,
	# or other ascii types, just pass an empty string
	# for the second param ala stealStuff(file_name,'',base_url)
	stealStuff(out_name,"b",img_url)


def MakeFindingChart(ra=198.40130,dec=8.09730,uncertainty=1.8,src_name='GRB090313',pos_label='XRT',survey='dss2red'):
    fc = qImage()
    fc.dss_grab(ra,dec,3.0,survey)
    fc.invert_and_resize(600)
    svgout = fc.overlay_finding_chart(ra,dec,uncertainty,src_name,pos_label)
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