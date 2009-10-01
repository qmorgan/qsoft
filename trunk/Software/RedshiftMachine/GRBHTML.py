#!/usr/bin/env python
# encoding: utf-8
"""
CreateSimpleHTML.py
Author: Adam N. Morgan
Created: Sept 22, 2009
	
This Program creates a very simple GRB HTML Document. Moves
everything to the correct directory.

"""
import sys
import os
import shutil
import time
from MiscBin import q

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

class GRBHTML:
    '''Creates a tiny block of HTML to put within the larger page.  
    
    '''
    def __init__(self,triggerid,base_dir):
        self.base_dir = base_dir
        if not os.path.exists(base_dir):
            print "Output Directory does not exist. Exiting"
            sys.exit(1)
        self.triggerid = triggerid
        self.create_folder()
        self.add_header()
    
    def create_folder(self):
        self.out_dir = self.base_dir + '/' + str(self.triggerid)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
    
    def copy_file(self, file_path):

        if os.path.exists(file_path):
            shutil.copy(file_path,self.out_dir)
        else:
            print "File %s does not exist" % (file_path)
                
    
    def add_header(self):
        self.html_block = '''
        <html><head><title>Swift Trigger %i</title></head>
        <body bgcolor="#FFFFFF" text="#000066">
        <center>
        <font size="+2">Swift Trigger %i</font><p>
        
        ''' % (self.triggerid, self.triggerid)
    
    def add_position_info(self,bat_pos=(123.45,54.32,30),xrt_pos=None,uvot_pos=None,reg_path=None):
        '''Requires bat_pos, xrt_pos, uvot_pos in format of (ra,dec,uncertainty), 
        where uncertainty is in arcseconds. 
        
        '''
        self.html_block += '''
        <hr width="50%">
        <b>Position Information:</b><p>
        '''
        if bat_pos: 
            bat_pos_sex = q.dec2sex((bat_pos[0],bat_pos[1]))
            self.best_pos = bat_pos
        if xrt_pos: 
            xrt_pos_sex = q.dec2sex((xrt_pos[0],xrt_pos[1]))
            self.best_pos = xrt_pos
        if uvot_pos: 
            uvot_pos_sex = q.dec2sex((uvot_pos[0],uvot_pos[1]))
            self.best_pos = uvot_pos
        
        if bat_pos != None:
            self.html_block += '''
            <b>BAT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
            ''' % (str(bat_pos[0]).rstrip('0'),str(bat_pos[1]).rstrip('0'),str(bat_pos[2]).rstrip('0'),\
                   bat_pos_sex[0],bat_pos_sex[1])
        if xrt_pos != None:
            self.html_block += '''
            <b>XRT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
            ''' % (str(xrt_pos[0]).rstrip('0'),str(xrt_pos[1]).rstrip('0'),str(xrt_pos[2]).rstrip('0'),\
            xrt_pos_sex[0],xrt_pos_sex[1])
        if uvot_pos != None:
            self.html_block += '''
            <b>XRT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br><br>
            (RA = %s, Dec = %s)<br>
            ''' % (str(uvot_pos[0]).rstrip('0'),str(uvot_pos[1]).rstrip('0'),str(uvot_pos[2]).rstrip('0'),\
            uvot_pos_sex[0],uvot_pos_sex[1])
        if reg_path != None:
            if os.path.exists(reg_path):
                self.copy_file(reg_path)
                self.html_block += '''
                <a href='./%s'>DS9 Region File</a>
        
                ''' % (os.path.basename(reg_path))
            else:
                print reg_path + ' does not exist.  Not including region file.'
    
    def add_finder_chart_info(self,fc_path):
        if os.path.exists(fc_path):
            fc_base = os.path.basename(fc_path)
            self.copy_file(fc_path)
            self.html_block += '''
            <hr width="50%%">
            <b>Finder Chart:</b><p>
            <a href='./%s'><img src='%s' alt="Finder Chart", title="Finder Chart",width=200, height=200></a>
        
            ''' % (fc_base,fc_base)
        else:
            print reg_path + ' does not exist.  Not including Finder Chart.'
    
    def add_telescope_info(self):
        too_str = ''
        if hasattr(self, 'best_pos'):
            too_str = '?ra=%f&dec=%f' % (self.best_pos[0],self.best_pos[1])
        self.html_block += '''
            <hr width="50%%">
            <b>Telescope ToO Information:</b><p>
            <a href='http://lyra.berkeley.edu/GRB/too/too.php%s'>GRAASP ToO Page</a>
        ''' % (too_str)
    
    def add_footer(self):
        update_time = time.ctime(time.time())
        self.html_block += '''
        <HR WIDTH="50%%">
        This page is updated automatically as more information arrives.<br>
        Last Updated: %s <P>
        <ADDRESS> Adam N. Morgan (qmorgan@gmail.com)</ADDRESS>
        </center>
        </html>
        ''' % (update_time)
        
    def export_html(self):
        '''
        Export the html_block string into the created folder.
        '''
        filename = self.out_dir + '/index.html'
        f = open(filename,'w')
        f.write(self.html_block)
        f.close()
    

def test():
    inst = GRBHTML(123416,'/Users/amorgan/Public/TestDir')
    inst.add_position_info(bat_pos=(13.45,-32.52,22),xrt_pos=(13.452,-32.622,2.5),reg_path='/Users/amorgan/Desktop/123456.reg')
    inst.add_telescope_info()
    inst.add_finder_chart_info('/Users/amorgan/Desktop/test.png')
    inst.add_footer()
    inst.export_html()

