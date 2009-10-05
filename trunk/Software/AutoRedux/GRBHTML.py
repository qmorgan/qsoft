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
            print "Output Directory %s does not exist. Exiting" % (base_dir)
            sys.exit(1)
        self.triggerid = triggerid
        self.create_folder()
        self.add_header()
    
    def create_folder(self):
        self.out_dir = self.base_dir + '/' + str(self.triggerid)
        self.out_dir_name = os.path.basename(self.out_dir)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
    
    def copy_file(self, file_path):

        if os.path.exists(file_path):
            shutil.copy(file_path,self.out_dir)
            newpath = self.out_dir + '/' + os.path.basename(file_path)
            return newpath
        else:
            print "File %s does not exist" % (file_path)
                
    
    def add_header(self):
        self.html_block = '''
        <html><head><title>Swift Trigger %s</title></head>
        <body bgcolor="#FFFFFF" text="#000066">
        <center>
        <font size="+2">Swift Trigger %s</font><p>
        
        ''' % (self.triggerid, self.triggerid)
    
    def add_timing_info(self,grb_time=None):
        '''Requires burst time as string, t90 as float or string'''
        if grb_time:
            self.html_block += '''
            <hr width="50%%">
            <b>Temporal Information:</b><p>
            <b>Burst Time:</b> %s
            
            ''' % (grb_time)
    
    def add_position_info(self,bat_pos=None,xrt_pos=None,uvot_pos=None,reg_path=None):
        '''Requires bat_pos, xrt_pos, uvot_pos in format of (ra,dec,uncertainty), 
        where uncertainty is in arcseconds. 
        
        '''
        self.html_block += '''
        <hr width="50%%">
        <b>Spatial Information:</b><p>
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
            <b>UVOT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
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
            self.reg_path = reg_path
            self.reg_name = os.path.basename(reg_path)
    
    def add_finder_chart_info(self,fc_path):
        if fc_path:
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
            self.fc_path = fc_path
            self.fc_name = os.path.basename(fc_path)
    
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
    


def MakeGRBPage(html_path='/Users/amorgan/Public/TestDir',triggerid='000000',\
                bat_pos=None,xrt_pos=None,uvot_pos=None,reg_path=None,\
                grb_time=None,fc_path=None):
    '''Make a GRB page given inputs and return the instance of the html.'''
    triggerid = str(triggerid)
    inst = GRBHTML(triggerid,html_path)
    inst.add_timing_info(grb_time)
    inst.add_position_info(bat_pos,xrt_pos,uvot_pos,reg_path)
    inst.add_telescope_info()
    inst.add_finder_chart_info(fc_path)
    inst.add_footer()
    inst.export_html()
    return inst

def MakeGRBIndex(collected_grb_dict,html_path='/Users/amorgan/Public/TestDir'):
    '''Takes a collected dictionary from CollectGRBInfo and creates an index
    html page for all the GRBs
    '''
    failed_grbs = []
    incl_files = ['reg_path','fc_path']
    incl_keys = ['triggerid_str','z']
    html_block = '''
    <html><head><title>List of GRB Pages</title></head>
    <body bgcolor="#FFFFFF" text="#000066">
    <font size="+2">List of GRB Pages</font><p>
    BETA!<p>
    
    <table border="1">
    ''' 
    table_columns = ('GRB','Region File','Finding Chart','Trigger #','Redshift')
    table_label = '<tr>'
    for column_name in table_columns:
        table_label += '<td align=center><b>%s</b></td>' % column_name
    table_label += '<tr>'
    
    html_block += table_label
    
    sortedkeys=collected_grb_dict.keys()
    sortedkeys.sort()
    
    table_entry_count=0
    
    for grb in reversed(sortedkeys):
        # Grab the folder name of the succesful Web Page creations
        grbdict = collected_grb_dict[grb]
        try:
            html_block += "<tr>"
            try:
                grbfolder = os.path.basename(grbdict['out_dir'])
                html_block += "<td><a href='./%s'>%s</a><br></td>" % (grbfolder,grb)
            except:
                grbfolder = None
                html_block += "<td>%s</td>" % (grb)
            for incl_file in incl_files:
                try:
                    file_path = grbfolder + '/' + os.path.basename(grbdict[incl_file])
                    html_block += "<td align=center><a href='./%s'>Link</td>" % file_path
                except:
                    file_path = None
                    html_block += "<td></td>"
            for incl_item in incl_keys:
                try:
                    html_block += "<td>%s</td>" % str(grbdict[incl_item])
                except:
                    pass
#            html_block += "<td>%s</td>" % (grbdict['triggerid_str'])
#            html_block += "<td>%f</td>" % (grbdict['z'])
            html_block += "</tr>"
            if table_entry_count == 20:
                html_block += table_label
                table_entry_count = 0
            table_entry_count += 1
        except:
            failed_grbs.append(grb)
    
    update_time = time.ctime(time.time())
    
    print failed_grbs
    
    html_block += '''
    </table>
    
    <HR WIDTH="50%%">
    This page is updated automatically as more information arrives.<br>
    Last Updated: %s <P>
    <ADDRESS> Adam N. Morgan (qmorgan@gmail.com)</ADDRESS>
    </center>
    </html>
    ''' % (update_time)

    filename = html_path + '/index.html'
    f = open(filename,'w')
    f.write(html_block)
    f.close()

def SortGRBDict(collected_grb_dict):
    keys = collected_grb_dict.keys()
    keys.sort()
    sorted_vals = map(collected_grb_dict.get,keys)
    sorted_tuple = (keys,sorted_vals)
    return sorted_tuple

def test():
    triggerid='543210'
    bat_pos=(13.45,-32.52,22)
    xrt_pos=(13.452,-32.623,2.6)
    uvot_pos=(13.4526,-32.6234,1.1)
    reg_path='/Users/amorgan/Desktop/123456.reg'
    grb_time='09/09/27 10:07:16.92'
    fc_path = '/Users/amorgan/Desktop/test.png'
    
    MakeGRBPage(triggerid=triggerid,bat_pos=bat_pos,xrt_pos=xrt_pos,uvot_pos=uvot_pos,
                reg_path=reg_path, grb_time=grb_time,fc_path=fc_path)
    # inst = GRBHTML('123416','/Users/amorgan/Public/TestDir')
    #     inst.add_position_info(bat_pos=(13.45,-32.52,22),xrt_pos=(13.452,-32.622,2.5),reg_path='/Users/amorgan/Desktop/123456.reg')
    #     inst.add_timing_info(grb_time='09/09/27 10:07:16.92')
    #     inst.add_telescope_info()
    #     inst.add_finder_chart_info('/Users/amorgan/Desktop/test.png')
    #     inst.add_footer()
    #     inst.export_html()
    
