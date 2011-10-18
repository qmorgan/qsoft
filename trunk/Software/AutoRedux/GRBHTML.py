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
from MiscBin import qErr
from Phot import extinction
from AutoRedux import qHTML

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
        
        self.qHTML = qHTML.qHTML(triggerid,base_dir)
        self.triggerid = triggerid
        self.qHTML.create_folder()
        title = "Swift Trigger %s" % str(self.triggerid)
        self.qHTML.create_header(title=title)
        self.successful_export = False
    
    def create_folder(self):
        self.out_dir = self.base_dir + '/' + str(self.triggerid)
        self.out_dir_name = os.path.basename(self.out_dir)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        
    
    def add_timing_info(self,grb_time=None):
        '''Requires burst time as string, t90 as float or string'''
        if grb_time:
            content = "Burst Time: " + str(grb_time)
            self.qHTML.add_post(title='Temporal Information',content=content)
        
    def add_position_info(self,bat_pos=None,xrt_pos=None,uvot_pos=None,reg_path=None):
        '''Requires bat_pos, xrt_pos, uvot_pos in format of (ra,dec,uncertainty), 
        where uncertainty is in arcseconds. 
        
        '''
        content = ""
        if bat_pos: 
            bat_pos_sex = q.dec2sex((bat_pos[0],bat_pos[1]))
            self.best_pos = bat_pos
            self.best_pos_type = 'BAT'
        if xrt_pos: 
            xrt_pos_sex = q.dec2sex((xrt_pos[0],xrt_pos[1]))
            self.best_pos = xrt_pos
            self.best_pos_type = 'XRT'
        if uvot_pos: 
            uvot_pos_sex = q.dec2sex((uvot_pos[0],uvot_pos[1]))
            self.best_pos = uvot_pos
            self.best_pos_type = 'UVOT'
        
        if bat_pos != None:
            content += '''
            <b>BAT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
            ''' % (str(bat_pos[0]).rstrip('0'),str(bat_pos[1]).rstrip('0'),str(bat_pos[2]).rstrip('0'),\
                   bat_pos_sex[0],bat_pos_sex[1])
        if xrt_pos != None:
            content += '''
            <b>XRT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
            ''' % (str(xrt_pos[0]).rstrip('0'),str(xrt_pos[1]).rstrip('0'),str(xrt_pos[2]).rstrip('0'),\
            xrt_pos_sex[0],xrt_pos_sex[1])
        if uvot_pos != None:
            content += '''
            <b>UVOT</b>: RA = %s, Dec = %s, Uncertainty: %s"<br>
            (RA = %s, Dec = %s)<br><br>
            ''' % (str(uvot_pos[0]).rstrip('0'),str(uvot_pos[1]).rstrip('0'),str(uvot_pos[2]).rstrip('0'),\
            uvot_pos_sex[0],uvot_pos_sex[1])
            
        if self.best_pos:
            try:
                source_name = 'swift_' + str(self.triggerid)
                Gal_EB_V = extinction.qExtinction(source_name,self.best_pos[0],self.best_pos[1])
            except:
                Gal_EB_V = 'Unknown'
            content += '''
            <b>Extinction</b>: E_(B-V) = %s<br><br>
            ''' % (str(Gal_EB_V))
        
        if reg_path != None:
            if os.path.exists(reg_path):
                self.qHTML.copy_file(reg_path)
                content += '''
                <a href='./%s'>DS9 Region File</a>
        
                ''' % (os.path.basename(reg_path))
            else:
                print reg_path + ' does not exist.  Not including region file.'
            self.reg_path = reg_path
            self.reg_name = os.path.basename(reg_path)
        if content:
            self.qHTML.add_post(title='Spatial Information',content=content)
            
    
    def add_finder_chart_info(self,fc_path):
        content = ''
        if fc_path:
            if os.path.exists(fc_path):
                fc_base = os.path.basename(fc_path)
                self.qHTML.copy_file(fc_path)
                content += '''
                <a href='http://fc.qmorgan.com/fcserver.py?ra=%f&dec=%f&uncertainty=%f&err_shape=combo&incl_scale=yes&size=AUTO&src_name=Swift_%s&pos_label=%s&cont_str=&survey=sdss'>SDSS Finding Chart</a> (May not be available)<br>
                DSS Finding Chart:<br>
                <a href='./%s'><img src='%s' alt="Finder Chart", title="Finder Chart",width=500></a>
                ''' % (self.best_pos[0],self.best_pos[1],self.best_pos[2],self.triggerid,self.best_pos_type,fc_base,fc_base)
            else:
                print fc_path + ' does not exist.  Not including Finder Chart.'
            self.fc_path = fc_path
            self.fc_name = os.path.basename(fc_path)
        if content:
            self.qHTML.add_post(title='Finding Charts',content=content)
            
    
    def add_telescope_info(self):
        too_str = ''
        content = ''
        if hasattr(self, 'best_pos'):
            too_str = '?ra=%f&dec=%f' % (self.best_pos[0],self.best_pos[1])
            react_str = 'http://www.srl.caltech.edu/~react/Swift%s.html' % (self.triggerid)
            datascope_str = '?position=%f+%f&size=0.05' % (self.best_pos[0],self.best_pos[1])            
            content += '''
                <a href='http://lyra.berkeley.edu/GRB/too/too.php%s'>GRAASP ToO Page</a><br>
                <a href='%s'>Caltech React Page</a><br>
                <a href='http://heasarc.gsfc.nasa.gov/cgi-bin/vo/datascope/jds.pl%s'>NVO DataScope Query</a><br>
                ''' % (too_str,react_str,datascope_str)
        content += '''
        <a href='http://gcn.gsfc.nasa.gov/other/%s.swift'>GCN Notices for this Trigger</a><br>
        ''' % (str(self.triggerid))
        if content:
            self.qHTML.add_post(title='Useful Links',content=content)
            
    
    def add_footer(self):
        update_time = time.ctime(time.time())
        footercontent = '''
        This page is updated automatically as more information arrives.<br>
        Last Updated: %s <P>
        <ADDRESS> Adam N. Morgan (qmorgan@gmail.com)</ADDRESS>
        </center>
        </html>
        ''' % (update_time)
        self.qHTML.create_footer(footercontent)
        

def MakeGRBPage(html_path='/home/amorgan/www/swift',triggerid='000000',\
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
    inst.qHTML.export_html()
    return inst

def MakeGRBIndex(collected_grb_dict,html_path='/home/amorgan/www/swift'):
    '''Takes a collected dictionary from CollectGRBInfo and creates an index
    html page for all the GRBs
    '''
    failed_grbs = []
    incl_files = ['reg_path','fc_path']
    incl_keys = ['triggerid_str','z','Q_hat']
    html_block = '''
    <html><head><title>Q's GRB Pages</title></head>
    <body background="http://static.tumblr.com/snnreod/Fx4l8ig9j/background_dark.jpg" bgcolor="#363636" text="#000000" link="#111111" alink="#111111" vlink="#333333">
    <center><font size="+2" color="#FFFFFF">Q's GRB Pages</font><p></center>
    
    <style type="text/css">
    table.sample {
    border-width: 1px;
    border-spacing: ;
    border-style: outset;
    border-color: white;
    border-collapse: collapse;
    background-color: rgb(255, 255, 240);
    }
    table.sample th {
    border-width: 3px;
    padding: 3px;
    border-style: solid;
    border-color: black;
    background-color: white;
    -moz-border-radius: ;
    }
    table.sample td {
    border-width: 3px;
    padding: 3px;
    border-style: solid;
    border-color: black;
    background-color: white;
    -moz-border-radius: ;
    }

    </style>
    
    <table class="sample" align="center">
    ''' 
    table_columns = ('GRB','Region File','Finding Chart','Trigger #','Redshift','Q_hat')
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
                    html_block += "<td></td>"
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
    bat_pos=(132.45,3.52,22)
    xrt_pos=(132.452,3.623,2.6)
    uvot_pos=(132.4526,3.6234,1.1)
    reg_path='/Users/amorgan/Desktop/123456.reg'
    grb_time='09/09/27 10:07:16.92'
    fc_path = '/Users/amorgan/Desktop/fcserver.png'
    
    MakeGRBPage(html_path='/Users/amorgan/Desktop/',triggerid=triggerid,bat_pos=bat_pos,xrt_pos=xrt_pos,uvot_pos=uvot_pos,
                reg_path=reg_path, grb_time=grb_time,fc_path=fc_path)
    # inst = GRBHTML('123416','/Users/amorgan/Public/TestDir')
    #     inst.add_position_info(bat_pos=(13.45,-32.52,22),xrt_pos=(13.452,-32.622,2.5),reg_path='/Users/amorgan/Desktop/123456.reg')
    #     inst.add_timing_info(grb_time='09/09/27 10:07:16.92')
    #     inst.add_telescope_info()
    #     inst.add_finder_chart_info('/Users/amorgan/Desktop/test.png')
    #     inst.add_footer()
    #     inst.export_html()
    
