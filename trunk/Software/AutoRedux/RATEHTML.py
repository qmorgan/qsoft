#!/usr/bin/env python
# encoding: utf-8
"""
RATEHTML.py
Author: Adam N. Morgan
Created: Jan 18 2012
	
A wrapper around qHTML to make RATE webpages

"""
import sys
import os
import shutil
import time
from MiscBin import q
from MiscBin import qErr
from AutoRedux import qHTML
from AutoRedux import GRBHTML
from RedshiftMachine import LoadDB

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

class RATEHTML:
    def __init__(self,basedir="/home/amorgan/www/",
                linklist=[('About RATE GRB-z',"./index.html"),
                        ('How it Works',"./howitworks.html"),
                        ('The Paper',"./paper.html"),
                        ('Contact Us',"./ratecontact.html"),
                        ('Swift GRB Pages',"http.//swift.qmorgan.com"),
                        ('RSS Feed',"http://swift.qmorgan.com/rss.xml"),
                        ('Finding Chart Generator',"http://fc.qmorgan.com"),
                        ("A. N. Morgan's Webpage","http://qmorgan.com/")]):
        '''Initialize a RATEHTML object'''
        self.linklist = linklist
        self.basedir = basedir
        
    def RemakeIndex(self):
        title = "RATE GRB-z"
        sidetitle = "Random forest Automated Triage Estimator for GRB Redshifts"
        logopath = "rategrbzlogoc_300.png"
        inst = qHTML.qHTML("rategrbz",self.basedir,create_folder=True,linklist=self.linklist)
        inst.create_header(title=title)
        inst.create_sidebar(logo_url=logopath)
        inst.create_footer()
        introcontent = '''
        <p>	As the discovery rate of potentially interesting transients dwarfs the 
        	available resources with which to perform follow-up observations, it 
        	becomes necessary to have tools to inform whether or not to 
        	spend resources on each new event. 	
        </p>
        <p>	We have developed the Random forest Automated Triage Estimator for GRB redshifts 
        	(RATE GRB-z) to assist in rapidly identifying 
        	high-redshift (<i>z</i> > 4) candidates using early-time metrics 
        	from the three telescopes onboard <i>Swift</i>.
        </p>
        <p> <strong>The method reduces this decision of whether or not to follow-up 
        	a new GRB to a single parameter: <i>Q</i>_hat. If a user wishes to devote 20% of 
        	telescope observing time to chasing high-<i>z</i> GRB candidates, then we 
        	recommend follow-up for all events with <i>Q</i>_hat < 0.2. </p></strong>

        <p> See <a class="page_link" href="./howitworks.html"><span>How it Works</span></a> for more
        	information on the methodology and expected performance metrics, or read 
        	<a class="page_link" href="http://arxiv.org/abs/1112.3654"><span>the paper</span></a>
            for the full details. The latest available predictions are shown below, and real-time
        	alerts are available by subscribing to our <a class="page_link" href="http://swift.qmorgan.com/rss.xml"><span>RSS Feed</span></a>.
        	We also provide our data products for all past <a class="page_link" href="http://swift.qmorgan.com"><i>Swift</i> GRBs</a>.
        	We welcome questions/comments/suggestions: please <a class="page_link" href="./ratecontact.html"><span>contact us</span></a>.
        </p>
        '''
        inst.add_post(title="Introduction",content=introcontent)
        
        db = LoadDB.LoadDB("GRB_full")
        GRB_table_html = GRBHTML.MakeGRBTable(db.dict,maxlength=15) 
        inst.add_post(title="Predictions for Latest GRBs",content = GRB_table_html)
        
        update_time = time.ctime(time.time())
        footer_content = '''
        This page is updated automatically as more information arrives.<br>
        Last Updated: %s <P>
        <ADDRESS> Adam N. Morgan (qmorgan@gmail.com)</ADDRESS>
        </center>
        </html>
        ''' % (update_time)
        inst.create_footer(footer_content)
        
        
        inst.export_html()