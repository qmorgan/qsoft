#!/usr/bin/env python
# encoding: utf-8
"""
qHTML.py
Author: Adam N. Morgan
Created: Oct 15, 2011
    
Creates an HTML document based on a tumblr post.

"""
import sys
import os
import shutil
import time
from MiscBin import q
from MiscBin import qErr
import glob

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'


class qHTML:
    '''The base block for a q website. The basic form for this website is a 
    folder of a certain name with an index.html file as the main window. 
    Included in the folder are relavant images, files, and the pickle file for 
    this qHTML object.
    '''
    
    def __init__(self,name,base_dir,create_folder=True,logo=None):
        '''If not create_folder, then dump all the files in the current base directory.'''
        self.base_dir = base_dir
        if not os.path.exists(base_dir):
            print "Output Directory %s does not exist. Exiting" % (base_dir)
            sys.exit(1)
        self.name = name
        self.posts = []
        self.plainhtmls = []
        if create_folder:
            self.create_folder()
        else:
            self.out_dir = self.base_dir
            self.out_dir_name = os.path.basename(self.out_dir)
        
        # Check to make sure we can find the relevant css files. If not, put them there.
        css_path = loadpath + '*.css'
        check_path1 = self.out_dir + '/*.css'
        check_path2 = self.base_dir + '/*.css'
        if not glob.glob(check_path1) and not glob.glob(check_path2):
            self.copy_file(css_path)
        
        self.create_header()
        self.create_footer()
        self.create_sidebar()
        self.successful_export = False

    def create_folder(self):
        self.out_dir = self.base_dir + '/' + str(self.name)
        self.out_dir_name = os.path.basename(self.out_dir)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
    
    def copy_file(self, file_path):
        try:
            cmd = "cp %s %s" % (file_path,self.out_dir)
            # having problems with shutil.copy
            # shutil.copy(file_path,self.out_dir)
            os.system(cmd)
            newpath = self.out_dir + '/' + os.path.basename(file_path)
            return newpath
        except:
            if not os.path.exists(file_path):
                errmesg = "File %s does not exist. Could not copy to %s" % (file_path,self.out_dir)
                qErr.qErr(errtitle=errmesg)
            else:
                errmesg = "Could not copy file %s to %s" % (file_path,self.out_dir)
                qErr.qErr(errtitle=errmesg)
            
    def create_header(self,title=''):
        self.header = '''
        <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
            "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
        <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
            <head>
                <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
                <title>%s</title>
                <link rel="icon" href="http://i.imgur.com/wjKnI.png"/>
                <link rel="alternate" type="application/rss+xml" title="RSS" href="http://qmorgan.com/rss"/>
                <meta name="description" content="Random forest Automated Triage Estimator for GRB Redshifts" />
                <meta name="viewport" content="width=775"/> <!-- iPhone -->
                <meta name="image:Background" content=""/>
                <meta name="image:Header" content=""/>
                <meta name="text:Google Analytics ID" content="" />          
                <meta name="text:Disqus Shortname" content="" />     
                <meta name="if:Show My Portrait" content="" />
                <meta name="if:Huge Photos" content="" />
                <meta name="if:Dark Layout" content="" />
                <link rel="apple-touch-icon" href="http://i.imgur.com/rWYK2.png"/> 
                <link rel="stylesheet" href="http://static.tumblr.com/snnreod/yt6l8jer7/screen.css" type="text/css" media="screen" charset="utf-8" />  
                <link rel="stylesheet" href="../screen_backup.css" type="text/css" media="screen" charset="utf-8" />  
                <link rel="stylesheet" href="screen_backup.css" type="text/css" media="screen" charset="utf-8" />  
                <link rel="stylesheet" href="../table.css" type="table/css" media="screen" charset="utf-8" />
                <link rel="stylesheet" href="table.css" type="table/css" media="screen" charset="utf-8" />
                
                <style type="text/css">

                  html, body {
                    background: #d1d2d1 url(http://static.tumblr.com/snnreod/DTFl8fjqb/background.jpg) top left repeat;
                  }  

                     /* Dark Layout aka SIMPLE BRUTALITY! */
                      html, body {
                      background: #363636 url(http://static.tumblr.com/snnreod/Fx4l8ig9j/background_dark.jpg) top left repeat;
                      }

                      #header h1 a.no_header {
                        color: #fff;
                        text-shadow: 0 2px 4px #000;

                      }

                      .post {
                        -moz-box-shadow:    0 1px 8px rgba(0,0,0,.9);
                        -webkit-box-shadow: 0 1px 8px rgba(0,0,0,.9);
                        box-shadow:         0 1px 8px rgba(0,0,0,.9);
                      }

                      #sidebar #description {
                        color: #e1e1e1;
                        text-shadow: 0 1px 1px #000;
                        text-shadow: 0 1px 1px rgba(0,0,0,.9) !important;
                        font-weight: bold !important;
                      }

                      #sidebar #description a {
                        text-shadow: 0 1px 1px #000;
                        text-shadow: 0 1px 1px rgba(0,0,0,.9);
                        color: #e1e1e1;
                      }

                      #sidebar #description a:hover {
                        color: #fff;
                      }

                      #sidebar  .page_link, #sidebar #search label, #sidebar #description a.page_link {
                        color: #e1e1e1;
                        text-shadow: 0 1px 1px #000;
                        font-weight: normal;
                        border-top: 1px solid #555555;
                        *border-top:1px solid #555555;
                        border-top: 1px solid rgba(255,255,255, .15);
                      }

                      #sidebar  .page_link:hover, #sidebar #description a.page_link:hover {
                        color: #fff;
                      }

                      #sidebar #search form #search_submit {
                        -moz-box-shadow:    0 1px 4px rgba(0,0,0,.7);
                        -webkit-box-shadow: 0 1px 4px rgba(0,0,0,.7);
                        box-shadow:         0 1px 4px rgba(0,0,0,.7); 
                        border: none;
                      }

                      #newer_posts, #older_posts {
                        color: #e1e1e1;
                        text-shadow: 0 1px 1px #000;
                        font-weight: normal;
                      }

                      #newer_posts:hover, #older_posts:hover {
                        color: #fff;
                      }

                      .page_of {
                        display: inline-block;
                        width: auto;
                        color: #bebebe;
                        font-size: .75em;
                        text-shadow: 0 1px 0px rgba(0,0,0,.3);  
                        background: rgba(0,0,0,.1);
                        border: 1px solid rgba(0,0,0,.12);
                        padding: 8px 10px;
                        margin: 0 14px 0 0;
                        position: relative;
                        top: 5px;
                        -moz-box-shadow:       0 1px 0px rgba(255,255,255,.15);
                        -webkit-box-shadow:    0 1px 0px rgba(255,255,255,.15);
                        box-shadow:            0 1px 0px rgba(255,255,255,.15); 
                        -moz-border-radius:    3px;
                        -webkit-border-radius: 3px;
                        border-radius:         3px;
                      }

                      #footer_info, #footer_info strong, #footer p, #footer p a {
                        color: #e1e1e1;
                        text-shadow: 0 1px 1px #000; 
                      }


                      .short_url {
                        text-shadow: 0 1px 0px rgba(0,0,0,.3);  
                        background: rgba(0,0,0,.1);
                        border: 1px solid rgba(0,0,0,.12);
                        *border: 1px solid #232323;
                        padding: 8px 10px;
                        color: #bebebe;
                        margin: 0 0 0 0;
                        -moz-box-shadow:       0 1px 0px rgba(255,255,255,.15);
                        -webkit-box-shadow:    0 1px 0px rgba(255,255,255,.15);
                        box-shadow:            0 1px 0px rgba(255,255,255,.15);
                        -moz-border-radius:    3px;
                        -webkit-border-radius: 3px;
                        border-radius:         3px;  
                        position: relative;
                        top: -30px;
                      }

                      .short_url a {
                        color: #fff;
                      }

                      /* */

                    </style>

            <meta http-equiv="x-dns-prefetch-control" content="off"/>
            </head>
            <body>
                <div id="content" class="clearfix">

                    <!-- begin #header -->      
                    <div id="header" class="clearfix">
                        <h1>
                              <a class="no_header">%s</a>
                        </h1>
                      </div>
                    <!-- end #header -->  
        ''' % (title, title)
        
        
    def create_sidebar(self,linkdict={"Home":"./index.html","A website":"http://www.google.com"},
        sidebar_title='',logo_url="http://i.imgur.com/rWYK2.png"):
        sidelinks = ''
        for title, link in linkdict.iteritems():
            sidelink = '''<a class="page_link" href="%s"><span>%s</span></a>
            ''' % (link,title)
            sidelinks += sidelink
            
        self.sidebar = '''
        <!-- begin #sidebar -->      
        <div id="sidebar">
            <div id="description">
                    <!-- logo for sidebar -->      
                 <img id="tumblr_portrait" src="%s" width="128" height="128" alt="Q"/>

                    <!-- title under sidebar logo -->      
                    %s
            </div>
                
                <!-- sidebar links -->      
                %s
        </div>
        <!-- end #sidebar -->
        
            <!-- begin .postcontainer -->
            <div class="postcontainer clearfix">        
                <!-- begin .post_inner_container -->      
                <div class="post_inner_container clearfix">
                ''' % (logo_url,sidebar_title,sidelinks)
                
    def create_footer(self,footercontent=""):
        self.footer = '''
                              <div id="footer">
                                <div id="footer_content" class="clearfix">
                                </div>                      
                                <div id="footer_info">
                                    <p><strong>Brutal Simplicity</strong> theme by <a href="http://kevin.tumblr.com">Kevin Burg</a></p>
                                  %s
                                  </div>
                            </div>
                        </div>
                        <!-- end .post_inner_container -->
                    </div> 
                    <!-- end .postcontainer -->
                </div>
            </body>
        </html>''' % (footercontent)
        
    def add_plain_html(self,content=''):
        plainhtml = '''%s
        ''' % content
        self.plainhtmls += plainhtml
        
    def add_post(self,title='',content=''):
        post = '''
        <!-- begin a post -->
        <div class="post clearfix">
            <div class="date_and_notes">
                <!-- Post Title (if desired) -->
                %s
                <!-- end post title -->
            </div>
            <div class="regular clearfix">
                <div class="post_semi_top">
                </div>
                <div class="regular_body">
                    <!-- Post Content -->
                    %s
                    <!-- end post content -->
                </div>
                 
            </div>
        </div>
        <!-- end a post -->
        
        ''' % (title,content)
        self.posts.append(post)
    
    def export_html(self):
        self.html_block = ""
        self.html_block += self.header
        self.html_block += self.sidebar
        for post in self.posts:
            self.html_block += post
        for plainhtml in self.plainhtmls:
            self.html_block += plainhtml
        self.html_block += self.footer

        try:
            filename = self.out_dir + '/index.html'
            f = open(filename,'w')
            f.write(self.html_block)
            f.close()
            self.successful_export = True
        except:
            self.successful_export = False