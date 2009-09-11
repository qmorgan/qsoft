#!/usr/bin/env python
# encoding: utf-8
"""
MakeDocs.py
Author: Adam Morgan

This Program uses pydoc to create documentation for 
"""
import sys
import os
import glob
import socket

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
docspath = os.environ.get("Q_DIR") + '/docs/'
softpath = os.environ.get("Q_DIR") + '/trunk/Software/'
hostname = socket.gethostname()

def MakeDocs():
    modlist = glob.glob(softpath + '*/')
    print modlist
    indexbody = ''
    cmd1 = 'cd %s' % (softpath)
    os.system(cmd1)


    for mod in modlist:
        modname=os.path.dirname(mod).split('/')[-1]
        print "Now preparing documentation for %s" % (modname)
        add_to_index = '''<a href="%s.html">%s</a><br>
        ''' % (modname,modname)
        indexbody += add_to_index
        
        cmd2 = 'pydoc -w %s' % (modname)
        os.system(cmd2)
    
        proglist = glob.glob(mod+'/*.py')
        for prog in proglist:
            progname=os.path.basename(prog)
            if progname[0:2] != '__':
                outname = modname+'.'+progname.rstrip('py').rstrip('.')
                
                add_to_index= '''&nbsp;&nbsp;&nbsp;-&nbsp;<a href="%s.html">%s</a><br>
                ''' % (outname,progname.rstrip('py').rstrip('.'))
                indexbody += add_to_index
                
                cmd2 = 'pydoc -w %s' % (outname)
                os.system(cmd2)
        

    cmd3 = 'mv *.html %s' % (docspath)
    os.system(cmd3)
    
    MakeIndex(indexbody)
    
    if hostname == 'goofy.Berkeley.EDU':
        cmd1 = 'cp %s/*.html /o/amorgan/public_html/pydoc/.' % (docspath)
        try:
            print 'Moving files to Public Html folder'
            os.system(cmd1)
        except:
            print 'Cannot copy to public html folder'
    
def MakeIndex(indexbody):
    
    
    indexhead=''' 
    <!doctype html PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
    <html><head><title>Python: q_soft packages</title>
    </head><body bgcolor="#f0f0f8">

    <table width="100%" cellspacing=0 cellpadding=2 border=0 summary="heading">
    <tr bgcolor="#7799ee">
    <td valign=bottom>&nbsp;<br>
    <font color="#ffffff" face="helvetica, arial">&nbsp;<br><big><big><strong>Q_SOFT</strong></big></big></font></td><td align=right valign=bottom><font color="#ffffff" face="helvetica, arial"><a href=".">index</a><br><a href="file:/goofy1/amorgan/q_soft/trunk/Software/">/goofy1/amorgan/q_soft/trunk/Software/</a></font></td></tr></table>
        <p></p>
    <p>
    <table width="100%" cellspacing=0 cellpadding=2 border=0 summary="section">
    <tr bgcolor="#aa55cc">
    <td colspan=3 valign=bottom>&nbsp;<br>
    <font color="#ffffff" face="helvetica, arial"><big><strong>Package Contents</strong></big></font></td></tr>

    <tr><td bgcolor="#aa55cc"><tt>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</tt></td><td>&nbsp;</td>
    	<td width="100%">
    		<table width="100%" summary="list"><tr>
    		<td width="25%" valign=top>
    
    '''
    indextail='''
    </td><td width="25%" valign=top></td><td width="25%" valign=top></td></tr></table></td></tr></table>
    </body></html>
    '''
    index_out = indexhead + indexbody + indextail
    f=open(docspath+'index.html','w')
    f.write(index_out)
    f.close()
    

