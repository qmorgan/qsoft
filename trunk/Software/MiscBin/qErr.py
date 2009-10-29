#!/usr/bin/env python
import sys, traceback
from AutoRedux import send_gmail
import os, time 

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'

def qErr(mail=True,mail_to='qmorgan@gmail.com'):
    '''Output the traceback of an error to a log file, and send the tail
    of the log file to someone by email, if desired.  Most useful in the 
    try/except, i.e. when you want information that something failed, but 
    don't want the program to crash.
    
    Usage:
    
    >>> a=[1,2,3]
    >>> try: 
    ...    print a[6]
    ... except:
    ...    qErr()
    ...
    Exception Triggered, traceback info forwarded to log file.
      
    the following will be appended to the log and emailed:
    ------------
    Wed Oct 28 16:47:46 2009
    Traceback (most recent call last):
      File "<stdin>", line 2, in <module>
    IndexError: list index out of range
    '''
    print "Exception Triggered, traceback info forwarded to log file."
    errfilepath = storepath + 'errlog.txt'
    
    errfile = open(errfilepath,"a")
    errfile.write('------------\n')
    errfile.write(time.ctime())
    errfile.write('\n')
    traceback.print_exc(file=errfile)
    errfile.write('\n')
    errfile.close()
    if mail:
        errfilemailpath = storepath + 'errlogmail.txt'
        err_mail_sub = 'An error has occured!'
        err_mail_text = 'Please see attached for error log.'
        
        operation = "tail -n 20 %s > %s" % (errfilepath, errfilemailpath)
        os.system(operation) 
        send_gmail.domail(mail_to,err_mail_sub,err_mail_text,[errfilemailpath])
