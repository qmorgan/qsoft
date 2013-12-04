#!/usr/bin/env python
# encodiing: utf-8
'''
upperlimitloop.py
Loops around ptel_astrophot.py to grab upper limits and zeropoints

'''

import os
import ptel_astrophot
import numpy

def ptelastroloop(f_name,iterations,position=(247.37058,40.835648)):
    index = 0
    jullist = []
    kullist = []
    hullist = []
    while index < iterations:
        a = ptel_astrophot.PTEL_data_block(f_name,ap=3,req_filts=['j','h','k'])
        a.get_phot_at_pos(pos=position,req_lim=3,force_photom=False)
        outstr = 'j %f %f h %f %f k %f %f' % \
            (a.out_results['j']['photometry']['zp'], a.out_results['j']['photometry']['limits'][3.0], \
            a.out_results['h']['photometry']['zp'], a.out_results['h']['photometry']['limits'][3.0], \
            a.out_results['k']['photometry']['zp'], a.out_results['k']['photometry']['limits'][3.0])
                # Output to file 
        jullist.append(a.out_results['j']['photometry']['limits'][3.0])
        hullist.append(a.out_results['h']['photometry']['limits'][3.0])
        kullist.append(a.out_results['k']['photometry']['limits'][3.0])
        system_command = 'echo %s >> %s.ulims' % (outstr,f_name)
        os.system(system_command)
        # Remove Pickle file before repeating
        os.system('rm -f *.pkl')
        index += 1 
    javg = numpy.average(jullist)
    havg = numpy.average(hullist)
    kavg = numpy.average(kullist)
    jstd = numpy.std(jullist)
    hstd = numpy.std(hullist)
    kstd = numpy.std(kullist)
    outstring = 'j %f %f h %f %f k %f %f' % (javg, jstd, havg, hstd, kavg, kstd)
    system_command = 'echo %s >> %s.ulims' % (outstr,f_name)
    print 'J Upper Limit: %s +/- %s' % (javg,jstd)
    print 'H Upper Limit: %s +/- %s' % (havg,hstd)
    print 'K Upper Limit: %s +/- %s' % (kavg,kstd)
    
    
    
    
