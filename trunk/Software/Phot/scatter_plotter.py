import matplotlib
import os
import numpy

#from RedshiftMachine import CollectGRBInfo
#grbdb = CollectGRBInfo.LoadDB('GRB_full') # load the database pickle file. if it doesn't exist yet, it will be created
#grbdb.Reload_DB() # takes values from the dictionary and regenerates them as arrays that are attributes of the grbdb object itself, making for easy plotting.  e.g. it will create an attribute called grbdb.Z, which is a dictionary all the redshift information: 

#storepath = os.environ.get("Q_DIR") + '/store/'

# GRB list at 3 minutes
GRB_list = ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '090709A', '070208', '071025']
j_3min = [8052.06, 3050.04, 324.251, 20082.0, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423]

GRB_with_z = ['080310', '080330', '090618', '080319C', '080607', '051109A', '070208', '071025']
z_list = [2.4266, 1.51, 0.54, 1.95, 3.036, 2.346, 1.165, 4.8]

GRB_with_z_limits = ['061126','080310', '080330', '090618', '080319C', '090530','080607', '051109A', '090709A','070208', '071025']
z_list_limits = [1.16, 2.4266, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8]

spectral_index = []
# For beta_solve: P1 = ([pj, ph, pk], [pj_err, ph_err, pk_err]) As per usual: pj: [a_a, a_b, tbreak, sh, flux]. Use P2 if there is second component.
GRB_071025_P = ([-1.576, 1.782, 574.66, -1, 6100], [1.242, -10.218, 1436.927, -1, 1394.88])

GRB_071025_P1 = (([-1.576, 1.782, 574.66, -1, 4544.1], [-1.576, 1.782, 574.66, -1, 6100.18], [-1.576, 1.782, 574.66, -1, 10109.1]), ([0.221, 0.347, 1, 1 ,1], [0.221, 0.347, 1, 1 ,1], [0.221, 0.347, 1, 1 ,1]))

GRB_071025_P2 = (([1.242, -10.218, 1436.927, -1, 1039.06], [1.242, -10.218, 1436.927, -1, 1394.88], [1.242, -10.218, 1436.927, -1, 2311.57]), ([0.057, 2.714, 1, 1 ,1], [0.057, 2.714, 1, 1 ,1], [0.057, 2,714, 1, 1 ,1]))

def update_z_all_GRBs(GRB_list, z_list, very_good_pickle_path='/Users/pierrechristian/qrepo/store/picklefiles/very_good_pickles/'):
    '''Updates the z values of all GRBs with known z'''
    from Phot import q_phot
    from MiscBin import qPickle
    from glob import glob

    all_GRB_dict = {}

    for index, GRB in enumerate(GRB_list):
        globstr = very_good_pickle_path + GRB + '*'
        pathstr = glob(globstr)[0]
        print pathstr
        result = qPickle.load(pathstr)
        q_phot.update_all_z(result, z_list[index], 0)
        GRB_dict = {GRB:result}
        all_GRB_dict.update(GRB_dict)
    return all_GRB_dict

def correct_time_all(all_GRB_dict):
    '''Perform time corrections to all GRBs with known z (done to the output of update_z_all_GRBs)'''
    for GRB in all_GRB_dict:
        for epoch in all_GRB_dict[GRB]:
            time = all_GRB_dict[GRB][epoch]['t_mid'][0]
            redshift = all_GRB_dict[GRB][epoch]['z'][0]
            corr_time = time/(1.+redshift)
            # Ask Adam: How do you deal with time uncertainty/exposure time/etc
            rest_time_dict = {'t_rest':(corr_time,all_GRB_dict[GRB][epoch]['t_mid'][1])}
            all_GRB_dict[GRB][epoch].update(rest_time_dict)

def plotall(all_GRB_dict, time_correct=True):
    '''Plots all the lightcurves in the same plot (1 for each band)'''
    import matplotlib
    import datetime
    matplotlib.pyplot.clf()
    # H:
    for GRB in all_GRB_dict:
        print '--------------'
        print '--------------'
        print 'Doing GRB %s' % (GRB)
        maglist = []
        mag_errlist = []
        upper_list = []
        timelist = []
        time_errlist = []
        time_upperlist = []
        time_upper_errlist = []
        for epoch in all_GRB_dict[GRB]:
            if 'h' in epoch:
                print '--------------'
                print 'Doing epoch %s' % (epoch)
                if 'targ_mag' in all_GRB_dict[GRB][epoch]: 
                    maglist += [all_GRB_dict[GRB][epoch]['targ_mag'][0]]
                    mag_errlist += [all_GRB_dict[GRB][epoch]['targ_mag'][1]]
                    if time_correct:
                        timelist += [all_GRB_dict[GRB][epoch]['t_rest'][0]]
                        time_errlist += [all_GRB_dict[GRB][epoch]['t_rest'][1]]
                    else:
                        timelist+= [all_GRB_dict[GRB][epoch]['t_mid'][0]]
                        time_errlist += [all_GRB_dict[GRB][epoch]['t_mid'][1]]                  
                elif 'upper_green' in all_GRB_dict[GRB][epoch]:
                    upper_list += [all_GRB_dict[GRB][epoch]['upper_green']]
                    if time_correct:
                        time_upperlist += [all_GRB_dict[GRB][epoch]['t_rest'][0]]
                        time_upper_errlist += [all_GRB_dict[GRB][epoch]['t_rest'][1]]
                    else:
                        time_upperlist += [all_GRB_dict[GRB][epoch]['t_mid'][0]]
                        time_upper_errlist += [all_GRB_dict[GRB][epoch]['t_mid'][1]]
                else:
                    print 'NO MAG OR ULIM FOUND, SKIPPING %s' % (all_GRB_dict[GRB][epoch])
                
        maglist = numpy.array(maglist)
        mag_errlist =  numpy.array(mag_errlist)
        upper_list = numpy.array(upper_list)
        timelist = numpy.array(timelist)
        time_errlist = numpy.array(time_errlist)
        time_upperlist = numpy.array(time_upperlist)
        time_upper_errlist = numpy.array(time_upper_errlist)

        if len(maglist):
            matplotlib.pyplot.errorbar(timelist, maglist, yerr=mag_errlist, xerr=time_errlist, \
                        marker = 'o', linestyle = 'None', label = GRB)
        if len(upper_list):
            matplotlib.pyplot.errorbar(time_upperlist, upper_list , xerr=time_upper_errlist,linestyle = 'None', marker = 'v')
             
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    matplotlib.pyplot.legend()
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    savepath ='./all_lightcurve_nocorrect.png'
    print 'lightcurve saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)
    matplotlib.pyplot.close()
    
def Time_corrected_plot():
    Burstlist = update_z_all_GRBs(GRB_with_z_limits, z_list_limits)
    correct_time_all(Burstlist)
    plotall(Burstlist)
#    return Burstlist

def redshift_corrections(GRBlist, z_list, spectral_index, photometry_result):
    import k_correct
    corrected_magnitude = []
    for index, burst in GRBlist:
        spectral_correction = k_correct.correct(z_list[index], spectral_index[index])
#        corrected_magnitude += 

 
def getproperty(GRBs, property):
    if property == 'fluence':
        dat = grbdb.fluence
    elif property == 't90':
        dat = grbdb.t90
    elif property == 'z':
        dat = grbdb.z

    burstlist = []
    errlist = []
    for my_burst in GRBs:
        for index, database_burst in enumerate(dat['names']):
            if my_burst == database_burst:
                burstlist += [dat['array'][index]]

    return burstlist

def scatterplot(x_prop_name, x_property, y_prop_name, y_property, x_err=None, y_err=None, name='scatter_plot'):
    '''Normal Scatter Plot'''

    if not x_err:
        x_err = list(numpy.zeros(len(x_property)))
    if not y_err:
        y_err = list(numpy.zeros(len(y_property)))
    
    matplotlib.pyplot.errorbar(numpy.log10(x_property), numpy.log10(y_property), yerr=y_err, xerr=x_err, \
                        marker = 'o', linestyle ='None', mfc = 'blue', mec = 'green', \
                        ecolor = 'blue')

    matplotlib.pyplot.xlabel(x_prop_name)
    matplotlib.pyplot.ylabel(y_prop_name)
    matplotlib.pyplot.legend()

    savepath = storepath + name + '.png'
    print 'lightcurve saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)    
    matplotlib.pyplot.close()
 
def testtry():
    a = getproperty(GRB_list, 'fluence')
    b = getproperty(GRB_list, 't90')

    scatterplot(j_3min,a)

def multibeuermann_test():
    import matplotlib
    from MiscBin import qPickle
    matplotlib.pyplot.clf()
#    s = s + f_c,f * ( (t/tbk_c,f) ^ (-sh_c * alpha_c,b) +  (t/tbk_c,f) ^ (-sh_c * alpha_c,a) )^(-1./sh_c)
    t = []
    photdict = qPickle.load('/Users/pierrechristian/qrepo/store/picklefiles/very_good_pickles/071025_goodpickle.data')
    for epoch in photdict:
        if 'h' in epoch:
            t += [float(photdict[epoch]['t_mid'][0])]
    t = numpy.array(t)
    print t
    s = -2.5*numpy.log10(6100 * ((t/574.66) ** (-1 * -1.782) + (t/574.66) ** (-1 * 1.576)) ** (-1./1)  +  1394.88  * ((t/1436.927) ** (- 1 * 10.218) + (t/1436.927) ** (-1 *-1.242)) ** (-1./1))
    matplotlib.pyplot.errorbar(t, s, linestyle = 'None', marker = 'o')
    print s
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    savepath ='./071025_multibeuermann_test.png'
    print 'lightcurve saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)
    matplotlib.pyplot.close()

def multibeuermann(t, p, sec_comp=False):    
    ''' p: [a_a, a_b, tbreak, sh, flux]'''
    #    s = s + f_c,f * ( (t/tbk_c,f) ^ (-sh_c * alpha_c,b) +  (t/tbk_c,f) ^ (-sh_c * alpha_c,a) )^(-1./sh_c)
    if sec_comp:
        s = -2.5*numpy.log10(p[4] * ((t/p[2]) ** (p[3] * (-1.*p[1])) + (t/p[2]) ** (p[3] * (-1.*p[0]))) ** (p[3]/1) + sec_comp[4]  * ((t/sec_comp[2]) ** (sec_comp[3] * (-1.*sec_comp[1])) + (t/sec_comp[2]) ** (sec_comp[3] *(-1.*sec_comp[0]))) ** (sec_comp[3]/1))
    else:
         s = -2.5*numpy.log10(p[4] * ((t/p[2]) ** (p[3] * p[1]) + (t/p[2]) ** (p[3] * p[0])) ** (p[3]/1))
    return s

def multibeuermann_test2():
    import matplotlib
    from MiscBin import qPickle
    matplotlib.pyplot.clf()
#    s = s + f_c,f * ( (t/tbk_c,f) ^ (-sh_c * alpha_c,b) +  (t/tbk_c,f) ^ (-sh_c * alpha_c,a) )^(-1./sh_c)
    t = []
    photdict = qPickle.load('/Users/pierrechristian/qrepo/store/picklefiles/very_good_pickles/071025_goodpickle.data')
    for epoch in photdict:
        if 'h' in epoch:
            t += [float(photdict[epoch]['t_mid'][0])]
    t = numpy.array(t)
    print t
    P1 = [-1.576, 1.782, 574.66, -1, 6100]
    P2 = [1.242, -10.218, 1436.927, -1, 1394.88]
    s = multibeuermann(t, P1, sec_comp=P2)
    matplotlib.pyplot.errorbar(t, s, linestyle = 'None', marker = 'o')
    print s
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    matplotlib.pyplot.xlabel('Time since Burst (s)')
    matplotlib.pyplot.ylabel('Mag')
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    savepath ='./071025_multibeuermann_test2.png'
    print 'lightcurve saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)
    matplotlib.pyplot.close()

def beta_solve(t, P1, P2=False):
    '''P1 = ([pj, ph, pk], [pj_err, ph_err, pk_err]) As per usual: pj: [a_a, a_b, tbreak, sh, flux]. Use P2 if there is second component.'''
    import numpy
    from scipy import optimize
    from scipy import log10
    from MiscBin import q
    from Phot import pofit
    
    matplotlib.pyplot.clf()
    #Ask Adam: difference between different bands: only in the flux multiplication?
    if P2:
        j_mag = multibeuermann(t, P1[0][0], sec_comp=P2[0][0])
        h_mag = multibeuermann(t, P1[0][1], sec_comp=P2[0][1])
        k_mag = multibeuermann(t, P1[0][2], sec_comp=P2[0][2])
    else:
        j_mag = multibeuermann(t, P1[0][0])
        h_mag = multibeuermann(t, P1[0][1])
        k_mag = multibeuermann(t, P1[0][2])

    k_flux = numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) # zeropoints are in erg/s/cm^2/Hz, convert to AB mag
    h_flux = numpy.log10(q.mag2flux(h_mag,0,1.024e-20))
    j_flux = numpy.log10(q.mag2flux(j_mag,0,1.594e-20)) 

    specmag = numpy.array([k_flux, h_flux, j_flux])

    # Calculating the errors in flux using the errors in the parameters
    if P2:
        j_mag_with_err = multibeuermann(t, (P1[0][0]+P1[1][0]), sec_comp=(P2[0][0]+P2[1][0]))
        h_mag_with_err = multibeuermann(t, (P1[0][1]+P1[1][1]), sec_comp=(P2[0][1]+P2[1][1]))
        k_mag_with_err = multibeuermann(t, (P1[0][2]+P1[1][2]), sec_comp=(P2[0][2]+P2[1][2]))
                                                                                                         
    else:
        j_mag_with_err = multibeuermann(t, P1[0][0]+P1[1][0])
        h_mag_with_err = multibeuermann(t, P1[0][1]+P1[1][1])
        k_mag_with_err = multibeuermann(t, P1[0][2]+P1[1][2])

    j_flux_err = numpy.log10(q.mag2flux(j_mag_with_err,0,1.594e-20)) 
    h_flux_err = numpy.log10(q.mag2flux(h_mag_with_err,0,1.024e-20))
    k_flux_err = numpy.log10(q.mag2flux(k_mag_with_err,0,6.667e-21))
        
    specmag_error = numpy.array([numpy.absolute(j_flux_err-j_flux), numpy.absolute(h_flux_err-h_flux), numpy.absolute(k_flux_err-k_flux)])# What should be the errors?

    p_freqs = numpy.array([1.394383526e+14,1.816923988e+14,2.398339664e+14]) # placeholder values, get the actual values!
    log_freqs = log10(p_freqs)

    print specmag
    log_specmag_err = numpy.array(specmag_error)/numpy.array(specmag) # why is this line here?
    
    results = pofit.fit(p_freqs, specmag, specmag_error, 'beta', name='spectral fit')
    return results

def beta_solve_test(t):
    result = beta_solve(t, GRB_071025_P1, GRB_071025_P2)
    return result
