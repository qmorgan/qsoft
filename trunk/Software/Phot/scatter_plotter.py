import matplotlib
import os
import numpy
from pylab import *
from numpy import array as arr
import cosmocalc

#from RedshiftMachine import CollectGRBInfo
#grbdb = CollectGRBInfo.LoadDB('GRB_full') # load the database pickle file. if it doesn't exist yet, it will be created
#grbdb.Reload_DB() # takes values from the dictionary and regenerates them as arrays that are attributes of the grbdb object itself, making for easy plotting.  e.g. it will create an attribute called grbdb.Z, which is a dictionary all the redshift information: 

#storepath = os.environ.get("Q_DIR") + '/store/'

# GRB list at 3 minutes
GRB_list = ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '090709A', '070208', '071025', '080319A', '090418A']
j_3min = [8052.06, 3050.04, 324.251, 20082.0, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423, 112.846, 126.616]

GRB_list =  ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '090709A', '070208', '071025']
beta_3min = []

GRB_with_z = ['080310', '080330', '090618', '080319C', '080607', '051109A', '070208', '071025']
z_list = [2.4266, 1.51, 0.54, 1.95, 3.036, 2.346, 1.165, 4.8]

# 090530, 090709A not known yet
GRB_with_z_limits = ['061126','080310', '080330', '090618', '080319C', '090530','080607', '051109A', '090709A','070208', '071025']
z_list_limits = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8]

spectral_index = []
# For beta_solve: P1 = ([pj, ph, pk], [pj_err, ph_err, pk_err]) As per usual: pj: [a_a, a_b, tbreak, sh, flux]. Use P2 if there is second component.
GRB_071025_P = ([-1.576, 1.782, 574.66, -1, 6100], [1.242, -10.218, 1436.927, -1, 1394.88])

GRB_071025_P1 = (([-1.576, 1.782, 574.66, -1, 4544.1], [-1.576, 1.782, 574.66, -1, 6100.18], [-1.576, 1.782, 574.66, -1, 10109.1]), ([0.221, 0.347, 1, 1 ,1], [0.221, 0.347, 1, 1 ,1], [0.221, 0.347, 1, 1 ,1]))

GRB_071025_P2 = (([1.242, -10.218, 1436.927, -1, 1039.06], [1.242, -10.218, 1436.927, -1, 1394.88], [1.242, -10.218, 1436.927, -1, 2311.57]), ([0.057, 2.714, 1, 1 ,1], [0.057, 2.714, 1, 1 ,1], [0.057, 2.714, 1, 1 ,1]))

def get_flux(burst_data,t,ftopt,band):
    '''
    do this in the burst#_fitting folder.
    burst_data = file path to the burst_data.dat file containing photometry
    t = time in seconds you wish to obtain the flux
    ftopt = file path to the ftopt fitting parameter file 
    band = string of desired band (J, H, or K)
    '''
    band = band.upper() #making sure band in upper case
    import pidly
    idl_path = '/Applications/itt/idl71/bin/idl'
    idl = pidly.IDL(idl_path)
    #Performing Fit
    
    IDL_command = "lcurve, '" + str(burst_data) + "', reffilt='J', ftopt='" + str(ftopt) +"', timeunit='sec', yr=[21,12], tspec=" + str(t)+"., captionfmt='(F6.2)', /residual"
    idl(IDL_command)

    #Read the filename
    filename = burst_data.rstrip('.dat') + 'datased.dat'
    f=open(filename,'r')
    lines = f.readlines()
    for line in lines:
        if band in line:
            flux = (line.split()[1], line.split()[2])

    return flux

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
    # J:
    for index, GRB in enumerate(all_GRB_dict):
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
            if 'j' in epoch:
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

        if index >= 6:
            if len(maglist):
                matplotlib.pyplot.errorbar(timelist, maglist, yerr=mag_errlist, xerr=time_errlist, \
                        marker = '+', linestyle = 'None', label = GRB)
            if len(upper_list):
                matplotlib.pyplot.errorbar(time_upperlist, upper_list , xerr=time_upper_errlist,linestyle = 'None', marker = 'v')
        else:
            if len(maglist):
                matplotlib.pyplot.errorbar(timelist, maglist, yerr=mag_errlist, xerr=time_errlist, \
                        marker = 'o', linestyle = 'None', label = GRB)
            if len(upper_list):
                matplotlib.pyplot.errorbar(time_upperlist, upper_list , xerr=time_upper_errlist,linestyle = 'None', marker = 'v')
             
    ax = matplotlib.pyplot.gca()
    ax.set_ylim(ax.get_ylim()[::-1]) # reversing the ylimits
    matplotlib.pyplot.xlabel('Time since Trigger (s)', fontsize=19)
    matplotlib.pyplot.ylabel('Mag (J)', fontsize=21)
    matplotlib.pyplot.legend()
    ax = matplotlib.pyplot.gca()
    ax.set_xscale('log')
    xlim([80, 2*10**5])
    savepath ='./all_lightcurve_nocorrect_j2.eps'
    print 'lightcurve saved to ' + savepath
    matplotlib.pyplot.savefig(savepath)
    matplotlib.pyplot.close()
    
def Time_corrected_plot():
    Burstlist = update_z_all_GRBs(GRB_with_z_limits, z_list_limits)
    correct_time_all(Burstlist)
    plotall(Burstlist)
#    return Burstlist

def makedict(GRB_list, very_good_pickle_path='/Users/pierrechristian/qrepo/store/picklefiles/very_good_pickles/'):
    ''' make dict for all_bursts() '''
    from Phot import q_phot
    from MiscBin import qPickle
    from glob import glob

    all_GRB_dict = {}

    for index, GRB in enumerate(GRB_list):
        globstr = very_good_pickle_path + GRB + '*'
        pathstr = glob(globstr)[0]
        print pathstr
        result = qPickle.load(pathstr)
        GRB_dict = {GRB:result}
        all_GRB_dict.update(GRB_dict)
    return all_GRB_dict

def all_bursts():
    #Burstlist = update_z_all_GRBs(GRB_list)
    GRB_list =  ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '070208', '071025']
    z_list_bad = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8, 0, 0]
    Burstlist = makedict(GRB_list)
    plotall(Burstlist, time_correct=False)
    return GRB_list

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
    from numpy import array as arr
    
    matplotlib.pyplot.clf()
    #Ask Adam: difference between different bands: only in the flux multiplication?
#    P1=numpy.array(P1)
#    if P2:
#        P2=numpy.array(P2)
    if P2:
        j_mag = multibeuermann(t, P1[0][0], sec_comp=P2[0][0])
        print 'j_mag ', j_mag
        h_mag = multibeuermann(t, P1[0][1], sec_comp=P2[0][1])
        print 'h_mag ', h_mag
        k_mag = multibeuermann(t, P1[0][2], sec_comp=P2[0][2])
        print 'k_mag ', k_mag
    else:
        j_mag = multibeuermann(t, P1[0][0])
        print 'j_mag ', j_mag
        h_mag = multibeuermann(t, P1[0][1])
        print 'h_mag ', h_mag
        k_mag = multibeuermann(t, P1[0][2])
        print 'k_mag ', k_mag

    k_flux = numpy.log10(q.mag2flux(k_mag,0,6.667e-21)) # zeropoints are in erg/s/cm^2/Hz, convert to AB mag
    h_flux = numpy.log10(q.mag2flux(h_mag,0,1.024e-20))
    j_flux = numpy.log10(q.mag2flux(j_mag,0,1.594e-20)) 

    specmag = numpy.array([k_flux, h_flux, j_flux])

    # Calculating the errors in flux using the errors in the parameters
    if P2:
        print 'P1 ', arr(P1[0][0])
        print 'P1+err ', arr(P1[0][0])+arr(P1[1][0])
        j_mag_with_err = multibeuermann(t, (arr(P1[0][0])+arr(P1[1][0])), sec_comp=list(arr(P2[0][0])+arr(P2[1][0])))
        print 'j_mag_with_err ', j_mag_with_err
        h_mag_with_err = multibeuermann(t, (arr(P1[0][1])+arr(P1[1][1])), sec_comp=list(arr(P2[0][1])+arr(P2[1][1])))
        print 'h_mag_with_err ', h_mag_with_err
        k_mag_with_err = multibeuermann(t, (arr(P1[0][2])+arr(P1[1][2])), sec_comp=list(arr(P2[0][2])+arr(P2[1][2])))
        print 'k_mag_with_err ', k_mag_with_err
                                                                                                         
    else:
        j_mag_with_err = multibeuermann(t, arr(P1[0][0])+arr(P1[1][0]))
        print 'j_mag_with_err ', j_mag_with_err
        h_mag_with_err = multibeuermann(t, arr(P1[0][1])+arr(P1[1][1]))
        print 'h_mag_with_err ', h_mag_with_err
        k_mag_with_err = multibeuermann(t, arr(P1[0][2])+arr(P1[1][2]))
        print 'k_mag_with_err ', k_mag_with_err
 
    j_flux_err = numpy.log10(q.mag2flux(j_mag_with_err,0,1.594e-20)) 
    h_flux_err = numpy.log10(q.mag2flux(h_mag_with_err,0,1.024e-20))
    k_flux_err = numpy.log10(q.mag2flux(k_mag_with_err,0,6.667e-21))
        
    specmag_error = numpy.array([numpy.absolute(j_flux_err-j_flux), numpy.absolute(h_flux_err-h_flux), numpy.absolute(k_flux_err-k_flux)])# What should be the errors?

    print 'j_flux ',j_flux
    print 'j_flux_err ', j_flux_err

    p_freqs = numpy.array([1.394383526e+14,1.816923988e+14,2.398339664e+14]) # placeholder values, get the actual values!
    log_freqs = log10(p_freqs)

    print specmag
    log_specmag_err = numpy.array(specmag_error)/numpy.array(specmag) # why is this line here?
    
    results = pofit.fit(p_freqs, specmag, specmag_error, 'beta', name='spectral fit')
    return results

def beta_solve_test(t):
    result = beta_solve(t, GRB_071025_P1, GRB_071025_P2)
    return result

def plot_lum():
    clf()
    j_3min = [8052.06, 3050.04, 324.251, 20082.0, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423, 112.846, 126.616]
    j_3min2 = [8052.06, 3050.04, 324.251, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423, 112.846, 126.616]
    j_3min3 = [3050.04, 324.251, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423, 112.846, 126.616]

    j_3min = [8052.06, 3050.04, 324.251, 20082.0, 1443.05, 1070.26, 1879.54, 3210.33, 312.932, 233.877, 714.423, 188.211, 1594, 57.29, 833466.82317]

    #convert to cgs from microjansky:
    j_3min = arr(j_3min)*10**(-29)

    #convert to AB magnitude:
    j_3min = -2.5*numpy.log10(j_3min) - 48.60
    
    hist(j_3min,13)
    xlabel('$m_j$', fontsize=28)
    ylabel('Number', fontsize=28)
    yticks(arr([0, 1., 2., 3., 4.]))
    ax = matplotlib.pyplot.gca()
    ax.set_xlim(ax.get_xlim()[::-1]) # reversing the xlimits
    savefig('Lum_dist.eps')

    clf()
 #   hist(j_3min,20,cumulative=True, histtype='step')
  #  hist(j_3min2,20,cumulative=True, histtype='step')
   # hist(j_3min3,20,cumulative=True, histtype='step')
    #ylim(0,14)
    #xlim(-1000,22000)
    #xlabel('J Flux at 3 Minutes (Micro Jansky)')
  #  savefig('lum_dist.eps')
    return j_3min

def betaplot():
    beta = [-1.35, -0.8, -0.96, -0.22, -1.73, -0.84, -3.48, -0.42, -3.81, -0.3, -1.7, -0.47]
    hist(beta, 10)

def plot_lum_rest():
    '''f_{rest,V} = f_{rest_corr}*[nu_V/ ((1+z)nu_J)]^beta for  flux \propto nu^beta and beta negative values'''
    clf()
    f_rest_corr = [2252.14, 1626.48, 403.717, 11783.2, 913.329, 549.616, 286.863, 990.110, 14.7689, 174.540, 1419.79, 149309.80115] 
    beta = [-1.35, -0.8, -0.96, -0.22, -1.73, -0.84, -3.48, -0.42, -3.81, -0.3, -1.7, -0.47]
    z_list_limits = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8, 0.9382]
    arrf = arr(f_rest_corr)
    arrb = arr(beta)
    arrz = arr(z_list_limits)
    nu_V = 5.444646098003629764065335753176043557e+14 
    nu_J = 2.398339664e+14
    f_rest_V = arrf * (nu_V/ ((1+arrz)*nu_J))**(beta)
    print 'f_rest_V in microjansky:'
    print f_rest_V

    #convert to cgs from microjansky:
    f_rest_V = f_rest_V*10**(-29)
    print 'f_rest_V in cgs:'
    print f_rest_V

    #get luminosity distance from cosmocalc (lambdaCDM: omega_M = 0.27 and omega_lambda=0.73)
    dist = []
    for redshift in z_list_limits:
        dist += [cosmocalc.cosmocalc(z=redshift)['DL_cm']]
    
    arrd = arr(dist)
    print 'dist:'
    print arrd

    L_rest_V = f_rest_V*4*numpy.pi*arrd**2./(1.+arrz)
    print 'L_rest_V:'
    print L_rest_V
    #convert to ABSOLUTE AB magnitude:
    parsec = 3.085677581e18 # cm
    F_10pc = L_rest_V/(4 * numpy.pi * (10*parsec)**2)    #flux density at 10 parsecs     
    Absol_Mag = -2.5*numpy.log10(F_10pc) - 48.60    #Absolute mag in AB mag
    
    hist(Absol_Mag,6)
    xlabel('$M_v$', fontsize=27)
    ylabel('Number', fontsize=28)
    yticks(arr([0, 1., 2., 3., 4.]))
    ax = matplotlib.pyplot.gca()
    ax.set_xlim(ax.get_xlim()[::-1]) # reversing the xlimits
    savefig('Lum_dist_rest.eps')

    print 'Done'
    
    return Absol_Mag
