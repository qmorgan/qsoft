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
