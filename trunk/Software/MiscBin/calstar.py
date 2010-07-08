from Phot import q_super_photometry
import matplotlib
from scipy import array
from MiscBin import t_mid
import os
import pylab

storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'

def magplot(reg, filelist):
    
    '''temporary comment: Plot magnitudes of calibration stars, needs q_super_photometry and t_mid.'''

    matplotlib.pyplot.clf()
    regpath = storepath + reg
    temppath = storepath + 'temp.reg'
    regfile = open(regpath, 'r')
    reglist = regfile.readlines()
    callist = []
    for line in reglist:
        if 'circle' in line:
            callist += [line]
        else:
            pass
    
    colornumber = len(callist)
    
    for index, star_reg in enumerate(callist):
        if os.path.exists(temppath):
            os.remove(temppath)
        datalist = []
        dataerrlist = []
        timelist = []
        timeerrlist = []
        colorstr = str(float((1/colornumber))*float(index + 1))
        colortuple = (colorstr, 0.5, 0)
        starname = 'star'+str(index)
        tempreg = open(temppath, 'w')
        tempreg.write('# Region file format: DS9 version 4.1\n')
        secondstr='global color=green dashlist=8 3 width=2 font="helvetica '+ \
                 '16 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '+ \
                 'delete=1 include=1 source=1\n'
        tempreg.write(secondstr)
        tempreg.write('fk5\n')
        tmp_str = star_reg[:-37]
        tempreg.write(tmp_str)
        tempreg.close()

        for image in filelist:
            print '**************************************'
            print 'Photometry of star' + str(index) 
            print 'doing image ' + image
            data = q_super_photometry.dophot(image, temppath)
            if 'targ_mag' in data:
                datalist += [data['targ_mag'][0]] 
                dataerrlist += [data['targ_mag'][1]]
                time = float(t_mid.t_mid(image))
                terr = float(t_mid.t_mid(image, delta = True))/2.
                timelist += [time]
                timeerrlist += [terr]
            else:
                pass

        datarr = array(datalist)
        daterrarr = array(dataerrlist)
        timarr = array(timelist)
        timerrarr = array(timeerrlist)

        matplotlib.pyplot.errorbar(timarr, datarr, yerr = daterrarr, label = starname, fmt='k.', color = colortuple) 
         
        matplotlib.pyplot.title('Calibration Stars Magnitude vs. t_mid')
        matplotlib.pyplot.xlabel('Time After Burst (s)')
        matplotlib.pyplot.ylabel('Magnitude')
        ax = matplotlib.pyplot.gca()
        ax.set_ylim(ax.get_ylim()[::-1])
        matplotlib.pyplot.legend()

    unique_name = (filelist[0].split('_'))[2]
    filepath = storepath + unique_name + '_calibration_stars.png'
    matplotlib.pyplot.savefig(filepath)

        
