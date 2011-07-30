from amara import bindery
import datetime
from matplotlib import pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import os
from matplotlib import rc
import copy
import cosmocalc
from AutoRedux.Signal import DownloadFile


overplot_high_z = True
#check out yuml.me for flow chart creator

#sent to me from JSB on 10/01/11

if not os.environ.has_key("Q_DIR"):
    print "You need to set the environment variable Q_DIR to point to the"
    print "directory where you have Q_DIR installed"
    sys.exit(1)
storepath = os.environ.get("Q_DIR") + '/store/'
loadpath = os.environ.get("Q_DIR") + '/load/'


# Import matplotlib run commands to set font
rc('font', family='serif')
default_filename=storepath+"grboxtxt.xml"
rc('text', usetex=True)

# xml_url = 'http://lyra.berkeley.edu/grbox/grboxtxt.php?form=submitted&starttime=700101&endtime=111231&sort=z&reverse=y&showindex=y&showt90=y&showepeak=y&showfluence=y&showra=y&showdec=y&showerr=y&showb=y&showebv=y&showz=y&showhostmag=y&showut=y&showxflux=y&showpeakmag=y&shownh=y&showsat=y&showclass=y&comments=y&xor=y&observatory=t&obsdate=2011-07-25&posfmt=sexc&xrtpos=best&format=xml'


def parse_grbox_xml(ignore_nonswift=True,filename=default_filename,cutoff_date=(2010,07,01),
                    exclude=[], ignore_short=True, ignore_xrf=False, ignore_long=False):
    
    a = bindery.parse(filename)

    # Parse the xml to grab the grb names and whether there have been certain 
    # types of detections (radio, optical, xray)
    names          = [x.xml_value.encode() for x in a.xml_select(u'/grbs/grb/@index')]
    has_xray  = [x.xml_value.encode() == 'y' for x in a.xml_select(u'/grbs/grb/greiner/@x')]
    has_opt   = [x.xml_value.encode() == 'y' for x in a.xml_select(u'/grbs/grb/greiner/@o')]
    has_radio = [x.xml_value.encode() == 'y' for x in a.xml_select(u'/grbs/grb/greiner/@r')]
    

    # Define lists to use later
    has_host_z = []
    ignores = ['hostphotz','photz']
    
    grbox_dict = {}
    
    # loop through all the GRB names in the XML file and grab all the redshifts
    # associated with each name (there may be more than one!) 
    for grbname in names:
        ## get all the redshifts
        if grbname in exclude:
            exclude.remove(grbname)
            continue
        print grbname
        # Obtain a list of tentative redshifts for every grbname
        instrument = [str(x) for x in a.xml_select(u'/grbs/grb[@index="%s"]/burst/instrument' % grbname)][0]
        grb_class = [str(x) for x in a.xml_select(u'/grbs/grb[@index="%s"]/burst/class' % grbname)]
        if not grb_class:
            grb_class = ['GRB']  #unlabeled classes are assumed to be just normal bursts
        grb_class = grb_class[0]
        tentative_z =[str(x) for x in a.xml_select(u'/grbs/grb[@index="%s"]/redshift/z' % grbname)]
        # if there is an "?" after the redshift number (e.g. <z>1.24?</z>), mark it
        uncertain_z = [str(x).find("?") != -1 for x in tentative_z]
        # Obtain a list of the redshift types
        z_type = [str(x) for x in a.xml_select(u'/grbs/grb[@index="%s"]/redshift/ztype' % grbname)]
        print grbname, [(tentative_z[i],uncertain_z[i],z_type[i]) for i in range(len(tentative_z))]
        use = False
        for i in range(len(tentative_z) - 1,-1,-1):
            use_this_z = True
            for ii in ignores:
                if z_type[i].find(ii) != -1:
                    # this type is to be ignored
                    use_this_z=False
            has_host = False
            if z_type[i].find("hostz") != -1:
                has_host = True
            
            if use_this_z and not uncertain_z[i]:
                use = True
                zz = tentative_z[i]
        if use:
           
            yr = int("19" + grbname[0:2]) if int(grbname[0]) > 8 else int("20" + grbname[0:2])
            mn = int(grbname[2:4])
            dy = int(grbname[4:6])
           
            if cutoff_date:
                cutoff_yr = int(cutoff_date[0])
                cutoff_mo = int(cutoff_date[1])
                cutoff_day = int(cutoff_date[2])
                if datetime.date(yr,mn,dy) > datetime.date(cutoff_yr,cutoff_mo,cutoff_day):
                    continue
            # if ignore_nonswift and datetime.date(yr,mn,dy) < datetime.date(2004,12,10):
               
            if ignore_nonswift and instrument != 'Swift':
                continue
            if ignore_short and grb_class == 'SHB':
                continue
            if ignore_xrf and grb_class == 'XRF':
                continue
            if ignore_long and grb_class == 'GRB':
                continue
               
            has_host_z.append(has_host)
           
           
            subdict = {grbname:{'date':datetime.date(yr,mn,dy),'class':grb_class,'instrument':instrument,'grbox_z':float(zz),'has_host_z':has_host}}
           
            grbox_dict.update(subdict)

    
    from pprint import pprint 
    
    pprint(grbox_dict)
    
    return grbox_dict


def Make_Z_Plot(grbox_dict):
    z_list = []
    zname_list = []
    date_list = []


    # 'un'zip all_grbs (is there a better way to do this? just assign rather than a for loop?)
    for key,value in grbox_dict.iteritems():
        date_list.append(value['date'])
        z_list.append(value['grbox_z'])
        zname_list.append(key)
    print '***All***'
    print ' '
    print z_list
    print len(z_list)


    ### Print out the most distant GRB as a function of time
    zmax = 0.0
    rr = []
    for key,value in grbox_dict.iteritems():
        if value['grbox_z'] > zmax:
            rr.append(key)
            zmax = value['grbox_z']
            print value['date'].year + value['date'].timetuple().tm_yday/365.0, zmax, "#  ", key
    

    ax = plt.subplot(111)
    n, bins, patches = plt.hist(plt.log10(z_list),bins=29,facecolor='grey',edgecolor='grey',alpha=0.95)

    # Define pre-swift burst index as bursts before 041210
    high_z_i = plt.where(plt.array(date_list) < datetime.date(2004,12,10))

    high_z_list = [z_list[i] for i in list(high_z_i[0])]
    #print high_z_list
    n, bins1, patches = plt.hist(plt.log10(high_z_list),bins=bins,facecolor='red',edgecolor='red',alpha=0.6)

    if overplot_high_z:
        high_z_list = [z for z in z_list if z > 4.0]
        n, bins1, patches = plt.hist(plt.log10(high_z_list),bins=bins,facecolor='red',edgecolor='red')


    ay = ax.twinx()

    argg = list(plt.ones(len(z_list)).cumsum().repeat(2))
    zz = copy.copy(z_list)
    zz.sort()
    tmp = list(plt.log10(zz).repeat(2))

    tmp.append(1)
    yy = [0]
    yy.extend(argg)
    
    
    ay.plot(tmp,yy,aa=True,linewidth = 4,color='black',alpha=0.95)

    argg = list(plt.ones(len(high_z_list)).cumsum().repeat(2))
    zz = copy.copy(high_z_list)
    zz.sort()
    tmp = list(plt.log10(zz).repeat(2))

    tmp.append(1)
    yy = [0]
    yy.extend(argg)


    ay.plot(tmp,yy,aa=True,linewidth = 2,color='#222222',alpha=0.75)
    ay.set_ylim((0,len(z_list)*1.05))
    ay.set_ylabel("Cumulative Number",fontsize=20)
    # formatter for bottom x axis 
    def ff(x,pos=None):
        if x < -1:
            return "%.2f" % (10**x)
        elif x < 0:
            return "%.1f" % (10**x)
        elif 10**x == 8.5:
            return "%.1f" % (10**x)
        else:
            return "%i" % (10**x)

    formatter = FuncFormatter(ff)
    ax.set_xticks([-2,-1,plt.log10(0.3),0,plt.log10(2),plt.log10(3),plt.log10(4),plt.log10(6),plt.log10(8.5)])
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlabel("Redshift ($z$)",fontsize=20)
    ax.set_ylabel("Number",fontsize=20)

    ax.set_xlim( (plt.log10(0.005),plt.log10(10)))

    ax2 = ax.twiny()
    xlim= ax.get_xlim()
    #ax2.set_xscale("log")
    ax2.set_xlim( (xlim[0], xlim[1]) )

    # Define function for plotting the top X axis; time since big bang in Gyr
    def rr(x,pos=None): 
        g = cosmocalc.cosmocalc(10.0**x, H0=71.0)
        if g['zage_Gyr'] < 1:
            return "%.2f" % g['zage_Gyr'] # Return 2 dec place if age < 1; e.g 0.62
        else:
            return "%.1f" % g['zage_Gyr'] # Return 1 dec place if age > 1; e.g. 1.5

    ax2.set_xticks([-1.91,-1.3,-0.752,-0.283,0.102,0.349,0.62,plt.log10(8.3)])

    formatter1 = FuncFormatter(rr)
    ax2.xaxis.set_major_formatter(formatter1)
    ax2.set_xlabel("Time since Big Bang (Gyr)",fontsize=20)

    #plt.bar(l,a['yy'],width=w,log=False)
    #ax.set_xscale("log",nonposx='clip')

    ## Now plot inset plot of GRBs greater than z=4.0

    axins = inset_axes(ax2,
                        width="30%", # width = 30% of parent_bbox
                        height="30%") # height : 1 inch)

    locator=axins.get_axes_locator()
    locator.set_bbox_to_anchor((-0.8,-0.45,1.35,1.35), ax.transAxes)
    locator.borderpad = 0.0

    high_z_list = [z for z in z_list if z > 4.0]
    if high_z_list: # if there are any high-z's, plot them
        n, bins, patches = plt.hist(plt.array(high_z_list),facecolor='red', edgecolor='red')
    axins.set_xlim(4.0,8.5)
    axins.set_xlabel("z")
    axins.set_ylabel("N")

    high_z_i = plt.where(plt.array(date_list) < datetime.date(2004,12,10))
    high_z_list = [z_list[i] for i in list(high_z_i[0]) if z_list[i] > 4.0]

    n, bins, patches = plt.hist(plt.array(high_z_list),bins=bins,facecolor='red',edgecolor='red')

    high_z_list = [z for z in z_list if z > 4.0]

    if high_z_list: # if there are any high-z's, plot them
        n, bins, patches = plt.hist(plt.array(high_z_list),facecolor='red',edgecolor='red')

    axins.set_xlim(4.0,9.0)
    #mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    #axins2.set_xlabel("Time since Big Bang [Gyr]",fontsize=20)
    ylabels = ax.get_yticklabels()
    plt.setp(ylabels, size=14, name='times', weight='light', color='k')

    xlabels = ax.get_xticklabels()
    plt.setp(xlabels, size=14, name='times', weight='light', color='k')

    xlabels = ax2.get_xticklabels()
    plt.setp(xlabels, size=14, name='times', weight='light', color='k')

    xlabels = ay.get_yticklabels()
    plt.setp(xlabels, size=14, name='times', weight='light', color='k')


    plt.draw()
    return z_list
#for i in range(len(z_list)):
#    print z_list[i], uncertain_z[i],ztype[i],names[i]

def MakePlotForRATEGRBz():
    grbdict = parse_grbox_xml(ignore_nonswift=True,filename=default_filename,cutoff_date=(2010,07,01),
                        exclude=['071020','050724','070612A','050126','050401','071010B','071003','091024','071010A','050603','090313',
                        '090726', '090313','070283','071028B','071010B','071010A','071003','061210','060814','060729','060614'], # this exclude file needs to be properly fleshed out due to observing constraints.
                        ignore_short=True, ignore_xrf=False, ignore_long=False)
    Make_Z_Plot(grbdict)
    return grbdict