'''Wrapper for the PTEL Early Afterglow Catalog final products'''

import matplotlib
import os
import numpy
from pylab import *
import numpy as np
from numpy import array as arr
import cosmocalc
import pandas as pd
from Phot import q_phot
from Phot import PhotParse
from MiscBin import qPickle
import glob
from GRB120119A import DustModel
from matplotlib import rc

rc('font', size=20, family='Times New Roman')

meta_dict = {'061126':{
                        'z':1.1588, #host
                        'zlim':1.1588,
                        'xmin':100,
                        'xmax':2e3,
                        'ymin':100,
                        'ymax':3e4
                        },
            # '080310': {
            #             'z':2.4266,
            #             'zlim':None,
            #             'xmin':1e3,
            #             'xmax':3e3,
            #             'ymin':6e2,
            #             'ymax':2e3
            #             },
            # '080330':{
            #            'z': 1.51,
            #            'zlim':None,
            #            'xmin':1000,
            #            'xmax':2e4,
            #            'ymin':40,
            #            'ymax':2e3
            #            },
            '090618':{
                       'z':0.54,
                       'zlim':None,
                       'xmin':100,
                       'xmax':3e3,
                       'ymin':1000,
                       'ymax':3e4
                       },
            '080319C':{
                        'z':1.95,
                        'zlim':None,
                        'xmin':100,
                        'xmax':2e3,
                        'ymin':100,
                        'ymax':5e3
                        },
            '090530':{
                       'z':1.266,
                       'zlim':1.6,
                       'xmin':100,
                       'xmax':5e3,
                       'ymin':50,
                       'ymax':2e3
                       },
            '080607':{
                       'z':3.036,
                       'zlim':None,
                       'xmin':100,
                       'xmax':5e3,
                       'ymin':30,
                       'ymax':3e4
                       },
            '051109A':{
                        'z':2.346,
                        'zlim':None,
                        'xmin':100,
                        'xmax':5e3,
                        'ymin':100,
                        'ymax':1e4
                        },
            '070208':{
                       'z':1.165,
                       'zlim':None,
                       'xmin':100,
                       'xmax':2e4,
                       'ymin':20,
                       'ymax':4e2
                       },
            '071025':{
                       'z':5.0, #approximate
                       'zlim':None,
                       'xmin':100,
                       'xmax':2e4,
                       'ymin':40,
                       'ymax':1e4
                       },
            '090709A':{
                        'z':1.8, #approximate
                        'zlim':None,
                        'xmin':100,
                        'xmax':10000,
                        'ymin':20,
                        'ymax':4e3
                        },
            '080319B':{
                        'z':0.9382, 
                        'zlim':None,
                        'xmin':30,
                        'xmax':2e4,
                        'ymin':100,
                        'ymax':300000000
            },
            '080319A':{
                        'z':1.6, #very approximate 
                        'zlim':None,
                        'xmin':100,
                        'xmax':3e3,
                        'ymin':30,
                        'ymax':1000
            },
            '080320':{
                        'z':None, 
                        'zlim':7,
                        'xmin':100,
                        'xmax':2e4,
                        'ymin':10,
                        'ymax':1000
            },
            '061222A':{
                        'z':2.08,
                        'zlim':None,
                        'xmin':100,
                        'xmax':3e3,
                        'ymin':10,
                        'ymax':1000
            },
            '080604':{
                        'z':1.416,
                        'zlim':None,
                        'xmin':100,
                        'xmax':2e4,
                        'ymin':8,
                        'ymax':200
            }
            
        }
        
extinction_grb_list = ['071025','090618','070208','080319C','090709A','051109A','090530','061126','080607','080319B','080319A','080320','061222A','080604']
extinction_list =     ['0.0635','0.0753','0.0140','0.0223', '0.0742', '0.1553', '0.0198','0.1568','0.0191','0.0095', '0.015',  '0.014', '0.099',  '0.05']

# extinction_grb_list = ['071025','090618','070208','080310','080319C','090709A','051109A','090530','061126','080607','080319B','080319A','080320','061222A','080604']
# # average from S & F = Schlafly & Finkbeiner 2011 (ApJ 737, 103) at a nearby calibration star
# extinction_list =     ['0.0635','0.0753','0.0140','0.0336','0.0223', '0.0742', '0.1553', '0.0198','0.1568','0.0191','0.0095', '0.015',  '0.014', '0.099',  '0.05']

fluxdict={"051109A":{"j3":3210.3        ,"h3":4024.4        ,"k3":4737.04        ,"jz3":990.11         ,"hz3":1241.1         ,"kz3":1460.9         },
          "061126" :{"j3":8052.0        ,"h3":12414.        ,"k3":16390.3        ,"jz3":2252.1         ,"hz3":3472.4         ,"kz3":4584.4         },
          "061222A":{"j3":np.nan        ,"h3":np.nan        ,"k3":np.nan         ,"jz3":np.nan         ,"hz3":np.nan         ,"kz3":np.nan         },
          "070208" :{"j3":233.87        ,"h3":324.37        ,"k3":np.nan         ,"jz3":174.54         ,"hz3":242.07         ,"kz3":np.nan         },
          "071025" :{"j3":714.42        ,"h3":959.06        ,"k3":1589.36        ,"jz3":1331.6         ,"hz3":1787.6         ,"kz3":2962.4         },
          # "080310" :{"j3":3928.5        ,"h3":4900.2        ,"k3":5997.86        ,"jz3":2070.6         ,"hz3":2583.1         ,"kz3":3161.9         },
          "080319A":{"j3":np.nan        ,"h3":np.nan        ,"k3":np.nan         ,"jz3":np.nan         ,"hz3":np.nan         ,"kz3":np.nan         },
          "080319B":{"j3":780000        ,"h3":990000        ,"k3":1230000        ,"jz3":106000         ,"hz3":129000         ,"kz3":160000         },
          "080319C":{"j3":1443.0        ,"h3":2114.1        ,"k3":3071.01        ,"jz3":913.57         ,"hz3":1338.4         ,"kz3":1944.2         },
          "080320" :{"j3":np.nan        ,"h3":np.nan        ,"k3":np.nan         ,"jz3":np.nan         ,"hz3":np.nan         ,"kz3":np.nan         },
          "080604" :{"j3":np.nan        ,"h3":np.nan        ,"k3":np.nan         ,"jz3":np.nan         ,"hz3":np.nan         ,"kz3":np.nan         },
          "080607" :{"j3":1856.4        ,"h3":4776.1        ,"k3":12726.5        ,"jz3":286.83         ,"hz3":737.94         ,"kz3":1966.3         },
          "090530" :{"j3":1070.2        ,"h3":1270.0        ,"k3":np.nan         ,"jz3":616.11         ,"hz3":731.10         ,"kz3":np.nan         },
          "090618" :{"j3":23000         ,"h3":23000         ,"k3":np.nan         ,"jz3":11000          ,"hz3":13000          ,"kz3":np.nan         },
          "090709A":{"j3":np.nan        ,"h3":470           ,"k3":1300           ,"jz3":np.nan         ,"hz3":160            ,"kz3":490            }
          }

dustdict={"051109A":{"beta":-0.70 ,"Av":0.00 ,"type":"smc"}, #AV assumed
         "061126" :{"beta":-0.93 ,"Av":0.00 ,"type":"smc"}, #took case where Av=0
         "061222A":{"beta":-0.60 ,"Av":5.00 ,"type":"smc"},
         "070208" :{"beta":-0.60 ,"Av":0.96 ,"type":"smc"},
         # "071025" :{"beta":0.96 ,"Av":1.09 ,"type":"maiolino"},
         "071025" :{"beta":-0.96 ,"Av":0.00 ,"type":"smc"}, #will correct manually later
         # "080310" :{"beta":-0.42 ,"Av":0.19 ,"type":"smc"},
         "080319A":{"beta":-0.60 ,"Av":0.25 ,"type":"smc"},
         "080319B":{"beta":-0.50 ,"Av":0.07 ,"type":"smc"},
         "080319C":{"beta":-0.98 ,"Av":0.59 ,"type":"smc"},
         "080320" :{"beta":np.nan,"Av":np.nan,"type":"smc"},
         "080604" :{"beta":np.nan,"Av":np.nan,"type":"smc"},
         "080607" :{"beta":-0.70 ,"Av":3.30 ,"type":"080607"},
         "090530" :{"beta":-0.60 ,"Av":0.00 ,"type":"smc"}, #beta and AV assumed
         "090618" :{"beta":-0.64 ,"Av":0.24 ,"type":"smc"},
         "090709A":{"beta":-0.70 ,"Av":3.40 ,"type":"smc"}
         }


df_z= pd.DataFrame(meta_dict).T

dfext=pd.DataFrame({'galebv':extinction_list},index=extinction_grb_list)

df_dust=pd.DataFrame(dustdict).T

df_grb = pd.merge(df_z,dfext,how='outer',left_index=True,right_index=True)

df_flux=pd.DataFrame(fluxdict).T

df_grb= pd.merge(df_grb,df_flux,how='outer',left_index=True,right_index=True)

df_grb= pd.merge(df_grb,df_dust,how='outer',left_index=True,right_index=True)

from MiscBin import qObs
jfilt=qObs.J
hfilt=qObs.H
kfilt=qObs.Ks

    

def getdustlaw(dusttype):
    from Modelling import ExtinctModel
    if dusttype == 'smc': return ExtinctModel.avgsmc()
    if dusttype == 'lmc': return ExtinctModel.avglmc()
    if dusttype == 'lmc2': return ExtinctModel.lmc2()
    if dusttype == 'mw': return ExtinctModel.avgmw()
    if dusttype == '080607': return ExtinctModel.fm080607()
    
from Modelling.ExtinctModel import CorrectFluxForGalExt 
from Modelling.ExtinctModel import CorrectLateTimeDust 

# correct J flux for galactic extinction
wavearr = np.zeros(len(df_grb)) + jfilt.wave_A 
j3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.j3,df_grb.j3)
jz3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.jz3,df_grb.jz3)
# correct for HOST galaxy extinction [assumes SMC with current iteration of CorrectLateTimeDust]
# gotta first find what filter to assume for each
wavedf=jfilt.wave_A/(1+df_grb.z)
jz3corrhostlist=[]
for ind in df_grb.index:
    print ind
    try:
        corr_item,ignore = CorrectLateTimeDust(df_grb.loc[ind]['Av'],wavedf[ind],jz3corr[ind],jz3corr[ind],dustlaw=df_grb.type.apply(getdustlaw)[ind])
        jz3corrhostlist.append(float(corr_item))
    except:
        print 'cannot append for index {}'.format(ind)
jz3corrhost=np.array(jz3corrhostlist)
jz3corrhost=pd.Series(jz3corrhost,name='jz3corrhost',index=df_grb.index)
# maiolino: at z=5
# attenuation is 1.14 mag J, 1.15 mag H, 0.99 mag K
# or fators of 2.67, 2.89, 2.48 if you prefur
jz3corrhost.loc['071025']=jz3corrhost.loc['071025']*2.67

offsetslist=[0.5,1.5,-1,0.5,2.5,-1,2.5,1.5,-1,-1,1.5,0.5,0.5,2.5]
offsets=pd.Series(offsetslist,name='offsets',index=df_grb.index)

# correct H flux for galactic extinction
wavearr = np.zeros(len(df_grb)) + hfilt.wave_A 
h3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.h3,df_grb.h3)
hz3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.hz3,df_grb.hz3)
# correct for HOST galaxy extinction [assumes SMC with current iteration of CorrectLateTimeDust]
# gotta first find what filter to assume for each
wavedf=hfilt.wave_A/(1+df_grb.z)
hz3corrhostlist=[]
for ind in df_grb.index:
    print ind
    try:
        corr_item,ignore = CorrectLateTimeDust(df_grb.loc[ind]['Av'],wavedf[ind],hz3corr[ind],hz3corr[ind],dustlaw=df_grb.type.apply(getdustlaw)[ind])
        hz3corrhostlist.append(float(corr_item))
    except:
        print 'cannot append for index {}'.format(ind)
hz3corrhost=np.array(hz3corrhostlist)
hz3corrhost=pd.Series(hz3corrhost,name='hz3corrhost',index=df_grb.index)
# maiolino: at z=5
# attenuation is 1.14 mag J, 1.15 mag H, 0.99 mag K
# or fators of 2.67, 2.89, 2.48 if you prefur
hz3corrhost.loc['071025']=hz3corrhost.loc['071025']*2.89

# correct K flux for galactic extinction
wavearr = np.zeros(len(df_grb)) + kfilt.wave_A 
k3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.k3,df_grb.k3)
kz3corr,ignore = CorrectFluxForGalExt(df_grb.galebv.astype('float'),wavearr,df_grb.kz3,df_grb.kz3)
# correct for HOST galaxy extinction [assumes SMC with current iteration of CorrectLateTimeDust]
# gotta first find what filter to assume for each
wavedf=kfilt.wave_A/(1+df_grb.z)
kz3corrhostlist=[]
for ind in df_grb.index:
    print ind
    try:
        corr_item,ignore = CorrectLateTimeDust(df_grb.loc[ind]['Av'],wavedf[ind],kz3corr[ind],kz3corr[ind],dustlaw=df_grb.type.apply(getdustlaw)[ind])
        kz3corrhostlist.append(float(corr_item))
    except:
        print 'cannot append for index {}'.format(ind)
kz3corrhost=np.array(kz3corrhostlist)
kz3corrhost=pd.Series(kz3corrhost,name='kz3corrhost',index=df_grb.index)
# maiolino: at z=5
# attenuation is 1.14 mag J, 1.15 mag H, 0.99 mag K
# or fators of 2.67, 2.89, 2.48 if you prefur
kz3corrhost.loc['071025']=kz3corrhost.loc['071025']*2.67



# Make brightness dist
df_3min=pd.DataFrame({"j":j3corr,"h":h3corr,"k":k3corr})

# #convert to cgs from microjansky:
# j_3min = arr(j_3min)*10**(-29)
# #convert to AB magnitude:
# j_3min = -2.5*numpy.log10(j_3min) - 48.60
conv2AB = lambda x: -2.5*np.log10(np.array(x)*10**(-29)) - 48.60

df_3min_mag = df_3min.apply(conv2AB)

# df_3min_mag.hist(bins=range(8,19,1),sharex=True,layout=(3,1))


kc = '#CC6677'
hc = '#117733'
jc = '#88CCEE'

#light filled areas
plt.hist(np.array(df_3min_mag.k),bins=np.linspace(df_3min_mag.k.min(),df_3min_mag.k.min()+10,10),histtype='stepfilled',alpha=0.1,color=kc,linewidth=0)
plt.hist(np.array(df_3min_mag.h),bins=np.linspace(df_3min_mag.h.min(),df_3min_mag.h.min()+10,10),histtype='stepfilled',alpha=0.1,color=hc,linewidth=0)
plt.hist(np.array(df_3min_mag.j),bins=np.linspace(df_3min_mag.j.min(),df_3min_mag.j.min()+10,10),histtype='stepfilled',alpha=0.1,color=jc,linewidth=0)

plt.hist(np.array(df_3min_mag.k),bins=np.linspace(df_3min_mag.k.min(),df_3min_mag.k.min()+10,10),histtype='step',color=kc,linewidth=3,label="$K$",linestyle='dashdot')
plt.hist(np.array(df_3min_mag.h),bins=np.linspace(df_3min_mag.h.min(),df_3min_mag.h.min()+10,10),histtype='step',color=hc,linewidth=3,label="$H$",linestyle='dashed')
plt.hist(np.array(df_3min_mag.j),bins=np.linspace(df_3min_mag.j.min(),df_3min_mag.j.min()+10,10),histtype='step',color=jc,linewidth=3,label="$J$")

plt.xlabel('Apparent Magnitude (AB)', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.title('$t=3$ min (observer frame)')
plt.yticks(arr([0, 1., 2., 3., 4.]))
ax = plt.gca()
# ax.set_xlim(ax.get_xlim()[::-1]) # reversing the xlimits
ax.set_xlim((19.0,8.0))
ax.set_ylim((0,3.5))
plt.legend(fontsize=18,frameon=False)



def plot_lum_rest():
    '''f_{rest,V} = f_{rest_corr}*[nu_V/ ((1+z)nu_J)]^beta for  flux \propto nu^beta and beta negative values'''
    clf()
    
        # 
        # those are values from the lightcurve, where I have 3minutes*(1+z) = t_0, and then I get the value of of the lightcurve at t_0
        # 
        # me: right but what function
        # 
        # Pierre: most of them are interpolations from Dan Perley's fits
        # 10:07 AM 
        # some actually have values at times that are very close to 3 minutes *(1+z) and for those I just use the actual data points
    
    # f_rest_corr = [2252.14, 1626.48, 403.717, 11783.2, 913.329, 549.616, 286.863, 990.110, 14.7689, 174.540, 1419.79, 149309.80115] 
    
    # beta = [-1.35, -0.8, -0.96, -0.22, -1.73, -0.84, -3.48, -0.42, -3.81, -0.3, -1.7, -0.47]
    beta = df_grb.beta
    # z_list_limits = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8, 0.9382]
    z_list_limits = df_grb.z
    
    arrb = arr(beta)
    arrz = arr(z_list_limits)
#    nu_V = 5.444646098003629764065335753176043557e+14 
#    nu_J = 2.398339664e+14
    nu_V = qObs.V.freq
    nu_J = qObs.J.freq
    nu_K = qObs.Ks.freq
    nu_H = qObs.H.freq
    
    abs_mag=[]
    for nu_obs, f_obs in [(nu_J,jz3corrhost),(nu_H,hz3corrhost),(nu_K,kz3corrhost)]:
        arrf = arr(f_obs)
        f_rest_V = arrf * (nu_V/ ((1+arrz)*nu_obs))**(beta)
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
        F_10pc = F_10pc.fillna(-1).astype('float')   # type problems..  
        Absol_Mag = -2.5*numpy.log10(F_10pc) - 48.60    #Absolute mag in AB mag
        abs_mag.append(Absol_Mag)

    plt.hist(np.array(abs_mag[1]),bins=np.linspace(-33.8,-25.5,9),histtype='stepfilled',color=hc,linewidth=0,alpha=0.15)
    plt.hist(np.array(abs_mag[0]),bins=np.linspace(-33.8,-25.5,9),histtype='stepfilled',color=jc,linewidth=0,alpha=0.15)
    plt.hist(np.array(abs_mag[2]),bins=np.linspace(-33.8,-25.5,9),histtype='stepfilled',color=kc,linewidth=0,alpha=0.15)
    
    plt.hist(np.array(abs_mag[1]),bins=np.linspace(-33.8,-25.5,9),histtype='step',color=hc,linewidth=11,label="$H$ detection")
    plt.hist(np.array(abs_mag[0]),bins=np.linspace(-33.8,-25.5,9),histtype='step',color=jc,linewidth=7,label="$J$ detection")
    plt.hist(np.array(abs_mag[2]),bins=np.linspace(-33.8,-25.5,9),histtype='step',color=kc,linewidth=3,label="$K$ detection")
    switch=-0.2
    for ind in Absol_Mag.index[::-1]:
        # plt.text(Absol_Mag.loc[ind],1.0,ind,rotation=90,ha='center',va='center',fontsize=12)
        # plt.annotate(ind,xy=(abs_mag[0].loc[ind],0),xytext=(abs_mag[0].loc[ind],0.6-switch),color=jc,arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color=jc),rotation=90,ha='center',va='center',fontsize=14)
        # plt.annotate(ind,xy=(abs_mag[1].loc[ind],0),xytext=(abs_mag[1].loc[ind],1.6-switch),color=hc,arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color=hc),rotation=90,ha='center',va='center',fontsize=14)
        # plt.annotate(ind,xy=(abs_mag[2].loc[ind],0),xytext=(abs_mag[2].loc[ind],2.6-switch),color=kc,arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color=kc),rotation=90,ha='center',va='center',fontsize=14)
        # 
        plt.annotate(ind,xy=(abs_mag[1].loc[ind],0),xytext=(abs_mag[1].loc[ind],offsets[ind]),color='dimgray',arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color='gray',lw=2),rotation=90,ha='center',va='center',bbox=dict(boxstyle="round", fc="whitesmoke",ec="gray",lw=2),fontsize=16)
        switch*=-1
    # hist(Absol_Mag,6,color="#88CCEE")
    
    
    
    plt.xlabel('Absolute Magnitude $M_V$', fontsize=20)
    plt.ylabel('Number', fontsize=20)
    plt.yticks(arr([0, 1., 2., 3., 4.,5.,6.]))
    ax = matplotlib.pyplot.gca()
    ax.set_xlim(ax.get_xlim()[::-1]) # reversing the xlimits
    ax.set_xlim((-24,-35))
    ax.set_ylim((0,3.5))
    plt.title('$t=3$ min (rest frame)',fontsize=20)
    plt.legend(loc='upper left',fontsize=18,frameon=False)
    
    plt.savefig('Lum_dist_rest.eps')
    plt.savefig('Lum_dist_rest.png')

    print 'Done'
    
    return abs_mag


def get_dict(burstlist):
    #Burstlist = update_z_all_GRBs(GRB_list)
    GRB_list =  ['061126', '080310', '080330', '090618', '080319C', '090530', '080607', '051109A', '070208', '071025']
    z_list_bad = [1.1588, 2.4274, 1.51, 0.54, 1.95, 1.6, 3.036, 2.346, 3.5, 1.165, 4.8, 0, 0]
    Burstlist = makedict(GRB_list,outputtext=True)
    return Burstlist

vgpp="/Volumes/MyPassport/PTELBACKUP2/picklefiles/very_good_pickles/"
vgpp="/Volumes/GRB/picklefiles/very_good_pickles/"
def makedict(df_grb, outputtext=False,
    very_good_pickle_path=vgpp):
    ''' make dict for all_bursts() '''
    all_GRB_dict = {}
    
    GRB_list = list(df_grb.index)
    cc=0
    for index, GRB in enumerate(GRB_list):
        
        globstr = very_good_pickle_path + GRB + '*'
        globresult = glob.glob(globstr)
        if len (globresult) == 0:
            print "Can't Find any files in {}. External HD plugged in?".format(globresult)
            raise(ValueError)
        pathstr = globresult[0]
        if not os.path.exists(pathstr):
            print "The path {} doesn't appear to exist. Wrong path?".format(pathstr)
        print pathstr
        result = qPickle.load(pathstr)
        galebv=df_grb.loc[GRB]['galebv']
        redshift=df_grb.loc[GRB]['z']
        if outputtext:      
            
            try:
                 q_phot.textoutput(result,name=GRB,galebv=galebv,redshift=redshift)
            except:
                 print "CANNOT MAKE OUTPUTTEXT FOR {}".format(GRB)
                 cc+=1
            if cc==3: raise ValueError
        GRB_dict = {GRB:result}
        all_GRB_dict.update(GRB_dict)
    return all_GRB_dict
    
# MAKE PARSEABLE PHOTOMETRY FILES
def make_outdict():
    outdict = makedict(df_grb, outputtext=True)
    return outdict
    
# TRANSPORT TABLES 
#todo

# load tables
catdir = "/Users/amorgan/qsoft/store/ptelcat/"
phottables=glob.glob(catdir+"*data.txt")

def _build_up_extraptime(objblock):
    # FIXME. Add the desired extrapolation times. like 3 min in rest frame.
    extraptime = [300.]
    return extraptime
    
def obj_loop(do_interp=False):
    for table in phottables:
        phot_obj=PhotParse.PhotParse(table)
        
        GRB_key = phot_obj.name.split(',')[-1] # eg 071025
        GRB_df = df_grb.loc[GRB_key]
        xlims = (GRB_df['xmin'],GRB_df['xmax']) # grab plot imits from DF
        ylims = (GRB_df['ymin'],GRB_df['ymax'])
        phot_obj.PlotLC(show=False,legend=False,xlimits=xlims,ylimits=ylims,figsize=(9,7))
        phot_obj.WriteCondensedTable()
    
        if do_interp:
            #build interpolation!
            extraptime=_build_up_extraptime(phot_obj)
            interpolate_list = ['PAIRITEL_J','PAIRITEL_H','PAIRITEL_K']
            try:
                objblock_interpolated=DustModel.BuildInterpolation(phot_obj,extraptime,interpolate_list,
                take_as_given_list=None,interp_type='smart',plot=True,value_lims=(18.5,10.5),error_lims=(0.00,0.30))
            except:
                print "UH OH. Can't do spline fit."
            #     continue