#schaeferparse.py
# myline='080516A >  3.15  0.5 55 >  3.53  1.24  3.42  1.24 122,123  >  3.55 0.3'

# Title: Redshift Catalog for Swift Long GRBs 
# Authors: Xiao L., Schaefer B.E. 
# Table: Our Redshifts and Spectroscopic Redshifts
# ================================================================================
# Byte-by-byte Description of file: apj385101t5_mrt.txt
# --------------------------------------------------------------------------------
#    Bytes Format Units   Label   Explanations
# --------------------------------------------------------------------------------
#    1-  7 A7     ---     GRB     GRB identification
#        9 A1     ---   l_zspec   Limit flag on zspec
#   11- 16 F6.3   ---     zspec   ? Spectroscopic redshift
#   18- 20 F3.1   ---   e_zspec   ? Uncertainty in zspec 
#   22- 23 A2     ---   r_zspec   Reference for zspec (1)
#       25 A1     ---   l_zbest   Limit flag on zbest 
#   27- 32 F6.3   ---     zbest   ? Best estimated redshift
#   34- 37 F4.2   ---   e_zbest   ? The 1{sigma} lower limit on zbest
#   39- 44 F6.3   ---   E_zbest   ? The 1{sigma} upper limit on zbest
#   46- 49 F4.2   ---     zphot   ? Upper limit from photometric redshift
#   51- 58 A8     ---   r_zphot   Reference for zphot (1)
#       60 A1     ---   l_zfinal  Limit flag on zfineal
#   62- 66 F5.2   ---     zfinal  Final redshift
#   68- 70 F3.1   ---   e_zfinal  ? Uncertainty in zfinal 
# --------------------------------------------------------------------------------

# In [43]: myline[0:7]
# Out[43]: '080516A'
# 
# In [44]: myline[8:9]
# Out[44]: '>'
# 
# In [45]: myline[8:10]
# Out[45]: '> '
# 
# In [46]: fmt
# Out[46]: '%7s%2s%6s%5s%3s%2s%6s%6s%6s%6s%8s%3s%6s%4s'
# 
# In [47]: myline[11:17]
# Out[47]: '3.15  '
# 
# In [48]: myline[18:23]
# Out[48]: '.5 55'
# 
# In [49]: myline[17:23]
# Out[49]: '0.5 55'
# 
# In [50]: myline[0:7]
# Out[50]: '080516A'
# 
# In [51]: myline[7:9]
# Out[51]: ' >'
# 
# In [52]: myline[9:15]
# Out[52]: '  3.15'
# 
# In [53]: myline[15:20]
# Out[53]: '  0.5'
# 
# In [54]: myline[20:23]
# Out[54]: ' 55'
# 
# In [55]: myline[23:25]
# Out[55]: ' >'
# 
# In [56]: myline[25:31]
# Out[56]: '  3.53'
# 
# In [57]: myline[31:37]
# Out[57]: '  1.24'
# 
# In [58]: myline[37:43]
# Out[58]: '  3.42'

import numpy
import pylab
from RedshiftMachine import LoadDB

f=file('schaefer11_nohead.txt','r')
lines = f.readlines()
mydict = {}
speczlist=[]
zbestlist=[]
qtrainlist=[]

mydb = LoadDB.LoadDB('GRB_short_removed')

for line in lines:
    grbname=line[0:7].strip()
    specislimit = line[8]
    spec_z = line[10:15]
    
    if grbname in mydb.dict:
        if 'Q_hat_train' in mydb.dict[grbname]:
            q_hat_train = mydb.dict[grbname]['Q_hat_train']
        else:
            print 'No q_hat_train for %s' % grbname
            q_hat_train = numpy.nan
    else:
        print '%s not in my dictionary' % grbname
        q_hat_train = numpy.nan
        
    try:
        spec_z = float(spec_z)
    except:
        spec_z = numpy.nan
    zbestislimit = line[24]
    z_best = line[26:31]
    try:
        z_best = float(z_best)
    except:
        z_best = numpy.nan
    
    speczlist.append(spec_z)
    zbestlist.append(z_best)
    qtrainlist.append(q_hat_train)
    mydict.update({grbname:{'q_hat_train':q_hat_train,'spec_z':spec_z,'spec_z_islim':specislimit,'z_best':z_best,'zbest_islim':zbestislimit}})
    

# pylab.scatter(zbestlist,speczlist)
# pylab.xlabel("Pseudo-Z")
# pylab.ylabel("Spectroscopic Redshift")

from matplotlib.ticker import FuncFormatter
from MiscBin import q

def ff(x,pos=None):
    if x < -1:
        return "%.2f" % (10**x)
    elif x < 0:
        return "%.1f" % (10**x)
    elif 10**x == 8.5:
        return "%.1f" % (10**x)
    else:
        return "%i" % (10**x)

ax=pylab.subplot(111)
formatter = FuncFormatter(ff)

ax.hist(q.RemoveNaN(numpy.log10(numpy.array(zbestlist))),bins=15)

ax.xaxis.set_major_formatter(formatter)
ax.set_xticks([-2,-1,numpy.log10(0.3),0,numpy.log10(2),numpy.log10(3),numpy.log10(4),numpy.log10(6),numpy.log10(12)])
ax.set_xlim((numpy.log10(0.005),numpy.log10(13)))
ax.set_ylabel("Number")
ax.set_xlabel("Pseudo-z")

aa=q.RemoveNaN(numpy.array(zbestlist))
fracgt4 = float(len(aa[aa>4]))/float(len(aa))
fracgt8 = float(len(aa[aa>8]))/float(len(aa))
print '%f percent of bursts have a pseudo-z of > 4.0' % (fracgt4*100.)
print '%f percent of bursts have a pseudo-z of > 8.0' % (fracgt8*100.)