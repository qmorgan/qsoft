# Adapted from http://www.scipy.org/Cookbook/Matplotlib/MultilinePlots

import math
from pylab import figure, show, setp
from numpy import sin, cos, exp, pi, arange
import numpy

# Draw a custom axis for the colorbar: 
# plt.colorbar(sc,cax=plt.gcf().add_axes((0.90,0.1,0.02,0.8)))


N = 3

A = numpy.array([1,2,3,4])
B = numpy.array([2,7,3,6])
C = numpy.array([6,3,4,1])

matr = numpy.array([A,B,C])
labels = ['A','B','C']

fig = figure()
t = arange(0.0, 2.0, 0.01)

## EDITABLE
right_buffer = 0.1
bottom_buffer = 0.1
## end editable

left_buffer = 0.1
top_buffer = 0.1
horiz_width = 1.0 - right_buffer - left_buffer
vert_width = 1.0 - top_buffer - bottom_buffer

hw = horiz_width/N # subplot width
vw = vert_width/N 

count = 0
for ii in arange(N):
    for jj in arange(N):
        qq = 0.1 + ii*hw
        rr = (0.9 - vw) - jj*vw
        ax1 = fig.add_axes([qq,rr,hw,vw])
        ax1.plot(matr[ii],matr[jj])
        
        if jj == N-1:
            ax1.set_xlabel(labels[ii])
        
        if ii == 0:
            ax1.set_ylabel(labels[jj])
            
        if ii != 0:
            setp(ax1.get_yticklabels(), visible=False)
                
        # turn off x ticklabels for all but the lower axes
        if jj != N-1:
            setp(ax1.get_xticklabels(), visible=False)
        
            

# ax1 =fig.add_axes([0.1, 0.7, 0.8, 0.2], **axprops)
# ax1.plot(t, s1)
# ax1.set_ylabel('S1', **yprops)
# 
# axprops['sharex'] = ax1
# axprops['sharey'] = ax1
# # force x axes to remain in register, even with toolbar navigation
# ax2 = fig.add_axes([0.1, 0.5, 0.3, 0.2], **axprops)
# 
# 
# 
# ax2.plot(t, s2)
# ax2.set_ylabel('S2', **yprops)
# 
# ax3 = fig.add_axes([0.1, 0.3, 0.8, 0.2], **axprops)
# ax3.plot(t, s4)
# ax3.set_ylabel('S3', **yprops)
# 
# ax4 = fig.add_axes([0.1, 0.1, 0.8, 0.2], **axprops)
# ax4.plot(t, s4)
# ax4.set_ylabel('S4', **yprops)
# 
# # turn off x ticklabels for all but the lower axes
# for ax in ax1, ax2, ax3:
#     setp(ax.get_xticklabels(), visible=False)

show()