# Adapted from http://www.scipy.org/Cookbook/Matplotlib/MultilinePlots

import math
from pylab import figure, show, setp
from numpy import sin, cos, exp, pi, arange
from Plotting.ColorScatter import ColorScatter as scatter
import numpy

# Draw a custom axis for the colorbar: 
# plt.colorbar(sc,cax=plt.gcf().add_axes((0.90,0.1,0.02,0.8)))


def GridPlot(data,labels=None,no_tick_labels=False):
    
    ## EDITABLE
    right_buffer = 0.1
    bottom_buffer = 0.1
    left_buffer = 0.1
    top_buffer = 0.1
    ## end editable
    
    matr = numpy.array(data)
    
    N = len(matr)
    if labels and len(labels) != N:
        print '''
        Length of labels array (%i) != length of data array (%i);
        not including labels.''' % (len(labels),N)
        labels = None

    fig = figure()

    point_size = numpy.ceil(500/(len(matr[0])*N))
    if point_size < 1:
        point_size = 1

    horiz_width = 1.0 - right_buffer - left_buffer
    vert_width = 1.0 - top_buffer - bottom_buffer

    hw = horiz_width/N # subplot width
    vw = vert_width/N 

    count = 0
    for ii in arange(N):
        for jj in arange(N):
            qq = left_buffer + ii*hw
            rr = ((1-top_buffer) - vw) - jj*vw
            ax1 = fig.add_axes([qq,rr,hw,vw])
            #ax1.scatter(matr[ii],matr[jj])
            scatter(matr[ii],matr[jj],s=point_size)
        
            if jj == N-1 and labels:
                ax1.set_xlabel(labels[ii])
        
            if ii == 0 and labels:
                ax1.set_ylabel(labels[jj])
            
            if no_tick_labels and ii == 0:
                setp(ax1.get_yticklabels(), visible=False)
            if ii != 0:
                setp(ax1.get_yticklabels(), visible=False)
                
            # turn off x ticklabels for all but the lower axes
            if no_tick_labels and jj == N-1:
                setp(ax1.get_xticklabels(), visible=False)
            if jj != N-1:
                setp(ax1.get_xticklabels(), visible=False)
        
            

    show()

def TestPlot():
    A = numpy.array([1,2,3,4])
    B = numpy.array([2,7,3,6])
    C = numpy.array([6,3,4,1])
    D = numpy.array([5,2,4,5])

    data = numpy.array([A,B,C,D])
    labels = ['A','3','C','D']

    GridPlot(data,labels=labels, no_tick_labels=True) 