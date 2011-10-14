# Adapted from http://www.scipy.org/Cookbook/Matplotlib/MultilinePlots

import math
from pylab import figure, show, setp, hist
from numpy import sin, cos, exp, pi, arange
from Plotting.ColorScatter import ColorScatter as scatter
import numpy
from Plotting.q_hist import histOutline

# Draw a custom axis for the colorbar: 
# plt.colorbar(sc,cax=plt.gcf().add_axes((0.90,0.1,0.02,0.8)))


def GridPlot(data,fig=None,zdata=None,labels=None,no_tick_labels=False,hist=True, 
            histbins=None, histrangelist=None, histloc='tl',color='black',
            colorbar=None, show=True, noalpha=False, **kwargs):
    '''hist text: tr, br, tl, bl, ce
    noalpha means dont use alpha for the boxes.
    '''
    
    histfill = True
    
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
        
    if not fig:
        fig = figure()
        oldfig = None
    else:
        oldfig = fig
    
    n_datapoints = len(matr[0])
    point_size = numpy.ceil(500/((n_datapoints)*N)) + 3
    if point_size < 1:
        point_size = 4

    horiz_width = 1.0 - right_buffer - left_buffer
    vert_width = 1.0 - top_buffer - bottom_buffer

    hw = horiz_width/N # subplot width
    vw = vert_width/N 
    
    if oldfig:
        old_axes = oldfig.get_axes()
        
    count = 0
    for ii in arange(N):
        for jj in arange(N):
            
            qq = left_buffer + ii*hw
            rr = ((1-top_buffer) - vw) - jj*vw
            
            if oldfig:
                ax1 = old_axes[count]
                old_xlims = ax1.get_xlim()
                old_ylims = ax1.get_ylim()
            else:
                ax1 = fig.add_axes([qq,rr,hw,vw])
            
            # ax1.scatter(matr[ii],matr[jj])
            if hist and ii==jj:
                if not oldfig:
                    ax1.scatter(matr[ii],matr[jj],s=0) #plotting blank to get x_lim correct
                    xlims = ax1.get_xlim()
                
                if not histrangelist or len(histrangelist) != N:
                    histrange = (numpy.nanmin(matr[ii]),numpy.nanmax(matr[ii]))
                    if len(histrangelist) != N: 
                        print 'Histrange list is not an N lenght list of tuples. Not using.'
                elif histrangelist:
                    histrange = histrangelist[ii]
                # get histogram
                if not histbins:
                    histbins = numpy.ceil(n_datapoints/4)
                
                (bins,nn) = histOutline(matr[ii],bins=histbins,range=histrange)
                if histfill:
                    ax1.hist(matr[ii],bins=histbins,range=histrange,histtype='stepfilled',color=color)
                else:    
                    ax1.plot(bins,nn,'k-',color=color)
                
                if not oldfig:
                    ax1.set_ylim((0,max(nn)))
                    ax1.set_xlim(xlims)
                else:
                    ax1.set_ylim(old_ylims)
                    ax1.set_xlim(old_xlims)
                if histloc:
                    pos_list = ['tl','bl','tr','br','ce']
                    if histloc not in pos_list:
                        print '%s not acceptable: see %s' % (histloc, str(pos_list))
                        histloc = raw_input('New histloc: ')
                    if histloc[0] == 't':
                        va = 'top'
                        vertloc = 0.85
                    elif histloc[0] == 'b': 
                        va = 'bottom'
                        vertloc = 0.15
                    elif histloc[0] == 'c':
                        va = 'center'
                        vertloc = 0.5
                        horiloc = 0.5
                    if histloc[1] == 'l':
                        ha = 'left'
                        horiloc = 0.15
                    elif histloc[1] == 'r':
                        ha = 'right'
                        horiloc = 0.85
                    
                    histtext = '''N=%i\nMax=%i''' % (numpy.sum(nn)/2,numpy.max(nn))
                    
                    if noalpha: 
                        bbox = dict(color='black',edgecolor='black',facecolor='white')
                    else:
                        bbox=dict(color='black',edgecolor='black',facecolor='white', alpha=0.7)
                    
                    ax1.text(horiloc,vertloc,histtext,color=color,fontsize=8,
                        horizontalalignment=ha, verticalalignment=va,
                        transform=ax1.transAxes, bbox=bbox)
                            
            elif ii != jj:    
                scatter(matr[ii],matr[jj],z=zdata,axis=ax1,s=point_size,color=color,colorbar=colorbar, **kwargs)
            else:
                pass
            
            if jj == N-1 and labels:
                ax1.set_xlabel(labels[ii])
        
            if ii == 0 and labels:
                ax1.set_ylabel(labels[jj])
            
            if no_tick_labels and ii == 0:
                setp(ax1.get_yticklabels(), visible=False)
            if ii != 0:
                setp(ax1.get_yticklabels(), visible=False)
            
            if oldfig:
                ax1.set_xlim(old_xlims)
                ax1.set_ylim(old_ylims)
                
            # turn off x ticklabels for all but the lower axes
            if no_tick_labels and jj == N-1:
                setp(ax1.get_xticklabels(), visible=False)
            if jj != N-1:
                setp(ax1.get_xticklabels(), visible=False)
            count += 1
    return fig

def TestPlot(fig=None):
    A = numpy.array([1,2,3,4,2,5,8,3,2,3,5,6])
    B = numpy.array([8,7,3,6,4,numpy.nan,9,3,7,numpy.nan,2,4])
    C = numpy.array([6,3,4,7,2,1,1,7,8,4,3,2])
    D = numpy.array([5,2,4,5,3,8,2,5,3,5,6,8])
    
    # A work around to get the histograms overplotted with each other to overlap correctly;
    histrangelist = [(numpy.nanmin(A),numpy.nanmax(A)),(numpy.nanmin(B),numpy.nanmax(B)),
                (numpy.nanmin(C),numpy.nanmax(C)),(numpy.nanmin(D),numpy.nanmax(D))]
    
    data = numpy.array([A,B,C,D])
    labels = ['A','3','C','D']

    fig = GridPlot(data,labels=labels, no_tick_labels=True, color='black', 
                    hist=True, histbins=3, histloc='tl', histrangelist=histrangelist, fig=None) 
    
    # Data of note to plot in different color
    A2 = numpy.array([1,2,3,4])
    B2 = numpy.array([8,7,3,6])
    C2 = numpy.array([6,3,4,7])
    D2 = numpy.array([5,2,4,5])
    data2 = numpy.array([A2,B2,C2,D2])
    
    fig = GridPlot(data2,labels=labels, no_tick_labels=True, color='red', 
                hist=True, histbins=3, histloc='tr', histrangelist=histrangelist, fig=fig) 
    
    return fig