import pylab
import scipy
import numpy
import matplotlib
from numpy import *
from scipy import *
from pylab import *
import scipy.interpolate as interpolate

def Bumps(x,y,z=None):
    pdict =  {'a':[3,2,1],'b':[2,3,2],'c':[1,1,3]}
    x = arange(len(pdict.values()[0]))
    for element in plist:
       pylab.plot(x,element)
    

def ColorScatterExample():
    x=list(arange(20))
    y=list(arange(20))
    z=scipy.rand(20)*50
    ColorScatter(x,y,z,zlabel='Test Label',cmap='jet',colorbar=True,xjitter=2)


def ColorScatter(x,y,z=None,zlabel=None,cmap='jet',colorbar=True,discrete=0,yjitter=None,\
    xjitter=None,marker='o',retjitter=False, **kwargs):
    '''set discrete to N for splitting up into N values
    yjitter sets a percent random jitter in the y direction to help distinguish
    overlapping values.  
    
    If retjitter=True, return a tuple of the jitter arrays.
    
    The jitter parameters can be either a number, to specify a randomly generated
    jitter (multiplied by the  value of the jitter parameter), or it can be an
    array of the length of the axis of pre-defined values to bump the array by.
    This is useful if you want to overplot the same jitter values twice.  
    
    E.g., if you have some special values that you want to overplot a ring around,
    separate the array into special and non-special values.  Plot the non-special
    values first as a normal scatter plot, then the special values as a normal 
    scatter plot (but save the jitter values with retjitter), and then overplot a 
    ring around the same values.
    
    COLORBAR STUFF: if you want to include a single colorbar for multiple sets
    of data on a single plot; set the vmin and vmax kwargs for each set of data
    to be the same.  This will set the range of acceptible colors. 
    
    In [158]: a1 = [1,3,5]
    In [159]: b1 = [1,3,5]
    In [160]: a2 = [2,4,6]
    In [161]: b2 = [2,4,6]
    In [162]: z1 = [1,5,10]
    In [163]: z2 = [15,20,25]
    In [164]: sc = pylab.scatter(a1,b1,c=z1,cmap=cmap,vmin=0,vmax=25)
    In [165]: sc = pylab.scatter(a2,b2,c=z2,cmap=cmap,vmin=0,vmax=25)
    In [166]: pylab.colorbar(sc)
    
    '''
        
    if yjitter:
        try:
            list(yjitter) # check if a list or array
            yjitter = numpy.array(yjitter)
            assert(len(yjitter)==len(y))
            y = y + yjitter
        except:
            y_lim_len = pylab.ylim()[1] - pylab.ylim()[0]
            yjitter = (pylab.random(len(y))-0.5)*yjitter*y_lim_len
            y = y + yjitter
    else:
        yjitter = pylab.zeros(len(y))
        
    if xjitter:
        try:
            list(xjitter) # check if a list or array
            xjitter = numpy.array(xjitter)
            assert(len(xjitter)==len(x))
            x = x + xjitter
        except:
            x_lim_len = pylab.xlim()[1] - pylab.xlim()[0]
            xjitter = (pylab.random(len(x))-0.5)*xjitter*x_lim_len
            x = x + xjitter
    else:
        xjitter = pylab.zeros(len(x))        
    
    if z != None:
        z = array(z)
        if discrete:
            try:
                N = int(discrete)
            except:
                N = 2
            cmap = cmap_discretize(pylab.cm.cool_r,N)
        sc = pylab.scatter(x,y,c=z,cmap=cmap,marker=marker,**kwargs)
        if colorbar: 
            cb = pylab.colorbar(sc)
            if zlabel:
                cb.set_label(zlabel)
    else:
        sc = pylab.scatter(x,y,marker=marker,**kwargs)
        
    if retjitter:
        return((list(xjitter),list(yjitter)))

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = linspace(0,1.,N)
    # N+1 indices
    indices = linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def GridDataExample():
    import numpy as np
    from matplotlib.mlab import griddata
    import matplotlib.pyplot as plt
    import numpy.ma as ma
    from numpy.random import uniform
    # This was stolen from somewhere
    # make up some randomly distributed data
    npts = 200
    x = uniform(-2,2,npts)
    y = uniform(-2,2,npts)
    z = x*np.exp(-x**2-y**2)
    # define grid.
    xi = np.linspace(-2.1,2.1,100)
    yi = np.linspace(-2.1,2.1,100)
    # grid the data.
    zi = griddata(x,y,z,xi,yi)
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(-2,2)
    plt.ylim(-2,2)
    plt.title('griddata test (%d points)' % npts)
    plt.show()
