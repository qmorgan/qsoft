# Adapted from http://www.scipy.org/Cookbook/Matplotlib/Interactive_Plotting
# on 1/26/10

import math

import pylab
import matplotlib
from pylab import scatter
from pylab import connect
from pylab import subplot

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    if xtol is None:
      xtol = ((max(xdata) - min(xdata))/float(len(xdata)))*2#/2 amorgan increased
    if ytol is None:
      ytol = ((max(ydata) - min(ydata))/float(len(ydata)))*2#/2
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = pylab.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """
    return the distance between two points
    """
    return math.hypot(x1 - x2, y1 - y2)

  def __call__(self, event):
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      if self.axis is None or self.axis==event.inaxes:
        annotes = []
        for x,y,a in self.data:
          if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
            annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote)
          for l in self.links:
            l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote,include_vals=True):
    """
    Draw the annotation on the plot
    """
    if (x,y) in self.drawnAnnotations:
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      if include_vals:
        t = axis.text(x,y, "(%3.2f, %3.2f) - %s"%(x,y,annote), )
      else:
        t = axis.text(x,y," - %s"%(annote), )  
      m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)


def test():
    x = range(10)
    y = range(10)
    annotes = ['a', 'b', 'c', 'd', 'GRB1010101', 'f', 'g', 'h', 'i', 'j']

    scatter(x,y)
    af =  AnnoteFinder(x,y, annotes)
    connect('button_press_event', af)
    
def testcolor():
    x = range(10)
    y = range(10)
    z = range(10)
    annotes = ['a', 'b', 'c', 'd', 'GRB1010101', 'f', 'g', 'h', 'i', 'j']

    from Plotting import ColorScatter
    ColorScatter.ColorScatter(x,y,z,cmap='cool')
    af =  AnnoteFinder(x,y, annotes)
    connect('button_press_event', af)
    
    
def AnnotatedSubPlot(xlist=[range(3),[1,3],range(3)],ylist=[range(3),[1,3],range(3)], 
            annotelist = [['a','b','c'],['a','d'],['a','b','d']],logx=False,logy=False,
            ynames=None,xnames=None, zlist=None):

    '''As long as the annotations are the same for each plot, then they will be 
    linked among the plots.
    '''
    
    def linkAnnotationFinders(afs):
      for i in range(len(afs)):
        allButSelfAfs = afs[:i]+afs[i+1:]
        afs[i].links.extend(allButSelfAfs)

    
    
    if len(xlist) != len(ylist):
        raise(ValueError('Length of lists does not match'))
    
    if zlist:
        from Plotting import ColorScatter
        if len(zlist) != len(ylist):
            raise(ValueError('Length of lists does not match'))
    
    if len(xlist) == 1: numsubplots = '11'
    if len(xlist) == 2: numsubplots = '21'
    if len(xlist) == 3: numsubplots = '31'
    if len(xlist) == 4: numsubplots = '22'
    if len(xlist) == 5: numsubplots = '23'
    if len(xlist) == 6: numsubplots = '23'
    if len(xlist) > 6: numsubplots = '33'
    if len(xlist) > 9:
        print 'Too many plots!'
        raise(ValueError('Too Many Plots!'))
    
    ind = 0
    aflist = []
    while ind < len(xlist):
        plind = str(ind + 1)
        subplotkey = numsubplots + plind
        pl = subplot(subplotkey)
        xx = xlist[ind]
        yy = ylist[ind]            
        annotes = annotelist[ind]
        if zlist:
            # if z dimension is included, add as a "color" dimension
            zz = zlist[ind]
            pl.scatter(xx,yy,c=zz,cmap='spring')
        else:
            pl.scatter(xx,yy)
        if logx:
            pl.set_xscale('log')
            pl.set_xlim(1e-6,1e6) # CHANGE TO GUESS AUTOMATICALLY
        if logy:
            pl.set_yscale('log')
            pl.set_ylim(1e-6,1e6) # CHANGE TO GUESS AUTOMATICALLY
        if xnames: pylab.xlabel(xnames[ind])
        if ynames: pylab.ylabel(ynames[ind])
        if len(annotes) != len(xx): print 'Annotes dont match!'
        if len(xx) != len(yy): print 'xx and yy dont match!'
        af = AnnoteFinder(xx,yy,annotes)
        aflist.append(af)
        connect('button_press_event',af)
        ind += 1
        
    # subplot(121)
    # scatter(x,y)
    # af1 = AnnoteFinder(x,y, annotes)
    # connect('button_press_event', af1)
    # 
    # subplot(122)
    # scatter(x2,y2)
    # af2 = AnnoteFinder(x2,y2, annotes)
    # connect('button_press_event', af2)

    linkAnnotationFinders(aflist)
