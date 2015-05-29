import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os
from matplotlib.patches import CirclePolygon
#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot

class ternaryfaces:
    def __init__(self, ax, ellabels=['A', 'B', 'C', 'D'], offset=None, nintervals=10., outlinealpha=0.4):
        self.outlinealpha=outlinealpha
        self.ternaryplot=TernaryPlot(ax, outline=False)
        self.ax=ax
        
        self.ax.set_xlim(-.1, 2.6)
        self.ax.set_ylim(-.1, 3.**.5/2+.1)
        self.cartendpts=numpy.float32([[0, 0], [.5, numpy.sqrt(3.)/2.], [1, 0]])
        self.ellabels=ellabels
        self.outline()
        self.nint=1.*nintervals
        self.delta=1./self.nint
        self.patch_xyc=lambda x, y, c, **kwargs:self.ax.add_patch(CirclePolygon((x, y),radius=self.delta/3.**.5,resolution=6, color=c, **kwargs))
        if offset is None:
            self.offset=self.delta
    
    def xy_skipind(self, x, y, skipind):
        x+=([0.5, 1, 1.5, 0.][skipind])
        if skipind%2==0:
            y=-1.*y+3.**.5/2
        return x, y
        
    def outline(self):
        for skipind in range(4):
            skipfirstline=skipind!=3
            for i, ep in enumerate(self.cartendpts):
                for ep2 in self.cartendpts[i+1:]:
                    if skipfirstline:
                        skipfirstline=False
                        continue
                    x, y=self.xy_skipind(numpy.array([ep[0], ep2[0]]), numpy.array([ep[1], ep2[1]]), skipind)
                    self.ax.plot(x, y, 'k-', alpha=self.outlinealpha)
        
    def label(self, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        for va, xst, y, inds in zip(['top', 'bottom'], [0, .5], [-self.offset, 3.**.5/2+self.offset], [[0, 2, 0], [1, 3, 1]]):
            for count, i in enumerate(inds):
                self.ax.text(xst+count*1., y, self.ellabels[i], ha='center', va=va, **kwargs)
    
    def toCart(self, quatcomps, skipinds=range(4)):#binary and ternary lines need to be plotted multiple times so returns  set of (inds,x,y)
        qc=numpy.array(quatcomps)
        #qc=qc[(qc==0.).sum(axis=1)>0]
#        x=numpy.empty(len(qc), dtype='float32')
#        y=numpy.empty(len(qc), dtype='float32')
        qindsfortern_skipind=[[1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]#sets the order of elements assuming mirror over y for skipA and skipC
        inds_x_y=[]
        for si in skipinds:
            inds=numpy.where(qc[:, si]==0.)[0]
            if len(inds)==0:
                continue
            xt, yt=self.ternaryplot.toCart(qc[inds][:, qindsfortern_skipind[si]])
            x, y=self.xy_skipind(xt, yt, si)
            inds_x_y+=[(inds, x, y)]
        return inds_x_y
    
    def scatter(self, quatcomps, c, skipinds=range(4), s='patch', **kwargs):
        if s=='patch':
            patchfcn=lambda x, y, c:self.patch_xyc(x, y, c, **kwargs)
        else:
            patchfcn=None
        inds_x_y=self.toCart(quatcomps, skipinds=skipinds)
        for inds, x, y in inds_x_y:
            if patchfcn is None:
                self.ax.scatter(x, y, c=c[inds], s=s, **kwargs)
            else:
                map(patchfcn, x, y, c[inds])
