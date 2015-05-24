import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot

class ternaryfaces_shells:
    def __init__(self, ax, ellabels=['A', 'B', 'C', 'D'], offset=0.05, nintervals=10.):
        self.nint=1.*nintervals
        self.delta=1./self.nint
        self.ternaryplot=TernaryPlot(ax, outline=False)
        self.ax=ax
        self.offset=offset
        #self.ax.set_xlim(-.1, 2.6)
        #self.ax.set_ylim(-.1, 3.**.5/2+.1)
        self.ax.set_ylim(-.1-3.**.5/4., .1+3.**.5/4.)
        self.cartendpts=numpy.float32([[0, 0], [.5, numpy.sqrt(3.)/2.], [1, 0]])
        self.ellabels=ellabels
        self.scalefcn=lambda nshell:(self.nint-4.*nshell)/self.nint
        shift=0.
        self.shift_nshell=[]
        for nshell in range(int(self.nint//4)+1):
            self.shift_nshell+=[shift]
            shift+=self.delta*2.+2.*self.scalefcn(nshell)
        self.ax.set_xlim(-.1, shift+self.delta+1.*self.scalefcn(nshell))
        
        self.s=numpy.diff(ax.transData.transform([0., self.delta]))[0]
        self.outline()
    
    def xy_skipind(self, x, y, skipind, nshell):
        
        if skipind%2==0:
            y=-1.*y+3.**.5/2
        y-=3.**.5/2/2.
        y*=self.scalefcn(nshell)
        x+=([0.5, 1, 1.5, 0.][skipind])
        x*=self.scalefcn(nshell)
        x+=self.shift_nshell[nshell]
        return x, y
        
    def outline(self):
        for nshell in range(int(self.nint//3)):
            for skipind in range(4):#skipind=3 done in ternaryplot
                for i, ep in enumerate(self.cartendpts):
                    for ep2 in self.cartendpts[i+1:]:
                        x, y=self.xy_skipind(numpy.array([ep[0], ep2[0]]), numpy.array([ep[1], ep2[1]]), skipind, nshell)
                        self.ax.plot(x, y, 'k-')
        
    def label(self, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        for va, xst, y, inds in zip(['top', 'bottom'], [0, .5], [-3.**.5/4.-self.offset, 3.**.5/4.+self.offset], [[0, 2, 0], [1, 3, 1]]):
            for count, i in enumerate(inds):
                self.ax.text(xst+count*1., y, self.ellabels[i], ha='center', va=va, **kwargs)
    
    def toCart(self, quatcomps, skipinds=range(4), nshell=0):#binary and ternary lines need to be plotted multiple times so returns  set of (inds,x,y)
        qc=numpy.array(quatcomps)
        #qc=qc[(qc==0.).sum(axis=1)>0]
#        x=numpy.empty(len(qc), dtype='float32')
#        y=numpy.empty(len(qc), dtype='float32')
        qindsfortern_skipind=[[1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]#sets the order of elements assuming mirror over y for skipA and skipC
        inds_x_y=[]
        for si in skipinds:
            inds=numpy.where(qc[:, si]==0.)[0]
            xt, yt=self.ternaryplot.toCart(qc[inds][:, qindsfortern_skipind[si]])
            x, y=self.xy_skipind(xt, yt, si, nshell)
            inds_x_y+=[(inds, x, y)]
        return inds_x_y
    
    def scatter(self, quatcomps, c, skipinds=range(4), s=None, **kwargs):
        if s is None:
            s=self.s
        quatcomps=numpy.int32(numpy.round(quatcomps*self.nint))
        for nshell in range(int(self.nint//4)+int(self.nint%4>0)):
            ba=((quatcomps==nshell).sum(axis=1, dtype='int32')>0)&((quatcomps>=nshell).prod(axis=1, dtype='int32')>0)
            self.shellcomps=quatcomps[ba]
            shellc=c[ba]
            self.shellcomps=(self.shellcomps-nshell)/(self.nint-4.*nshell)
            inds_x_y=self.toCart(self.shellcomps, skipinds=skipinds, nshell=nshell)
            for inds, x, y in inds_x_y:
                self.ax.scatter(x, y, c=shellc[inds], **kwargs)
        if self.nint%4==0:
            ba=(quatcomps==self.nint//4).prod(axis=1, dtype='int32')>0
            if True in ba:
                self.shellcomps=quatcomps[ba]#only 1 comp but might be duplicated
                shellc=c[ba]
                for cv in shellc:
                    self.ax.scatter(self.shift_nshell[-1], 0, c=cv, **kwargs)
#            if nshell==0:
#                break
