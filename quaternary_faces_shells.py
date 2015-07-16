import matplotlib.cm as cm
import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os
from matplotlib.patches import CirclePolygon
#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot

class ternaryfaces_shells:
    def __init__(self, ax, ellabels=['A', 'B', 'C', 'D'], offset=None, nintervals=10., outlinealpha=0.2, patchscale=1.):
        self.outlinealpha=outlinealpha
        self.nint=1.*nintervals
        self.delta=1./self.nint
        self.ternaryplot=TernaryPlot(ax, outline=False)
        self.ax=ax
        self.ax.set_aspect(1.)
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
        
        self.patch_xyc=lambda x, y, c, **kwargs:self.ax.add_patch(CirclePolygon((x, y),radius=patchscale*self.delta/3.**.5,resolution=6, color=c, **kwargs))
        if outlinealpha>0:
            self.outline()
        if offset is None:
            self.offset=self.delta
    
    def xy_skipind(self, x, y, skipind, nshell):
        
        if skipind%2==0:
            y=-1.*y+3.**.5/2
        y-=3.**.5/2/2.
        y*=self.scalefcn(nshell)
        x+=([0.5, 1, 1.5, 0.][skipind])
        x*=self.scalefcn(nshell)
        x+=self.shift_nshell[nshell]
        return x, y
        
    def outline(self, **kwargs):
        for nshell in range(int(self.nint//4)+int(self.nint%4>0)):
            for skipind in range(4):
                skipfirstline=skipind!=3
                for i, ep in enumerate(self.cartendpts):
                    for ep2 in self.cartendpts[i+1:]:
                        if skipfirstline:
                            skipfirstline=False
                            continue
                        x, y=self.xy_skipind(numpy.array([ep[0], ep2[0]]), numpy.array([ep[1], ep2[1]]), skipind, nshell)
                        self.ax.plot(x, y, 'k-', alpha=self.outlinealpha, **kwargs)
        
    def label(self, onlyterns=False, allelements=False, primeelements=False, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        for va, xst, y, inds in zip(['top', 'bottom'], [0, .5], [-3.**.5/4.-self.offset, 3.**.5/4.+self.offset], [[0, 2, 0], [1, 3, 1]]):
            for count, i in enumerate(inds):
                self.ax.text(xst+count*1., y, self.ellabels[i], ha='center', va=va, **kwargs)
        if onlyterns:
            return
        for nshell in range(1, int(self.nint//4)+int(self.nint%4>0)):
            for va, xst, y, inds in zip(['top', 'bottom'], [0, .5], [-3.**.5/4.*self.scalefcn(nshell)-self.offset, 3.**.5/4.*self.scalefcn(nshell)+self.offset], [[2, 0], [3, 1]]):
                for count, i in enumerate(inds):
                    l=self.ellabels[i]
                    if primeelements:
                        l+="$'$"*nshell
                    else:
                        l+=(r'$_{%d}$' %int(round(100*(1.-3*nshell*self.delta))))
                    x=(xst+(count+1)*1.)*self.scalefcn(nshell)+self.shift_nshell[nshell]
                    if allelements:
                        temp=copy.copy(self.ellabels)
                        temp.pop(i)
                        l+=''.join([el+(r'$_{%d}$' %int(round(100*(nshell*self.delta)))) for el in temp])
                    self.ax.text(x, y, l, ha='center', va=va, **kwargs)
    
    def toCart(self, quatcomps, skipinds=range(4), nshell=0):#binary and ternary lines need to be plotted multiple times so returns  set of (inds,x,y)
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
            x, y=self.xy_skipind(xt, yt, si, nshell)
            inds_x_y+=[(inds, x, y)]
        return inds_x_y
    
    def scatter(self, quatcomps, c, skipinds=range(4), s='patch', **kwargs):
        if s=='patch':
            patchfcn=lambda x, y, c:self.patch_xyc(x, y, c, **kwargs)
        else:
            patchfcn=None
        quatcomps=numpy.int32(numpy.round(quatcomps*self.nint))
        for nshell in range(int(self.nint//4)+int(self.nint%4>0)):
            ba=((quatcomps==nshell).sum(axis=1, dtype='int32')>0)&((quatcomps>=nshell).prod(axis=1, dtype='int32')>0)
            self.shellcomps=quatcomps[ba]
            shellc=c[ba]
            self.shellcomps=(self.shellcomps-nshell)/(self.nint-4.*nshell)
            inds_x_y=self.toCart(self.shellcomps, skipinds=skipinds, nshell=nshell)
            for inds, x, y in inds_x_y:
                if patchfcn is None:
                    self.ax.scatter(x, y, c=shellc[inds], s=s, **kwargs)
                else:
                    map(patchfcn, x, y, shellc[inds])
        if self.nint%4==0: #single point with no frame
            ba=(quatcomps==self.nint//4).prod(axis=1, dtype='int32')>0
            if True in ba:
                self.shellcomps=quatcomps[ba]#only 1 comp but might be duplicated
                shellc=c[ba]
                
                if patchfcn is None:
                    for cv in shellc:
                        self.ax.scatter(self.shift_nshell[-1], 0, c=cv, s=s, **kwargs)
                else:
                    [patchfcn(self.shift_nshell[-1], 0, cv) for cv in shellc]
                    
    def quatscatter(self, quatcomps, c, skipinds=range(4), azim=-60, elev=30, alphaall=.2, alphashell=1., fontsize=14, outline=True,  **kwargs):
        numsubs=int(self.nint//4)+1
        quatcomps=numpy.int32(numpy.round(quatcomps*self.nint))
        for nshell in range(int(self.nint//4)+int(self.nint%4>0)):
            ba=((quatcomps==nshell).sum(axis=1, dtype='int32')>0)&((quatcomps>=nshell).prod(axis=1, dtype='int32')>0)
            shellcomps=quatcomps[ba]
            shellc=c[ba]
            
            q=QuaternaryPlot((1, numsubs, nshell+1), outline=outline)
            if alphaall>0:
                q.scatter(quatcomps*1./self.nint,c=c, alpha=alphaall, **kwargs)
            if alphashell>0:
                q.scatter(shellcomps*1./self.nint,c=shellc, alpha=alphashell, **kwargs)
            if fontsize>0:
                q.label(ha='center', va='center', fontsize=fontsize)
            q.set_projection(azim=azim, elev=elev)

        if self.nint%4==0: #single point with no frame
            ba=(quatcomps==self.nint//4).prod(axis=1, dtype='int32')>0
            if True in ba:
                shellcomps=quatcomps[ba]#only 1 comp but might be duplicated
                shellc=c[ba]
                q=QuaternaryPlot((1, numsubs, numsubs), outline=outline)
                q.scatter(quatcomps*1./self.nint,c=c, alpha=alphaall, **kwargs)
                q.scatter(shellcomps*1./self.nint,c=shellc, alpha=alphashell, **kwargs)
                if fontsize>0:
                    q.label(ha='center', va='center', fontsize=fontsize)
                q.set_projection(azim=azim, elev=elev)
    def quatplot3D(self, quatcomps, c, skipinds=range(4), azim=-60, elev=30, alphaall=.2, alphashell=1., fontsize=14, outline=True,  **kwargs):
        numsubs=int(self.nint//4)+1
        quatcomps=numpy.int32(numpy.round(quatcomps*self.nint))
        for nshell in range(int(self.nint//4)+int(self.nint%4>0)):
            ba=((quatcomps==nshell).sum(axis=1, dtype='int32')>0)&((quatcomps>=nshell).prod(axis=1, dtype='int32')>0)
            shellcomps=quatcomps[ba]
            shellc=c[ba]
            
            q=QuaternaryPlot((1, numsubs, nshell+1), outline=outline)
            if alphaall>0:
                q.plot3D(quatcomps*1./self.nint,c, alpha=alphaall, **kwargs)
            if alphashell>0:
                q.plot3D(shellcomps*1./self.nint,shellc, alpha=alphashell, **kwargs)
            if fontsize>0:
                q.label(ha='center', va='center', fontsize=fontsize)
            q.set_projection(azim=azim, elev=elev)

        if self.nint%4==0: #single point with no frame
            ba=(quatcomps==self.nint//4).prod(axis=1, dtype='int32')>0
            if True in ba:
                shellcomps=quatcomps[ba]#only 1 comp but might be duplicated
                shellc=c[ba]
                q=QuaternaryPlot((1, numsubs, numsubs), outline=outline)
                q.plot3D(quatcomps*1./self.nint,c, alpha=alphaall, **kwargs)
                q.plot3D(shellcomps*1./self.nint,shellc, alpha=alphashell, **kwargs)
                if fontsize>0:
                    q.label(ha='center', va='center', fontsize=fontsize)
                q.set_projection(azim=azim, elev=elev)
