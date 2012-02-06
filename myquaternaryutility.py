import pylab 
import matplotlib.cm as cm
import numpy
from mpl_toolkits import mplot3d
#from myternaryutility import TernaryPlot
class QuaternaryPlot:
    """send a matplitlib Axis and a ternary plot is made with the utility functions. everything fractional"""
    def __init__(self, ax_subplottriplet, offset=.08, minlist=[0., 0., 0., 0.], ellabels=['A', 'B', 'C', 'D'], allowoutofboundscomps=True):
        self.offset=offset
        if isinstance(ax_subplottriplet, tuple) or isinstance(ax_subplottriplet, int):
            self.ax=pylab.subplot(ax_subplottriplet, projection='3d')
        else:
            self.ax=ax_subplottriplet
        self.allowoutofboundscomps=allowoutofboundscomps
        minlist=numpy.float32(minlist)
        self.rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
        if not ellabels is None:
            for el, r in zip(ellabels, self.rangelist):
                print 'range of %s is %.2f to %.2f' %((el,)+tuple(r))
        self.ax.set_axis_off()
        self.ax.set_aspect('equal')
        self.ax.figure.hold('True')
        self.ax.set_xlim(-.10, 1.10)
        self.ax.set_ylim(-.10, 1.10)
        self.cartendpts=numpy.float32([[0, 0, 0], [.5, numpy.sqrt(3.)/2., 0], [1, 0, 0], [.5, .5/numpy.sqrt(3.), numpy.sqrt(2./3.)]])
        self.ellabels=ellabels
        self.outline()
        self.mappable=None
    
    def set_projection(self, azim=None, elev=None):
        self.ax.view_init(azim=azim, elev=elev)
    
    def processterncoord(self, terncoordlist, removepoints=True):
        terncoordlist=numpy.float32(terncoordlist)
        if len(terncoordlist.shape)==1:
            terncoordlist=numpy.float32([terncoordlist])
        if removepoints and not self.allowoutofboundscomps:
            terncoordlist=numpy.float32([t for t in terncoordlist if (not removepoints) or numpy.all(t>=self.rangelist[:, 0]) and numpy.all(t<=self.rangelist[:, 1])])
        return terncoordlist
        
            
    def afftrans(self, terncoordlist):
        terncoordlist=self.processterncoord(terncoordlist)
        diff=self.rangelist[:, 1]-self.rangelist[:, 0]
        mn=self.rangelist[:, 0]
        return numpy.float32([(tc-mn)/diff for tc in terncoordlist])

    def invafftrans(self, terncoordlist):
        terncoordlist=self.processterncoord(terncoordlist, removepoints=False)
        diff=self.rangelist[:, 1]-self.rangelist[:, 0]
        mn=self.rangelist[:, 0]
        return numpy.float32([c*diff+mn for c in terncoordlist])
        
    def toCart(self, terncoordlist, affine=True):
        'Given an array of triples of coords in 0-100, returns arrays of Cartesian x- and y- coords'
        terncoordlist=self.processterncoord(terncoordlist)
        if affine:
            aff_tcl=self.afftrans(terncoordlist)
        else:
            aff_tcl=terncoordlist
        a=aff_tcl[:, 0]
        b=aff_tcl[:, 1]
        c=aff_tcl[:, 2]
        d=aff_tcl[:, 3]
        x=1.-a-b/2.-d/2.
        y=b/2.*(3.**.5)+d/2./(3.**.5)
        z=d*(2.**.5)/(3.**.5)
        return (x, y, z)

    def toComp(self, xycoordlist, process=True, affine=True):
        'Given an array of triples of coords in 0-100, returns arrays of Cartesian x- and y- coords'
#        print '*', xycoordlist
        xycoordlist=numpy.float32(xycoordlist)
        if len(xycoordlist.shape)==1:
            xycoordlist=numpy.float32([xycoordlist])
        x=xycoordlist[:, 0]
        y=xycoordlist[:, 1]
        z=xycoordlist[:, 2]
        d=numpy.sqrt(3./2.)*z
        b=y*2./numpy.sqrt(3.)-z/numpy.sqrt(6.)
        a=1.-x-b/2.-d/2.
        c=1.-a-b-d
        if affine:
            terncoordlist=self.invafftrans(numpy.float32([a, b, c, d]).T)
        else:
            terncoordlist=numpy.float32([a, b, c, d]).T
#        print 'a', a
#        print 'b', b
#        print 'c', c
#        print numpy.float32([a, b, c]).T
#        print 'tcl', terncoordlist
        if process:
            terncoordlist=self.processterncoord(terncoordlist)
#        print 'ptcl', terncoordlist
        return terncoordlist
        
    def scatter(self, terncoordlist, **kwargs):
        'Scatterplots data given in triples, with the matplotlib keyword arguments'
        
        (xs, ys, zs) = self.toCart(terncoordlist)
        self.mappable=self.ax.scatter(xs, ys, zs, **kwargs)

    def plot(self, terncoordlist, descriptor, **kwargs):
        (xs, ys) = self.toCart(terncoordlist)
        self.ax.plot(xs, ys, descriptor, **kwargs)

    def color_comp_calc(self, terncoordlist, rangelist=None):#could be made more general to allow for endpoint colors other than RGB
        if rangelist is None:
            rangelist=self.rangelist
        return numpy.array([[(c-minc)/(maxc-minc) for c, (minc, maxc) in zip(tc, rangelist)] for tc in terncoordlist])
        
    def colorcompplot(self, terncoordlist, descriptor, colors=None, hollow=False, **kwargs):
        (xs, ys) = self.toCart(terncoordlist)
        if colors is None:
            colors=self.color_comp_calc(terncoordlist)
        for col, x, y in zip(colors, xs, ys):
            if hollow:
                self.ax.plot([x], [y], descriptor, markeredgecolor=col, markerfacecolor='None',  **kwargs)
            else:
                self.ax.plot([x], [y], descriptor, color=col, **kwargs)

    def colorbar(self, label='', **kwargs):
        'Draws the colorbar and labels it'
        if self.mappable is None:
            print 'no mappable to create colorbar'
            return
        else:            
            #cb=pylab.colorbar()
            self.ax.figure.subplots_adjust(right=0.85)
            self.cbax=self.ax.figure.add_axes([0.86, 0.1, 0.04, 0.8])
            f=self.ax.figure.colorbar
            try:
                cb=self.ax.figure.colorbar(self.mappable, cax=self.cbax, **kwargs)
            except:
                cb=self.ax.figure.colorbar(self.mappable, cax=self.cbax)
            try:
                cb.set_label(label, **kwargs)
            except:
                cb.set_label(label)
            return cb
    
    def compdist(self, c1, c2):
        return ((c1-c2)**2).sum()/2.**.5
    
    def compdist_cart(self, c1, c2):
        return self.compdist(self.toCart([c1])[0], self.toCart([c2])[0])
            
    def line(self, begin, end, fmt='k-',  **kwargs):
        (xs, ys) = self.toCart([begin, end])
        self.ax.plot(xs, ys, fmt, **kwargs)

    def outline(self):
        for i, ep in enumerate(self.cartendpts):
            for ep2 in self.cartendpts[i+1:]:
                self.ax.plot([ep[0], ep2[0]], [ep[1], ep2[1]], [ep[2], ep2[2]], 'k-')

    def label(self, fmtstr='%.2f', takeabs=True, ternarylabels=False, hidezerocomp=False, rndzero=.0001, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        temp1=numpy.ones(8, dtype='float32')/3.
        temp1[3]=0.
        temp2=numpy.zeros(8, dtype='float32')
        temp2[3]=1.
        labcoords=[(1.+self.offset)*numpy.array(self.toCart(temp2[3-i:3-i+4], affine=False)).T[0]-self.offset*numpy.array(self.toCart(temp1[3-i:3-i+4], affine=False)).T[0] for i in range(4)]
        for i, ((x, y, z), endpt, t) in enumerate(zip(labcoords, self.cartendpts, self.ellabels)):
            c=self.toComp(endpt, process=False)[0]
            if takeabs:
                c=numpy.abs(c)
            cs=None
            ternarylabelsx=ternarylabels or (numpy.abs(c)>rndzero).sum()>1
            #print c, c!=0, (c!=0).sum()>1
            if not ternarylabelsx:
                cs=t
            elif not self.ellabels is None:
                f=fmtstr
                cs=''.join([('%s$_{'+f+'}$') %t for t in zip(self.ellabels, c) if not (hidezerocomp and ((f %numpy.abs(t[1]))==(f %0.)))])
                #cs=(r'%s$_{'+f+r'}$%s$_{'+f+r'}$%s$_{'+f+r'}$') %tuple([t[ind] for t in zip(self.ellabels, c) for ind in range(2)])
            if not cs is None:
                self.ax.text(x, y, z, cs, zdir=None, **kwargs)
                
#    def label(self, fmtstr='%.2f', takeabs=True, ternarylabels=False, hidezerocomp=False, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
#        hal=['right', 'center', 'left']
#        val=['top', 'bottom', 'top']
#        xdel=[-1.*self.offset, 0, self.offset]
#        ydel=[0, self.offset, 0]
#        for i, ((x, y), ha, va, t, xd, yd) in enumerate(zip(self.cartendpts, hal, val, self.ellabels, xdel, ydel)):
#            c=self.toComp([x, y], process=False)[0]
#            if takeabs:
#                c=numpy.abs(c)
#            cs=None
#            ternarylabels=ternarylabels or (c!=0).sum()>1
#            #print c, c!=0, (c!=0).sum()>1
#            if not ternarylabels:
#                cs=t
#            elif not self.ellabels is None:
#                f=fmtstr
#                cs=''.join([('%s$_{'+f+'}$') %t for t in zip(self.ellabels, c) if not (hidezerocomp and ((f %numpy.abs(t[1]))==(f %0.)))])
#                #cs=(r'%s$_{'+f+r'}$%s$_{'+f+r'}$%s$_{'+f+r'}$') %tuple([t[ind] for t in zip(self.ellabels, c) for ind in range(2)])
#            if not cs is None:
#                self.ax.text(x+xd, y+yd, cs, ha=ha, va=va, **kwargs)
        
    def grid(self, nintervals=4, fmtstr='%0.2f', takeabs=True, ternarylabels=False, printticklabels=True, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        lstyle = {'color': '0.6',
         #'dashes': (1, 1),
         'linewidth': 1.}
        rot=[60, 0, 300]
        hal=['right', 'left', 'center']
        val=['center', 'center', 'top']
        xdel=[-1.*self.offset, self.offset, 0]
        ydel=[0, 0, -1.*self.offset]
        side=[1, 2, 0]
        if isinstance(printticklabels, bool):
            if printticklabels:
                printticklabels=[True]*(nintervals-1)
            else:
                printticklabels=[False]*(nintervals-1)
        elif isinstance(printticklabels, list) and not isinstance(printticklabels, bool):
            printticklabels=[i in printticklabels for i in range(nintervals-1)]
        n=nintervals
        ep=self.cartendpts
        for i, j, k, r, ha, va, xd, yd, s in zip([0, 1, 2], [1, 2, 0], [2, 0, 1], rot, hal, val, xdel, ydel, side):
            for m, b in zip(range(1, n), printticklabels):
                x, y=((n-m)*ep[i]+m*ep[j])/n
                xe, ye=((n-m)*ep[k]+m*ep[j])/n
                self.ax.plot([x, xe], [y, ye], **lstyle)
                if not b:
                    continue
                c=self.toComp([x, y], process=False)[0]
                if takeabs:
                    c=numpy.abs(c)                    
                cs=None
                ternarylabels=ternarylabels or numpy.all(c!=0)
                #print c, c!=0, numpy.all(c!=0)
                if not ternarylabels:
                    cs=fmtstr %c[s]
                elif not self.ellabels is None:
                    f=fmtstr
                    cs=(r'%s$_{'+f+r'}$%s$_{'+f+r'}$%s$_{'+f+r'}$') %tuple([t[ind] for t in zip(self.ellabels, c) for ind in range(2)])
                if not cs is None:
                    self.ax.text(x+xd, y+yd, cs, ha=ha, va=va, **kwargs)
    
    def patch(self,limits, **kwargs): 
        '''Fill the area bounded by limits.
              Limits format: [[bmin,bmax],[lmin,lmax],[rmin,rmax]]
              Other arguments as for pylab.fill()'''
        coords = []
        bounds = [[1,-1,1],[1,0,-1],[-1,0,0],[1,-1,0],[1,1,-1],[-1,1,0],[0,-1,0],
                  [0,1,-1],[-1,1,1],[0,-1,1],[0,0,-1],[-1,0,1]]
        for pt in bounds:     #plug in values for these limits
            for i in [0,1,2]:
                if pt[i] == 1: 
                    pt[i] = limits[i][1]
                else:
                    if pt[i] == 0:pt[i] = limits[i][0]
            for i in [0,1,2]:
                if pt[i] == -1: pt[i] = 99 - sum(pt) 
            if self.satisfies_bounds(pt, limits): coords.append(pt) 
        coords.append(coords[0]) #close the loop
        xs, ys = self.toCart(coords)
        self.ax.fill(xs, ys, **kwargs) 


    def text(self, loctriple, word, **kwargs):
        
        (x, y) = self.toCart([loctriple])
        self.ax.text(x[0], y[0], word, **kwargs)

    def show(self):
        
        self.ax.legend(loc=1)
        self.ax.set_xlim(-.10, 1.10)
        self.ax.set_ylim(-.10, 1.00)

    def plotabcprojection(self, terncoordlist, **kwargs):
        terncoordlistproj=[numpy.append(c[:3]/(c[:3].sum()), self.rangelist[3][0]) for c in terncoordlist]
        (xs, ys, zs) = self.toCart(terncoordlistproj)
        i=numpy.argmin(xs+ys)
        print terncoordlist[i], terncoordlist[i].sum()
        print terncoordlistproj[i]
        print xs[i], ys[i], zs[i]
        self.mappable=self.ax.scatter(xs, ys, zs, **kwargs)
        

