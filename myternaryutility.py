import pylab
import matplotlib.cm as cm
import numpy
from colorsys import hsv_to_rgb, rgb_to_hsv

class TernaryPlot:
    def __init__(self, ax_subplottriplet, offset=.02, minlist=[0., 0., 0.], ellabels=['A', 'B', 'C'], allowoutofboundscomps=True, outline=True):
        self.offset=offset

        self.cartendpts=numpy.float32([[0, 0], [.5, numpy.sqrt(3.)/2.], [1, 0]])
        self.ellabels=ellabels
        if not ax_subplottriplet is None:
            if isinstance(ax_subplottriplet, int):
                self.ax=pylab.subplot(ax_subplottriplet)
            elif isinstance(ax_subplottriplet, tuple):
                a, b, c=ax_subplottriplet
                self.ax=pylab.subplot(a, b, c)
            else:
                self.ax=ax_subplottriplet

            self.allowoutofboundscomps=allowoutofboundscomps
            minlist=numpy.float32(minlist)
            self.rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
            
            self.prepax(outline=outline)
            self.mappable=None
    
    def prepax(self, outline=True):
        self.ax.set_axis_off()
        self.ax.set_aspect('equal')
        self.ax.figure.hold('True')
        self.ax.set_xlim(-.10, 1.10)
        self.ax.set_ylim(-.10, 1.10)
        if outline:
            self.outline()
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

    def toCart(self, terncoordlist):
        'Given an array of triples of coords in 0-100, returns arrays of Cartesian x- and y- coords'
        terncoordlist=self.processterncoord(terncoordlist)
        aff_tcl=self.afftrans(terncoordlist)
        cartxs = 1.-aff_tcl[:, 0]-aff_tcl[:, 1]/2.
        cartys = numpy.sqrt(3) * aff_tcl[:, 1] / 2.0
        return (cartxs, cartys)

    def toComp(self, xycoordlist, process=True):
        'Given an array of triples of coords in 0-100, returns arrays of Cartesian x- and y- coords'
#        print '*', xycoordlist
        xycoordlist=numpy.float32(xycoordlist)
        if len(xycoordlist.shape)==1:
            xycoordlist=numpy.float32([xycoordlist])
        b=xycoordlist[:, 1]*2./numpy.sqrt(3.)
        a=1.-xycoordlist[:, 0]-b/2.
        c=1.-a-b
        terncoordlist=self.invafftrans(numpy.float32([a, b, c]).T)
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

        (xs, ys) = self.toCart(terncoordlist)
        self.mappable=self.ax.scatter(xs, ys, **kwargs)

#    def plot(self, terncoordlist, descriptor, **kwargs):
#        (xs, ys) = self.toCart(terncoordlist)
#        self.ax.plot(xs, ys, descriptor, **kwargs)

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

    def colorbar(self, label='', axrect=[0.86, 0.1, 0.04, 0.8], **kwargs):
        'Draws the colorbar and labels it'
        if self.mappable is None:
            print 'no mappable to create colorbar'
            return
        else:
            self.ax.figure.subplots_adjust(right=axrect[0]-.01)
            self.cbax=self.ax.figure.add_axes(axrect)
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
                self.ax.plot([ep[0], ep2[0]], [ep[1], ep2[1]], 'k-')

    def label(self, fmtstr='%.2f', takeabs=True, ternarylabels=False, hidezerocomp=False, **kwargs):#takeabs is to avoid a negative sign for ~0 negative compositions
        hal=['right', 'center', 'left']
        val=['top', 'bottom', 'top']
        xdel=[-1.*self.offset, 0, self.offset]
        ydel=[0, self.offset, 0]
        for i, ((x, y), ha, va, t, xd, yd) in enumerate(zip(self.cartendpts, hal, val, self.ellabels, xdel, ydel)):
            c=self.toComp([x, y], process=False)[0]
            if takeabs:
                c=numpy.abs(c)
            cs=None
            ternarylabels=ternarylabels or (c!=0).sum()>1
            #print c, c!=0, (c!=0).sum()>1
            if not ternarylabels:
                cs=t
            elif not self.ellabels is None:
                f=fmtstr
                cs=''.join([('%s$_{'+f+'}$') %t for t in zip(self.ellabels, c) if not (hidezerocomp and ((f %numpy.abs(t[1]))==(f %0.)))])
                #cs=(r'%s$_{'+f+r'}$%s$_{'+f+r'}$%s$_{'+f+r'}$') %tuple([t[ind] for t in zip(self.ellabels, c) for ind in range(2)])
            if not cs is None:
                self.ax.text(x+xd, y+yd, cs, ha=ha, va=va, **kwargs)

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
                ternarylabels=ternarylabels or numpy.all(c>1.e-6)
                #ternarylabels=ternarylabels or numpy.all(c!=0)
                #print c, c!=0, numpy.all(c!=0)
                if not ternarylabels:
                    cs=fmtstr %c[s]
                elif not self.ellabels is None:
                    f=fmtstr
                    cs=(r'%s$_{'+f+r'}$%s$_{'+f+r'}$%s$_{'+f+r'}$') %tuple([t[ind] for t in zip(self.ellabels, c) for ind in range(2)])
                if not cs is None:
                    self.ax.text(x+xd, y+yd, cs, ha=ha, va=va, **kwargs)

    def patch(self,coords, limits=[], **kwargs):
        '''Fill the area bounded by limits.
              Limits format: [[bmin,bmax],[lmin,lmax],[rmin,rmax]]
              Other arguments as for pylab.fill()'''
#        coords = []
#        bounds = [[1,-1,1],[1,0,-1],[-1,0,0],[1,-1,0],[1,1,-1],[-1,1,0],[0,-1,0],
#                  [0,1,-1],[-1,1,1],[0,-1,1],[0,0,-1],[-1,0,1]]
#        for pt in bounds:     #plug in values for these limits
#            for i in [0,1,2]:
#                if pt[i] == 1:
#                    pt[i] = limits[i][1]
#                else:
#                    if pt[i] == 0:pt[i] = limits[i][0]
#            for i in [0,1,2]:
#                if pt[i] == -1: pt[i] = 99 - sum(pt)
#            if self.satisfies_bounds(pt, limits): coords.append(pt)
#        coords.append(coords[0]) #close the loop
        xs, ys = self.toCart(coords)
        self.ax.fill(xs, ys, **kwargs)


    def text(self, loctriple, word, **kwargs):

        (x, y) = self.toCart([loctriple])
        self.ax.text(x[0], y[0], word, **kwargs)

    def show(self):

        self.ax.legend(loc=1)
        self.ax.set_xlim(-.10, 1.10)
        self.ax.set_ylim(-.10, 1.00)

    def rgb_comp(self, terncoordlist, affine=True):
        if affine:
            aff_tcl=self.afftrans(terncoordlist)
        else:
            aff_tcl=terncoordlist
        return aff_tcl

    def plotpoints_rgb(self, terncoordlist, affine=True, **kwargs):
        cols=self.rgb_comp(terncoordlist, affine)
        for comp, c in zip(terncoordlist, cols):
            self.scatter([comp], color=c, **kwargs)
        return cols

    def complex_to_rgb_grid(self, complex_data, invert=False):
        from numpy import angle, max, pi, sin, zeros
        phase = angle(complex_data)
        amplitude = abs(complex_data)
        amplitude = amplitude/max(max(amplitude))
        A = zeros((complex_data.shape[0], complex_data.shape[1], 3))
        A[:,:,0] = .5*(sin(phase)+1)*amplitude
        A[:,:,1] = .5*(sin(phase+pi/2)+1)*amplitude
        A[:,:,2] = .5*(-sin(phase)+1)*amplitude
        if(invert):
            return 1-A
        else:
            return A

    def rgb_compdiff(self, compdiffarr, maxcompdist=None):
        sat = ((compdiffarr**2).sum(axis=1)/2.)**.5
        huelist=[0. if cd.sum()==0. else rgb_to_hsv(*(cd/cd.sum()))[0] for cd in numpy.abs(compdiffarr)]
        if maxcompdist is None:
            sat_norm=sat/max(sat)
        else:
            sat_norm=sat/maxcompdist
            sat_norm[sat_norm>1.]=1.
        rgbarr=numpy.array([hsv_to_rgb(h, s, 1) for h, s in zip(huelist, sat_norm)])
        return rgbarr
        
    def hsdiffplot(self, terncomps, terncomps2, descriptor='o', **kwargs):
        comps=numpy.float64(terncomps)
        comps2=numpy.float64(terncomps2)
        compsdiff=comps2-comps
        rgb_arr=self.rgb_compdiff(compsdiff)
        compdist = ((compsdiff**2).sum(axis=1)/2.)**.5
        
        self.colorcompplot(comps, descriptor=descriptor, colors=rgb_arr, hollow=False, markeredgecolor='none', **kwargs)

        # color wheel axes
        self.ax.figure.subplots_adjust(left=.05, right=.7)
        self.cwax=self.ax.figure.add_axes([0.6, 0.45, 0.3, 0.45], projection='polar')
        N = 1024
        x = numpy.linspace(-1, 1, N)
        y = numpy.linspace(-1, 1, N)
        X,Y = numpy.meshgrid(x,y)
        R = numpy.sqrt(X*X + Y*Y)
        PHI = numpy.arctan2(Y, X) - numpy.pi/2
        colorgrid=self.complex_to_rgb_grid(R*numpy.exp(-1j*PHI)  * (R<1), invert=True)
        self.cwax.imshow(colorgrid, extent=[0,2*numpy.pi, 0,1024])
        self.cwax.set_rgrids([1,N/3,2*N/3], angle=45)
        self.cwax.set_xticks([numpy.pi/2, 7*numpy.pi/6, 11*numpy.pi/6])
        self.cwax.set_yticks([N/3, 2*N/3, N])
        self.cwax.set_xticklabels(['%s' % ('G'),\
                                    '%s' % ('R'),\
                                    '%s' % ('B')])
        self.cwax.set_yticklabels([\
                                    '%.3f' % (max(compdist)/3.),\
                                    '%.3f' % (2.*max(compdist)/3.),\
                                    '%.3f' % (max(compdist))])
