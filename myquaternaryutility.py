import pylab 
import matplotlib.cm as cm
import numpy
from mpl_toolkits import mplot3d
#from myternaryutility import TernaryPlot
class QuaternaryPlot:
    """send a matplitlib Axis and a ternary plot is made with the utility functions. everything fractional"""
    def __init__(self, ax_subplottriplet=None, offset=.08, minlist=[0., 0., 0., 0.], ellabels=['A', 'B', 'C', 'D'], allowoutofboundscomps=True):
        self.allowoutofboundscomps=allowoutofboundscomps
        minlist=numpy.float32(minlist)
        self.rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
        
        if not ax_subplottriplet is None:
            self.offset=offset
            if isinstance(ax_subplottriplet, int):
                self.ax=pylab.subplot(ax_subplottriplet, projection='3d')
            elif isinstance(ax_subplottriplet, tuple):
                a, b, c=ax_subplottriplet
                self.ax=pylab.subplot(a, b, c, projection='3d')
            else:
                self.ax=ax_subplottriplet
            
#            if not ellabels is None:
#                for el, r in zip(ellabels, self.rangelist):
#                    print 'range of %s is %.2f to %.2f' %((el,)+tuple(r))
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
        if len(terncoordlist)==0:
            print 'no data for scatter plot'
            return
        (xs, ys, zs) = self.toCart(terncoordlist)
        self.mappable=self.ax.scatter(xs, ys, zs, **kwargs)

    def scalarmap(self, vals, norm, cmap):
        self.mappable=cm.ScalarMappable(norm=norm, cmap=cmap)
        self.mappable.set_array(vals)
        return [self.mappable.to_rgba(v) for v in vals]
    def plotbycolor(self, terncoordlist, cols,**kwargs):
        (xs, ys, zs) = self.toCart(terncoordlist)
        for xv, yv, zv, c in zip(xs, ys, zs, cols):
            self.ax.plot3D([xv], [yv], [zv], color=c, markeredgecolor=c, **kwargs)

#    def color_comp_calc(self, terncoordlist, rangelist=None):#could be made more general to allow for endpoint colors other than RGB
#        if rangelist is None:
#            rangelist=self.rangelist
#        return numpy.array([[(c-minc)/(maxc-minc) for c, (minc, maxc) in zip(tc, rangelist)] for tc in terncoordlist])
#        
#    def colorcompplot(self, terncoordlist, descriptor, colors=None, hollow=False, **kwargs):
#        (xs, ys) = self.toCart(terncoordlist)
#        if colors is None:
#            colors=self.color_comp_calc(terncoordlist)
#        for col, x, y in zip(colors, xs, ys):
#            if hollow:
#                self.ax.plot([x], [y], descriptor, markeredgecolor=col, markerfacecolor='None',  **kwargs)
#            else:
#                self.ax.plot([x], [y], descriptor, color=col, **kwargs)

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
        return (((c1-c2)**2).sum())**.5/2.**.5
    
    def cartdist_comp(self, c1, c2):#gives same answer as compdist.. I think
        x1=numpy.array(self.toCart([c1])).T[0]
        x2=numpy.array(self.toCart([c2])).T[0]
        return (((x1-x2)**2).sum())**.5

    def line(self, begin, end, fmt='k-',  **kwargs):
        (xs, ys, zs) = self.toCart([begin, end])
        self.ax.plot(xs, ys, zs, fmt, **kwargs)

    def outline(self):
        for i, ep in enumerate(self.cartendpts):
            for ep2 in self.cartendpts[i+1:]:
                self.ax.plot([ep[0], ep2[0]], [ep[1], ep2[1]], [ep[2], ep2[2]], 'k-')

    def singlelabeltext(self, c, takeabs=True, hidezerocomp=True, mult=1, fmtstr='%.3f'):
        f=fmtstr
        if takeabs:
            c=numpy.array(c)
            c=numpy.abs(c)
        cs=''.join([('%s$_{'+f+'}$') %(el, x*mult) for el, x in zip(self.ellabels, c) if not (hidezerocomp and ((f %numpy.abs(x*mult))==(f %0.)))])
        return cs
        
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
        

    def rgb_comp(self, terncoordlist, affine=True):
        cmy_cmyk=lambda a:a[:3]*(1.-a[3])+a[3]
        rgb_cmy=lambda a:1.-a
        rgb_cmyk=lambda a:rgb_cmy(cmy_cmyk(a))

        if affine:
            aff_tcl=self.afftrans(terncoordlist)
        else:
            aff_tcl=terncoordlist
        return numpy.array([rgb_cmyk(numpy.array(a)) for a in aff_tcl])
        
    def plotpoints_rgb(self, terncoordlist, affine=True, **kwargs):
        cols=self.rgb_comp(terncoordlist, affine)
        for comp, c in zip(terncoordlist, cols):
            self.scatter([comp], color=c, **kwargs)
        
    
    def filterbydistancefromline(self, terncoordlist, compend1, compend2, critdist, betweenpoints=True,  affine=False, invlogic=False, returnall=False): #not sure fi affine transformation makes sense here but haven't through through it
        #see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        xyzarr=numpy.array(self.toCart(terncoordlist, affine=affine)).T
        xyz1=numpy.array(self.toCart([compend1], affine=affine)).T[0]
        xyz2=numpy.array(self.toCart([compend2], affine=affine)).T[0]
        lineparameter=None
        distfromlin=numpy.array([numpy.linalg.norm(numpy.cross(xyz2-xyz1, xyz1-xyz))/numpy.linalg.norm((xyz2-xyz1)) for xyz in xyzarr])
        if betweenpoints:
            lineparameter=numpy.array([-numpy.inner(xyz2-xyz1, xyz1-xyz)/numpy.linalg.norm((xyz2-xyz1))**2 for xyz in xyzarr])
            if invlogic:
                inds=numpy.where(numpy.logical_not((distfromlin<=critdist) & (lineparameter>=0) & (lineparameter<=1)))[0]
            else:
                inds=numpy.where((distfromlin<=critdist) & (lineparameter>=0) & (lineparameter<=1))[0]
        else:
            if invlogic:
                inds=numpy.where(numpy.logical_not(distfromlin<=critdist))[0]
            else:
                inds=numpy.where(distfromlin<=critdist)[0]
        if returnall:
            if lineparameter is None:
                lineparameter=numpy.array([-numpy.inner(xyz2-xyz1, xyz1-xyz)/numpy.linalg.norm((xyz2-xyz1))**2 for xyz in xyzarr])
            return inds, distfromlin, lineparameter
        else:
            return inds
        
    def filterbydistancefromplane(self, terncoordlist, compvert0, compvert1, compvert2, critdist, withintriangle=True,  affine=False, invlogic=False, returnall=False): #not sure fi affine transformation makes sense here but haven't through through it
        xyzarr=numpy.array(self.toCart(terncoordlist, affine=affine)).T
        xyz0=numpy.array(self.toCart([compvert0], affine=affine)).T[0]
        xyz1=numpy.array(self.toCart([compvert1], affine=affine)).T[0]
        xyz2=numpy.array(self.toCart([compvert2], affine=affine)).T[0]
        nhat=numpy.cross(xyz1-xyz0, xyz2-xyz1)
        nhat/=numpy.linalg.norm(nhat)
        
        distfromplane=numpy.array([numpy.abs(numpy.linalg.norm(numpy.inner(nhat, xyz1-xyz))) for xyz in xyzarr])
        intriangle=None
        if returnall or withintriangle:
            xphat=xyz1-xyz0
            xphat/=numpy.linalg.norm(xphat)
            yphat=numpy.cross(nhat, xphat)
            xyparr=numpy.array([[numpy.inner(xyz-xyz0, xphat), numpy.inner(xyz-xyz0, yphat)] for xyz in xyzarr]) #this makes xyz0 the origin, x axis points to xyz1
            xyp_verts=numpy.array([[numpy.inner(xyz-xyz0, xphat), numpy.inner(xyz-xyz0, yphat)] for xyz in [xyz0, xyz1, xyz2]])
        if withintriangle:
            intriangle=numpy.array([self.point_wrt_polygon(xyp, xyp_verts) for xyp in xyparr])
            if invlogic:
                inds=numpy.where(numpy.logical_not((distfromplane<=critdist) & (intriangle==1)))[0]
            else:
                inds=numpy.where((distfromplane<=critdist) & (intriangle==1))[0]
        else:
            if invlogic:
                inds=numpy.where(numpy.logical_not(distfromplane<=critdist))[0]
            else:
                inds=numpy.where(distfromplane<=critdist)[0]
        if returnall:
            if intriangle is None:
                intriangle=numpy.array([self.point_wrt_polygon(xyp, xyp_verts) for xyp in xyparr])
            return inds, distfromplane, xyparr, xyp_verts,intriangle#xyparr is array of x,y projections into the selected plan with xyz0 as the origin
        else:
            return inds
            
    def plotfomalonglineparameter(self, ax, lineparameter, fom, compend1=None, compend2=None, lineparticks=numpy.linspace(0, 1, 4), labelfmtstr='%.3f', ticklabelkwargdict={}, **kwargs):
        sortinds=numpy.argsort(lineparameter)
        ax.plot(lineparameter[sortinds], fom[sortinds], **kwargs)
        if not lineparticks is None:
            tl=[]
            for i in lineparticks:
                c=compend1+(compend2-compend1)*i
                tl+=[self.singlelabeltext(c, fmtstr=labelfmtstr)]
            ax.xaxis.set_ticks(lineparticks)
            ax.xaxis.set_ticklabels(tl, **ticklabelkwargdict)
    
    def plotfominselectedplane(self, ax, xyparr, fom, xyp_verts=None, vertcomps_labels=None, vertlw=1., labelfmtstr='%.3f', **kwargs):
        ax.scatter(xyparr[:, 0], xyparr[:, 1], c=fom, **kwargs)
        xyp_verts=list(xyp_verts)
        if not vertlw is None:
            for xyp0, xyp1 in zip(xyp_verts, xyp_verts[1:]+[xyp_verts[0]]):
                ax.plot([xyp0[0], xyp1[0]], [xyp0[1], xyp1[1]], 'k-', lw=vertlw)
        if not vertcomps_labels is None:
            for xyp, c, ha, va in zip(xyp_verts, vertcomps_labels, ['right', 'left', 'center'], ['top', 'center', 'bottom']):
                if c is None:
                    continue
                ax.text(xyp[0], xyp[1], self.singlelabeltext(c, fmtstr=labelfmtstr), ha=ha, va=va, fontsize=16)
        ax.set_axis_off()
        ax.set_aspect('equal')

    def point_wrt_polygon(self, xy, xyarr_vert, inside=True,  perimeter=True, outside=False, tol=1.e-10):
        x=xy[0]
        y=xy[1]
        
        n = len(xyarr_vert)
        insidetest =False

        p1x,p1y = xyarr_vert[0]
        for i in range(n+1):
            p2x,p2y = xyarr_vert[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            insidetest = not insidetest
            p1x,p1y = p2x,p2y
        ans=(inside and insidetest) or (outside and not insidetest)
        if perimeter:
            xyarr_vert_cyc=numpy.concatenate([xyarr_vert, [xyarr_vert[0]]])
            pertest=numpy.any([numpy.cross(v2-v1, xy-v1)**2.<=(tol*(numpy.linalg.norm(v2-v1)*numpy.linalg.norm(xy-v1)))**2. for v1, v2 in zip(xyarr_vert_cyc[1:], xyarr_vert_cyc[:-1])])
            ans=ans or pertest
        return ans
