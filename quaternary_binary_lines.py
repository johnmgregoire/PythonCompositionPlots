import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myquaternaryutility import QuaternaryPlot



class binarylines:
    def __init__(self, ax, insetax, ellabels=['A', 'B', 'C', 'D'], offset=0.02, numcomppts=21, view_azim=-159, view_elev=30, **kwargs):
        self.ax=ax
        self.insetax=insetax
        self.ellabels=ellabels
        self.stpq=QuaternaryPlot(insetax, ellabels=ellabels, offset=offset)
        comppairs=[]
        a=numpy.linspace(0, 1, 21)
        count=-1
        for i in range(4):
            for j in range(i+1, 4):
                count+=1
                b=numpy.zeros((numcomppts, 4), dtype='float64')
                b[:, i]=a
                b[:, j]=1.-a
                comppairs+=[(c1, c2) for c1, c2 in zip(b[:-1], b[1:])]
        for (c1, c2) in comppairs:
            self.stpq.line(c1, c2, fmt='-', c=self.stpq.rgb_comp([(c1+c2)/2.])[0], **kwargs)
        self.stpq.set_projection(azim=view_azim, elev=view_elev)
        self.stpq.label()

    
    def plotbinaryfom(self, comps, fom, **kwargs):
        cb=comps>.001

        ms=['<','>','^','v','s','D']
        
        count=-1
        for i in range(4):
            for j in range(i+1, 4):
                count+=1
                k, l=tuple(set(range(4))-set([i, j]))
                barr=numpy.array([numpy.logical_not(b[k]|b[l]) for b in cb]) #numpy.logical_xor(b[i], b[j])&
                if not numpy.any(barr):
                    continue
                cmps=comps[barr]
                inds=numpy.argsort(cmps[:, j])
                cmps=cmps[inds]
                cols=self.stpq.rgb_comp(cmps)
                ys=fom[barr][inds]
                for count2, (c, col, y) in enumerate(zip(cmps, cols, ys)):
                    if count2==len(ys)//2:
                        self.ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor=col, label='%s,%s' %(self.ellabels[i], self.ellabels[j]), **kwargs)
                    else:
                        self.ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor=col, **kwargs)
                        #self.ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor='None')
                for count3, (c1, col1, y1, c2, col2, y2) in enumerate(zip(cmps[:-1], cols[:-1], ys[:-1], cmps[1:], cols[1:], ys[1:])):
                    col=numpy.array([col1, col2]).mean(axis=0)
                    self.ax.plot([c1[j], c2[j]], [y1, y2], '-', c=col, **kwargs)
        
    def binarylineslegend(self, **kwargs):
        try:
            self.ax.legend(**kwargs)
        except:
            pass


