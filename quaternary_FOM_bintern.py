import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot

def make4ternaxes(ellabels=['A', 'B', 'C', 'D']):
    
    fig=pylab.figure(figsize=(12, 5))
    
    w=.21
    h=.9
    axl=[]
    for i in range(4):
        axl+=[fig.add_axes([.02+i*(w+.01), (1.-h)/2., w, h])]

    stpl=[]
    for i, ax in zip([3, 2, 1, 0], axl):
        eltemp=copy.copy(ellabels)
        eltemp.pop(i)
        stp=TernaryPlot(ax, ellabels=eltemp)
        stp.label()
        stpl+=[stp]
    return axl, stpl



def scatter_4axes(comps, fom, stpl, **kwargs):
    for i, stp in zip([3, 2, 1, 0], stpl):
        d=comps[:, i]
        l=[0, 1, 2, 3]
        l.pop(i)
        ci=numpy.array(l)
        inds=numpy.where(d==0.)[0]
        if len(inds)>0:
            stp.scatter(comps[inds][:, ci], c=fom[inds], **kwargs)

def plotbinarylines_axandinset(ellabels=['A', 'B', 'C', 'D'], mainax=[.3, .12, .6, .83], insetax=[0, .7, .2, .3], numcomppts=21, view_azim=-159, view_elev=30, **kwargs):
    fig=pylab.figure(figsize=(8, 5))
    ax=fig.add_axes(mainax)
    ax2=fig.add_axes(insetax, projection='3d')
    stpq=QuaternaryPlot(ax2, ellabels=ellabels)
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
        stpq.line(c1, c2, fmt='-', c=stpq.rgb_comp([(c1+c2)/2.])[0], **kwargs)
    stpq.set_projection(azim=view_azim, elev=view_elev)
    stpq.label()
    return ax, ax2
    
def plotbinarylines_quat(ax, comps, fom, ellabels=['A', 'B', 'C', 'D'], legloc=7, **kwargs):
    cb=comps>.001
    qtemp=QuaternaryPlot(None)
    ms=['<','>','^','v','s','D']
    
    #    for i in range(4):
    #        barr=cb[:, i]>.999
    #        if not numpt.any(barr):
    #            continue
    #        ax.plot(c[j], y, ms[count], c=c, ms=ms, markeredgecolor='None', label='%s,%s' %(ellabels[i], ellabels[j]), **kwargs)
        
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
            cols=qtemp.rgb_comp(cmps)
            ys=fom[barr][inds]
            for count2, (c, col, y) in enumerate(zip(cmps, cols, ys)):
                if count2==len(ys)//2:
                    ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor=col, label='%s,%s' %(ellabels[i], ellabels[j]), **kwargs)
                else:
                    ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor=col, **kwargs)
                    #ax.plot(c[j], y, marker=ms[count], c=col, markeredgecolor='None')
            for count3, (c1, col1, y1, c2, col2, y2) in enumerate(zip(cmps[:-1], cols[:-1], ys[:-1], cmps[1:], cols[1:], ys[1:])):
                col=numpy.array([col1, col2]).mean(axis=0)
                ax.plot([c1[j], c2[j]], [y1, y2], '-', c=col, **kwargs)
    try:
        ax.legend(loc=legloc)
    except:
        pass
