import matplotlib.cm as cm
import numpy
import pylab
import h5py, operator, copy, os


#pylab.rc('font',**{'family':'serif''serif':['Times New Roman']})
#pylab.rcParams['font.family']='serif'
#pylab.rcParams['font.serif']='Times New Roman'
pylab.rc('font', family='serif', serif='Times New Roman')
#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot


def make30ternaxes(ellabels=['A', 'B', 'C', 'D'], fig=None):
    if fig is None:
        fig=pylab.figure(figsize=(14, 8))

    axl=[]

    delw=.02
    delh=.04
    yscalefct=[1.]*6#[.7, .7, .7, .7, .7, .7]
    xscalefct=[1., .9, .67, .67, .67, .67]
    npercol=numpy.array([3,4,5,6,6,6])
    colw=(1./npercol)
    colw[0]*=.85
    colw[1]*=.9
    colw=(colw/colw.sum())*(.9-len(npercol)*delw)
    #colw/=colw.sum()

    plotcount=0
    cumwidth=0
    for nc, cw, xsf, ysf in zip(npercol, colw, xscalefct, yscalefct):
        w=xsf*cw
        h=ysf*((1.-delh*(nc+1))/nc)
        cumwidth+=cw+delw
        for ic in range(nc):
            axl+=[fig.add_axes([cumwidth-cw/2.-w/2., 1.-(delh+ic*(h+delh))-h, w, h])]


    stpl=[]
    xpos=[.27]*30
    #    xpos[0:3]=[.35, .29, .28]
    #    xpos[-1]=.26
    for count, (ax, xp) in enumerate(zip(axl, xpos)):
        stp=TernaryPlot(ax, ellabels=ellabels[:3], offset=.03)
        stp.label(fontsize=15)#,fontdict={'fontname':'Times New Roman'})
        stpl+=[stp]
        
        
        if count<29:
            stp.ax.text(xp, .8, '%s$_{%.2f}$' %(ellabels[3], (count*.0333)), ha='right', va='center', fontsize=17)
        else:
            stp.ax.text(xp, .8, '%s$_{%.2f-%d}$' %(ellabels[3], (count*.0333), 1), ha='right', va='center', fontsize=17)
    return axl, stpl



def scatter_30axes(comps, fom, stpl, s=18, cb=False, cbrect=(.91, .3, .03, .4), cblabel='', **kwargs):
    abc=comps[:, :3]
    abc[abc.sum(axis=1)==0.]=numpy.array([1., 1., 1.])/3.
    abc=numpy.array([c/c.sum() for c in abc])
    d=comps[:, 3]
    d30=numpy.round(d*30.)
    d30[d30==30.]=29.
    dlims=numpy.array([0., 1., 2., 3.])
    marks=[('o', 1., 1.), ('D', .9, .7),('s', .8, .5)]
    sl=s*numpy.array([2.3, 1., .58, .38, .38, .42, .38, .46, .8, 1., 1.5]+[1.5]*20)
    for i, (stp, sv) in enumerate(zip(stpl, sl)):
        inds=numpy.where((d30==i))[0]
        #print a, b, len(inds)
        if len(inds)>0:
            stp.scatter(abc[inds], c=fom[inds], marker='o', s=20, **kwargs)
    if cb:
        cbax=stp.ax.figure.add_axes(cbrect)
        if 'extend' in kwargs.keys():
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'], extend=kwargs['extend'])
        else:
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'])
        sm.set_array(fom)
        cb=stp.ax.figure.colorbar(sm, cax=cbax)
        cb.set_label(cblabel, fontsize=18)
#make30ternaxes()
#pylab.show()
