import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os


#pylab.rc('font',**{'family':'serif''serif':['Times New Roman']})
#pylab.rcParams['font.family']='serif'
#pylab.rcParams['font.serif']='Times New Roman'
pylab.rc('font', family='serif', serif='Times New Roman')
#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot


def make10ternaxes(ellabels=['A', 'B', 'C', 'D'], fig=None, fontsize=17):
    if fig is None:
        fig=pylab.figure(figsize=(12, 8))
        
    ax_xc=[]
    ax_yc=[]
    xcdel=[.18, .19, .065, .1, .04, .05, .055, .03, .02, .02]
    ax_yc=[.49, .6, .38, .64, .48, .34, .66, .53, .42, .32]
    for i in range(10):
        if i==0:
            ax_xc+=[xcdel[i]]
        else:
            ax_xc+=[ax_xc[-1]+xcdel[i]]
        #ax_yc+=[.5+((i%2)*2.-1.)*((i>0)*.1+.072*i/10)]

    shape1=numpy.array([.35, 1.])
    scales=[.82, 0.51, 0.39, 0.3, 0.22, 0.2, 0.17, 0.14, 0.11, 0.06]
    axl=[]
    for sc, xc, yc in zip(scales, ax_xc, ax_yc):
        w, l=shape1*sc
        axl+=[fig.add_axes([xc-w/2, yc-l/2, w, l])]


    stpl=[]
    xpos=[.27]*10
    xpos[0:3]=[.38, .36, .33]
    xpos[-1]=.18
    for count, (ax, xp) in enumerate(zip(axl, xpos)):
        stp=TernaryPlot(ax, ellabels=ellabels[:3], offset=.03)
        if not fontsize is None:
            stp.label(fontsize=fontsize)#,fontdict={'fontname':'Times New Roman'})
        stpl+=[stp]
        
        if not fontsize is None:
            if count<9:
                stp.ax.text(xp, .8, '%s$_{%.2f-%.2f}$' %(ellabels[3], (count*.1), ((count+1)*.1)-.01), ha='right', va='center', fontsize=fontsize)
            else:
                stp.ax.text(xp, .8, '%s$_{%.2f-%d}$' %(ellabels[3], (count*.1), 1), ha='right', va='center', fontsize=fontsize)
    return axl, stpl



def scatter_10axes(comps, fom, stpl, s=18, cb=False, cbrect=(.85, .3, .04, .4), cblabel='', **kwargs):# for colorbar must pass kwargs norm and cmap and optionally cblabel
    abc=comps[:, :3]
    abc[abc.sum(axis=1)==0.]=numpy.array([1., 1., 1.])/3.
    abc=numpy.array([c/c.sum() for c in abc])
    d=comps[:, 3]
    d30=numpy.round(d*30.)
    dlims=numpy.array([0., 1., 2., 3.])
    marks=[('o', 1., 1.), ('D', .9, .7),('s', .8, .5)]
    sl=s*numpy.array([2.3, 1., .55, .4, .4, .45, .4, .5, .8, 1., 1.5])
    scplots=[]
    for i, (stp, sv) in enumerate(zip(stpl, sl)):
        dl=dlims+(i*3.)
        if i==9:
            dl[-1]+=.01
        for a, b, (m, sf, al) in zip(dl, dl[1:], marks):
            inds=numpy.where((d30>=a) & (d30<b))[0]
            #print a, b, len(inds)
            if len(inds)>0:
                scplots+=[stp.scatter(abc[inds], c=fom[inds], marker=m, s=sv*sf, alpha=al, **kwargs)]
    if cb:
        cbax=stp.ax.figure.add_axes(cbrect)
        if 'extend' in kwargs.keys():
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'], extend=kwargs['extend'])
        else:
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'])
        sm.set_array(fom)
        cb=stp.ax.figure.colorbar(sm, cax=cbax)
        cb.set_label(cblabel, fontsize=18)


