import matplotlib.cm as cm
import numpy
import pylab
import h5py, operator, copy, os

#os.chdir('C:/Users/Gregoire/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot


def make10ternaxes(ellabels=['A', 'B', 'C', 'D']):

    ax_xc=[]
    ax_yc=[]
    xcdel=[.22, .2, .1, .04, .04, .04, .03, .03, .03, .03]
    ax_yc=[.5, .6, .36, .65, .49, .33, .68, .55, .43, .31, .2]
    for i in range(10):
        if i==0:
            ax_xc+=[xcdel[i]]
        else:
            ax_xc+=[ax_xc[-1]+xcdel[i]]
        #ax_yc+=[.5+((i%2)*2.-1.)*((i>0)*.1+.072*i/10)]

    shape1=numpy.array([.35, 1.])
    fig=pylab.figure(figsize=(12, 8))

    axl=[]
    for i, xc, yc in zip(range(1, 11), ax_xc, ax_yc):
        w, l=shape1/i
        axl+=[fig.add_axes([xc-w/2, yc-l/2, w, l])]


    stpl=[]
    for count, ax in enumerate(axl):
        stp=TernaryPlot(ax, ellabels=ellabels[:3])
        stp.label()
        stpl+=[stp]
        if count<9:
            stp.ax.text(.3, .8, '%s$_{%.2f-%.2f}$' %(ellabels[3], (count*.1), ((count+1)*.1)-.01), ha='right', va='center')
        else:
            stp.ax.text(.3, .8, '%s$_{%.2f-%d}$' %(ellabels[3], (count*.1), 1), ha='right', va='center')
    return axl, stpl



def scatter_10axes(comps, fom, stpl, **kwargs):
    abc=comps[:, :3]
    abc[abc.sum(axis=1)==0.]=numpy.array([1., 1., 1.])/3.
    abc=numpy.array([c/c.sum() for c in abc])
    d=comps[:, 3]
    for i, stp in enumerate(stpl):
        j=i+1
        if i==9:
            j+=1
        inds=numpy.where((d>=(i*.1)) & (d<(j*.1)))[0]
        if len(inds)>0:
            stp.scatter(abc[inds], c=fom[inds], **kwargs)


