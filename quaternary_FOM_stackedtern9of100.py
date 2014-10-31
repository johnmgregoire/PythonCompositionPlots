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


def make9of100ternaxes(ellabels=['A', 'B', 'C', 'D'], fig=None):
    if fig is None:
        fig=pylab.figure(figsize=(8, 6))
        
    axl=[]
    for i in [1, 4, 7, 2, 5, 8, 3, 6, 9]:
        
        axl+=[fig.add_subplot(3, 3, i)]
    fig.subplots_adjust(left=0.05, right=.8, bottom=.05, top=.95, hspace=.08, wspace=.08)

    stpl=[]
    xpos=[.27]*9

    for count, (ax, xp) in enumerate(zip(axl, xpos)):
        stp=TernaryPlot(ax, ellabels=ellabels[:3], offset=.03)
        stp.label(fontsize=15)#,fontdict={'fontname':'Times New Roman'})
        stpl+=[stp]
        
        
        #if count<4:
        stp.ax.text(xp, .8, '%s$_{%.2f}$' %(ellabels[3], (count*.01)), ha='right', va='center', fontsize=15)
        #else:
        #    stp.ax.text(xp, .8, '%s$_{%.2f-%d}$' %(ellabels[3], (count*.2), 1), ha='right', va='center', fontsize=17)
    return axl, stpl



def scatter_9of100axes(comps, fom, stpl, s=14, cb=False, cbrect=(.85, .3, .04, .4), cblabel='', **kwargs):# for colorbar must pass kwargs norm and cmap and optionally cblabel
    abc=comps[:, :3]
    abc[abc.sum(axis=1)==0.]=numpy.array([1., 1., 1.])/3.
    abc=numpy.array([c/c.sum() for c in abc])
    d=comps[:, 3]
    d100=numpy.round(d*100.)

    scplots=[]
    for i, (stp) in enumerate(stpl):
        inds=numpy.where(d100==i)[0]
        #print a, b, len(inds)
        if len(inds)>0:
            scplots+=[stp.scatter(abc[inds], c=fom[inds], marker='H', s=s, alpha=1, **kwargs)]
    if cb:
        cbax=stp.ax.figure.add_axes(cbrect)
        if 'extend' in kwargs.keys():
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'], extend=kwargs['extend'])
        else:
            sm=cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap'])
        sm.set_array(fom)
        cb=stp.ax.figure.colorbar(sm, cax=cbax)
        cb.set_label(cblabel, fontsize=18)


