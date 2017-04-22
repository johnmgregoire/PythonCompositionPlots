import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

from quaternary_faces_shells import ternaryfaces_shells

from myquaternaryutility import QuaternaryPlot



intervs=10
compsint=[[b, c, (intervs-a-b-c), a] for a in numpy.arange(0,intervs+1)[::-1] for b in numpy.arange(0,intervs+1-a) for c in numpy.arange(0,intervs+1-a-b)][::-1]
print len(compsint)
comps=numpy.float32(compsint)/intervs




pylab.figure()
stpquat=QuaternaryPlot(111)
cols=stpquat.rgb_comp(comps)
stpquat.scatter(comps, c=cols, s=1200, edgecolors='none')
stpquat.label()


pylab.figure()
ax=pylab.gca()
tf=ternaryfaces_shells(ax, nintervals=intervs)
tf.label()

#inds_x_y=tf.toCart(comps)
tf.scatter(comps, cols, skipinds=[0, 1, 2, 3], s='patch')

#pylab.figure(figsize=(12, 4))
#tf.quatscatter(comps, cols, s=200, fontsize=0, azim=72, elev=20, edgecolor='none', outline=True)
#pylab.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
#pylab.savefig('//htejcap.caltech.edu/share/home/users/hte/catalysts on BVO/IJonFTOandBVOsummaries/compdemo.svg')

pylab.figure(figsize=(12, 4))
#tf.quatscatter(comps, cols, s=200, fontsize=0, azim=72, elev=20, edgecolor='none', outline=False, alphaall=0.1, alphashell=1)
pylab.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
tf.quatplot3D(comps, cols, marker='o', markersize=14, fontsize=0, azim=72, elev=20, markeredgecolor='none', outline=False, alphaall=0.1, alphashell=1)
#pylab.savefig('//htejcap.caltech.edu/share/home/users/hte/catalysts on BVO/IJonFTOandBVOsummaries/compdemo.png', dpi=400)
pylab.show()
