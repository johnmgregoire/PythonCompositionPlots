import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

from quaternary_faces_shells import ternaryfaces_shells

from myquaternaryutility import QuaternaryPlot



intervs=20
compsint=[[b, c, (intervs-a-b-c), a] for a in numpy.arange(0,intervs+1)[::-1] for b in numpy.arange(0,intervs+1-a) for c in numpy.arange(0,intervs+1-a-b)][::-1]
print len(compsint)
comps=numpy.float32(compsint)/intervs




pylab.figure()
stpquat=QuaternaryPlot(111)
cols=stpquat.rgb_comp(comps)
stpquat.scatter(comps, c=cols, s=100, edgecolors='none')
stpquat.label()


pylab.figure()
ax=pylab.gca()
tf=ternaryfaces_shells(ax, nintervals=intervs)
tf.label()

#inds_x_y=tf.toCart(comps)
tf.scatter(comps, cols, skipinds=[0, 1, 2, 3], s='patch')


pylab.show()
