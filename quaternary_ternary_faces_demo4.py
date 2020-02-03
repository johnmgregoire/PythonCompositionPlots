import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

from quaternary_ternary_faces import ternaryfaces

from myquaternaryutility import QuaternaryPlot



gridi=10
comps_10full=[(a*1./gridi, b*1./gridi, c*1./gridi, (gridi-a-b-c)*1./gridi) for a in numpy.arange(0,1+gridi) for b in numpy.arange(0,1+gridi-a) for c in numpy.arange(0,1+gridi-a-b)]
comps_10full=[c for c in comps_10full if c[1]>=0.5 and c[0]>0 and c[2]>0 and c[3]==0]
comps_10full=list(set(comps_10full))
print len(comps_10full)


#comps_10full=[[1., 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

comps_10full=numpy.array(comps_10full)


pylab.figure()
stpquat=QuaternaryPlot(111)
#cols=stpquat.rgb_comp(comps_10full)




pylab.figure()
ax=pylab.gca()
tf=ternaryfaces(ax, nintervals=gridi)
tf.label()

#in first figure
cols=tf.ternaryplot.rgb_comp(comps_10full[:, [1, 0, 2]])
stpquat.scatter(comps_10full, c=cols, s=20, edgecolors='none')
stpquat.label()

#now second figure
#inds_x_y=tf.toCart(comps_10full)
tf.scatter(comps_10full, cols, skipinds=[0, 1, 2, 3], s='patch')


pylab.show()
