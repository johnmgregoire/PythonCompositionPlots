import matplotlib.cm as cm
import numpy
import pylab
import h5py, operator, copy, os

from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot

from quaternary_FOM_stackedtern2 import *

axl, stpl=make10ternaxes()

gridi=30
comps_10full=[(a*1./gridi, b*1./gridi, c*1./gridi, (gridi-a-b-c)*1./gridi) for a in numpy.arange(0,1+gridi) for b in numpy.arange(0,1+gridi-a) for c in numpy.arange(0,1+gridi-a-b)]
comps_10full=list(set(comps_10full))
print len(comps_10full)
#plotpoints_cmyk
comps_10full=numpy.array(comps_10full)

pylab.figure()
stpquat=QuaternaryPlot(111)
cols=stpquat.rgb_comp(comps_10full)
stpquat.scatter(comps_10full, c=cols, s=20, edgecolors='none')
scatter_10axes(comps_10full, cols, stpl, s=20, edgecolors='none', cmap=cm.jet, norm=None, cb=True)
stpquat.label()

pylab.savefig('stackedtern_quat.png')
pylab.figure(axl[0].figure.number)
pylab.savefig('stackedtern.png')

pylab.show()
