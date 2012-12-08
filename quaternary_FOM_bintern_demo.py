import matplotlib.cm as cm
import numpy
import pylab
import h5py, operator, copy, os

from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot

from quaternary_FOM_bintern import *



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
stpquat.label()

axl, stpl=make4ternaxes()
scatter_4axes(comps_10full, cols, stpl, edgecolors='none')
#fom=numpy.random.rand(len(comps_10full))
fom=numpy.array([(c*numpy.arange(1, 5)).sum() for c in comps_10full])

#pylab.figure()
#axbin=pylab.subplot(111)
axbin, axbininset=plotbinarylines_axandinset(linewidth=2)
plotbinarylines_quat(axbin, comps_10full, fom, markersize=10)

pylab.figure(stpquat.ax.figure.number)
pylab.savefig('bintern_quat.png')
pylab.figure(axl[0].figure.number)
pylab.savefig('ternfaces.png')
pylab.figure(axbin.figure.number)
pylab.savefig('binlines.png')

pylab.show()
