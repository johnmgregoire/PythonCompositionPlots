import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

from quaternary_binary_lines import binarylines

gridi=30
comps_10full=[(a*1./gridi, b*1./gridi, c*1./gridi, (gridi-a-b-c)*1./gridi) for a in numpy.arange(0,1+gridi) for b in numpy.arange(0,1+gridi-a) for c in numpy.arange(0,1+gridi-a-b)]
comps_10full=list(set(comps_10full))
print len(comps_10full)
#plotpoints_cmyk
comps_10full=numpy.array(comps_10full)



fom=numpy.array([(c*numpy.arange(1, 5)).sum() for c in comps_10full])

fig=pylab.figure()
ax=fig.add_axes([.3, .12, .6, .83])#add_subplot(111)
insetax=fig.add_axes([0, .7, .2, .3], projection='3d')
bl=binarylines(ax, insetax)
bl.plotbinaryfom(comps_10full, fom)
bl.binarylineslegend(loc=7)

pylab.show()
