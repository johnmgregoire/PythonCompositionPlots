from myternaryutility import TernaryPlot
import matplotlib.cm as cm
import numpy
import pylab, copy
from colorsys import rgb_to_hsv
from colorsys import hsv_to_rgb

pylab.figure(figsize=(6, 3))
ax=pylab.gca()
#stp = TernaryPlot(ax, ellabels=['Au', 'Cu', 'Si']) 
stp = TernaryPlot(ax, ellabels=['A', 'B', 'C'])
stp.grid(nintervals=10, printticklabels=[4])
stp.label(fontsize=12)

comps=numpy.random.rand(50, 3)
comps/=comps.sum(axis=1)[:, numpy.newaxis]
#compdist=(numpy.random.rand(len(comps), 3)-0.5)/5
comps2=copy.copy(comps)
comps2[:, 2]+=.5
comps2/=comps2.sum(axis=1)[:, numpy.newaxis]


#compsdiff=comps2-comps
#
#terncoord=numpy.float64(comps)
#terncoord2=numpy.float64(comps2)
#
#
#
#sat = ((compsdiff**2).sum(axis=1)/2.)**.5
#
#huelist=[0. if cd.sum()==0. else rgb_to_hsv(*(cd/cd.sum()))[0] for cd in numpy.abs(compsdiff)]
#
#sat_norm=sat/max(sat)
#
#rgblist=[hsv_to_rgb(h, s, 1) for h, s in zip(huelist, sat_norm)]

#rgb_arr=stp.complex_to_rgb(ang, sat_norm)


stp.hsdiffplot(comps, comps2)
#
#
pylab.show()

print 'done'
