import pylab, numpy
from myquaternaryutility import QuaternaryPlot

q=QuaternaryPlot(111)
#t=numpy.linspace(0,1.,5)
#comps=[[a,b,c,d] for a in t for b in t for c in t for d in t if a+b+c+d==1.]
#comps=numpy.float32(comps)

t=numpy.linspace(0,1.,30)

comps=[[a,b,1.-a-b-(2.*a**2+b),2.*a**2+b] for a in t for b in t[:10] if a+b+(2.*a**2+b)<=1.]
comps=numpy.float32(comps)

x, y, z=q.toCart(comps)

q.scatter(comps,c=comps[:,3])
#q.ax.scatter(x, y, z, c=z)

q.plotabcprojection(comps, c=(.4, .4, .4))
q.label(ha='center', va='center', fontsize=16)
pylab.show()
