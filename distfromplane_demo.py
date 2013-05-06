import pylab, numpy
from myquaternaryutility import QuaternaryPlot

q=QuaternaryPlot(211)
q2=QuaternaryPlot(212)
#t=numpy.linspace(0,1.,5)
#comps=[[a,b,c,d] for a in t for b in t for c in t for d in t if a+b+c+d==1.]
#comps=numpy.float32(comps)

t=numpy.linspace(0,1.,30)

comps=[[a,b,1.-a-b-(2.*a**2+b),2.*a**2+b] for a in t for b in t[:10] if a+b+(2.*a**2+b)<=1.]
comps=numpy.float32(comps)

examplenum=0

if examplenum==0:
    compvert2=numpy.array([0.125, .125, .6, .15])
    compvert0=numpy.array([.2, .2, 0., .6])
    compvert1=numpy.array([1., 0., 0., 0])
    critdist=.04
    withintriangle=False
elif examplenum==1:
    compvert2=numpy.array([0.125, .125, .6, .15])
    compvert0=numpy.array([.2, .2, 0., .6])
    compvert1=numpy.array([1., 0., 0., 0])
    critdist=.04
    withintriangle=True
    
q.scatter(comps,c=comps[:,3])

q.label(ha='center', va='center', fontsize=16)
q.set_projection(azim=-17, elev=-6)

inds, distfromplane, xyparr, xyp_verts,intriangle=q2.filterbydistancefromplane(comps, compvert0, compvert1, compvert2, critdist, withintriangle=withintriangle, invlogic=False, returnall=True)
indsnot=q2.filterbydistancefromplane(comps, compvert0, compvert1, compvert2, critdist, withintriangle=withintriangle, invlogic=True)
print len(inds), ' points'
q2.scatter(comps[inds],c=comps[inds,3])
q2.scatter(comps[indsnot],c='grey', marker='.', s=5)
q2.line(compvert0, compvert1)
q2.line(compvert1, compvert2)
q2.line(compvert2, compvert0)


q2.label(ha='center', va='center', fontsize=16)
q2.set_projection(azim=-17, elev=-6)

pylab.figure()
ax=pylab.subplot(111)

q2.plotfominselectedplane(ax, xyparr[inds], comps[inds, -1], xyp_verts=xyp_verts, vertcomps_labels=[compvert0, compvert1, compvert2], s=20)


pylab.show()
