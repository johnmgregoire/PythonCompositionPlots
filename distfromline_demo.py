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

examplenum=1

if examplenum==0:
    compend1=numpy.array([.2, .4, .3, .1])
    compend2=numpy.array([.2, 0, .2, .6])
    critdist=.1
    betweenpoints=False
elif examplenum==1:
    compend1=numpy.array([0.1, .1, .8, 0])
    compend2=numpy.array([.2, .2, 0., .6])
    critdist=.05
    betweenpoints=False
elif examplenum==2:
    compend1=numpy.array([0.125, .125, .6, .15])
    compend2=numpy.array([.2, .2, 0., .6])
    critdist=.05
    betweenpoints=True
elif examplenum==3:
    compend1=numpy.array([0.125, .125, .6, .15])
    compend2=numpy.array([.2, .2, 0., .6])
    critdist=.05
    betweenpoints=False

q.scatter(comps,c=comps[:,3])

q.label(ha='center', va='center', fontsize=16)
q.set_projection(azim=-17, elev=-6)

inds, distfromlin, lineparameter=q2.filterbydistancefromline(comps, compend1, compend2, critdist, betweenpoints=betweenpoints, invlogic=False, returnall=True)
indsnot=q2.filterbydistancefromline(comps, compend1, compend2, critdist, betweenpoints=betweenpoints, invlogic=True)
print len(inds), ' points'
q2.scatter(comps[inds],c=comps[inds,3])
q2.scatter(comps[indsnot],c='grey', marker='.', s=5)
q2.line(compend1, compend2)


q2.label(ha='center', va='center', fontsize=16)
q2.set_projection(azim=-17, elev=-6)

pylab.figure()

lp=lineparameter[inds]
argsinds=numpy.argsort(lp)
pylab.plot(lp[argsinds], comps[inds,3][argsinds], 'b-')
pylab.plot(lp[argsinds], comps[inds,3][argsinds], 'bo')
ax=pylab.gca()
lineparticks=numpy.linspace(0, 1, 4)
tl=[]
for i in lineparticks:
    c=compend1+(compend2-compend1)*i
    tl+=[q.singlelabeltext(c)]
ax.xaxis.set_ticks(lineparticks)
ax.xaxis.set_ticklabels(tl)
pylab.show()
