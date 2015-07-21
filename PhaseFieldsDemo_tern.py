import pylab, numpy, copy, itertools
from myternaryutility import TernaryPlot


def interwithinsegs(p0, p1, p2, p3):
    x0, y0=p0
    x1, y1=p1
    x2, y2=p2
    x3, y3=p3
    d=(x0-x1)*(y2-y3)-(y0-y1)*(x2-x3);
    x=((x2-x3)*(x0*y1-y0*x1)-(x0-x1)*(x2*y3-y2*x3))/d
    y=((y2-y3)*(x0*y1-y0*x1)-(y0-y1)*(x2*y3-y2*x3))/d
    
    betweentest=lambda a, b, c: (a>=min(b, c)) and (a<=max(b, c))
    return numpy.array([x, y]), betweentest(x, x0, x1) and betweentest(x, x2, x3) and betweentest(y, y0, y1) and betweentest(y, y2, y3)
    


q=TernaryPlot(111)


#define these to be modified for each end member
z=numpy.zeros(3, dtype='float64')
ctr2=numpy.ones(3, dtype='float64')/2.
endmembers=[]
lineendpairs=[]
#iterate over 4 end members and draw a line from there to center of opposing face, e.g. (0,.33,.33,.33)
for i in range(3):
    a=copy.copy(z)
    a[i]=1.
    b=copy.copy(ctr2)
    b[i]=0.
    q.line(a, b, fmt='b-')
    q.scatter([b], c='b', s=15)
    endmembers+=[a]
    lineendpairs+=[[a, b]]
#convert the end members and pairs of endpts to cartesian
xy_lineendpairs=[numpy.array(q.toCart(ls)).T for ls in lineendpairs]
xy_endmembers=numpy.array(q.toCart(endmembers)).T

#choose the composition of a phase and draw the trivial phase field lines
phcomp=numpy.array([.5, .3, .2])
q.scatter([phcomp], c='r', s=20)
for i in range(3):
    a=copy.copy(z)
    a[i]=1.
    q.line(a, phcomp, fmt='r-')

# iterate over all 4 phase field triangular boundaries (triangle defined by 3 points, the phase p0 and 2 end members p1,p2) and all 4 composition lines. find intersections
p0=numpy.array(q.toCart([phcomp])).T[0]
xy_intr_dlist=[]
for countends, (p1) in enumerate(xy_endmembers):
    for countlines, (l0, l1) in enumerate(xy_lineendpairs):

        xy_intr, interbool=interwithinsegs(p0, p1, l0, l1)
        if interbool:
            xy_intr_dlist+=[dict({}, xy_intr=xy_intr, xy_lineends=(l0, l1), index_endmems=countends, index_lineends=countlines)]
            q.scatter(q.toComp([xy_intr]), c='g', s=20)



q.label(fontsize=16)
pylab.show()



