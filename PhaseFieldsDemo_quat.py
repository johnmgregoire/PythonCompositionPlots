import pylab, numpy, copy, itertools
from myquaternaryutility import QuaternaryPlot
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def samesside(p,q,a,b):
    cp1=numpy.cross(b-a,p-a)
    cp2=numpy.cross(b-a,q-a)
    return numpy.dot(cp1,cp2)>0.
def intriangle(p, trianglepoints):
    return not (False in [samesside(p,*triperm) for triperm in itertools.permutations(trianglepoints)])


q=QuaternaryPlot(111)


#define these to be modified for each end member
z=numpy.zeros(4, dtype='float64')
ctr3=numpy.ones(4, dtype='float64')/3.
endmembers=[]
lineendpairs=[]
#iterate over 4 end members and draw a line from there to center of opposing face, e.g. (0,.33,.33,.33)
for i in range(4):
    a=copy.copy(z)
    a[i]=1.
    b=copy.copy(ctr3)
    b[i]=0.
    q.line(a, b, fmt='b-')
    q.scatter([b], c='b', s=15)
    endmembers+=[a]
    lineendpairs+=[[a, b]]
#convert the end members and pairs of endpts to cartesian
xyz_lineendpairs=[numpy.array(q.toCart(ls)).T for ls in lineendpairs]
xyz_endmembers=numpy.array(q.toCart(endmembers)).T

#choose the composition of a phase and draw the trivial phase field lines
phcomp=numpy.array([.1, .3, .2, .4])
q.scatter([phcomp], c='r', s=20)
for i in range(4):
    a=copy.copy(z)
    a[i]=1.
    q.line(a, phcomp, fmt='r-')

# iterate over all 4 phase field triangular boundaries (triangle defined by 3 points, the phase p0 and 2 end members p1,p2) and all 4 composition lines. find intersections
p0=numpy.array(q.toCart([phcomp])).T[0]
xyz_intr_dlist=[]
for countends, (p1, p2) in enumerate(itertools.combinations(xyz_endmembers, 2)):
    for countlines, (l0, l1) in enumerate(xyz_lineendpairs):
        t = numpy.cross(p1-p0, p2 - p0)
        #this is where the line intersects the plane, although it may not be within the triangle
        xyz_intr= l0 + (numpy.dot(t, p0 - l0) / numpy.dot(t, l1 - l0)) * (l1 - l0)
        #check within triangle phase field boundary and if so draw a green dot at itnersection and a translucent triangle
        xyz_triverts=numpy.array([p0, p1, p2])
        if intriangle(xyz_intr,xyz_triverts):
            xyz_intr_dlist+=[dict({}, xyz_intr=xyz_intr, xyz_triverts=xyz_triverts, xyz_lineends=(l0, l1), index_endmems=countends, index_lineends=countlines)]
            q.scatter(q.toComp([xyz_intr]), c='g', s=20)

            verts = [[p0, p1, p2]]

            collection = Poly3DCollection(verts, linewidths=0, alpha=0.2)
            face_color = [0.5, 0.5, 0.5] 
            collection.set_facecolor(face_color)
            q.ax.add_collection3d(collection)

q.label(ha='center', va='center', fontsize=16)
q.set_projection(azim=-17, elev=-6)
pylab.show()



