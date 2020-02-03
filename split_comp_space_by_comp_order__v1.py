import matplotlib.cm as cm
import numpy
import pylab
import operator, copy, os

from quaternary_faces_shells import ternaryfaces_shells

from myquaternaryutility import QuaternaryPlot



intervs=10
compsint=[[b, c, (intervs-a-b-c), a] for a in numpy.arange(0,intervs+1)[::-1] for b in numpy.arange(0,intervs+1-a) for c in numpy.arange(0,intervs+1-a-b)][::-1]
print len(compsint)
comps=numpy.float32(compsint)/intervs




pylab.figure()
stpquat=QuaternaryPlot(111)
cols=stpquat.rgb_comp(comps)



pylab.figure()
ax=pylab.gca()
tf=ternaryfaces_shells(ax, nintervals=intervs, outlinealpha=1)
#tf.label()

#inds_x_y=tf.toCart(comps)

if 1:#only quat comps
    plotinds=numpy.where((comps>0.).prod(axis=1))[0]
elif 0:#only A,B,C ternnaryComps
    plotinds=numpy.where((comps[:, :3]>0.).prod(axis=1)*(comps[:, 3]==0.))[0]
elif 0:#only A,C binaryComps
    plotinds=numpy.where((comps[:, 0]>0.)*(comps[:, 1]==0.)*(comps[:, 2]>0.)*(comps[:, 3]==0.))[0]

comps=comps[plotinds]
cols=cols[plotinds]


nshell_inds_x_y=tf.get_plot_coords(comps)


patchfcn=lambda x, y, c:tf.patch_xyc(x, y, c)

output_table=[]
for nshell, inds, x, y in nshell_inds_x_y:
    map(patchfcn, x, y, cols[inds])
    for xv, yv, i in zip(x, y, inds):
        output_table+=[list(comps[i])+[nshell, xv, yv]]#if only exporting quat comps then no duplicates

lines=['A,B,C,D,nshell,x,y']
lines+=['%.2f,%.2f,%.2f,%.2f,%d,%.6f,%.6f' %tuple(tup) for tup in output_table]

print '\n'.join(lines)




pylab.show()
