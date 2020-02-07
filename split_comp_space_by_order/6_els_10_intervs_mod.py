import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import CirclePolygon
from matplotlib.collections import PatchCollection
from itertools import combinations as comb
from collections import defaultdict
from glob import glob

ANALYSIS_DIR = 'L:/processes/analysis/eche'
plate_ids = ['5052', '5070', '5151', '5137', '5156', '5162', '5176']
pre_anas = ['20191126.143112', '20191126.143939', '20191126.144647', '20191126.145624', '20191126.153131', '20191126.160608', '20191126.164123']
pre_ints = [2, 4]
post_anas = ['20191126.171630', '20191126.173907', '20191126.180144', '20191126.182425', '20191126.184655', '20191203.090307', '20191203.090925']
post_ints = [2, 4, 6]
ana_fomnames = ['Jmin.mAcm2']

intervs=10


make_patch_xy=lambda x, y, **kwargs:CirclePolygon((x, y),radius=1./intervs/3.**.5,resolution=6, **kwargs)

d={}

d[1]={}
d[1]['comps']=np.int32([[intervs,0,0,0]])
d[1]['xy']=np.array([[0.,0.]])

df = pd.read_csv(r'ACcomps.csv',delimiter=',')
d[2]={}
d[2]['comps']=np.int32(np.round(intervs*np.array([df[el].values for el in ['A','C','B','D']]))).T#switch b and c so 1st 2 channels used
d[2]['xy']=np.array([df[el].values for el in ['x','y']]).T
d[2]['xy'][:,1]=0. #the points read are form ternary plot so set y to 0

df = pd.read_csv(r'ABCcomps.csv',delimiter=',')
d[3]={}
d[3]['comps']=np.int32(np.round(intervs*np.array([df[el].values for el in ['A','B','C','D']]))).T
d[3]['xy']=np.array([df[el].values for el in ['x','y']]).T

df = pd.read_csv(r'ABCDcomps.csv',delimiter=',')
d[4]={}
d[4]['comps']=np.int32(np.round(intervs*np.array([df[el].values for el in ['A','B','C','D']]))).T
d[4]['xy']=np.array([df[el].values for el in ['x','y']]).T

compsint=[[ee, dd, b, c, (intervs-a-b-c-dd-ee), a] for a in np.arange(0,intervs+1)[::-1] for b in np.arange(0,intervs+1-a) for c in np.arange(0,intervs+1-a-b) for dd in np.arange(0,intervs+1-a-b-c) for ee in np.arange(0,intervs+1-a-b-c-dd) if [ee, dd, b, c, (intervs-a-b-c-dd-ee), a].count(0)>=2] [::-1]


#larger fontsize
xystarts_comporder=[[]]#empty list of el order 0
xystarts_comporder+=[[[.5,b] for b in 4.5-0.9*np.arange(6) ]]
xystarts_comporder+=[[[1.2,b] for b in 4.5-0.4*np.arange(15) ]]
xystarts_comporder+=[[[a,b] for a in 2.5+np.arange(4)*1.1 for b in 4.3-.95*np.arange(5)]]
xystarts_comporder+=[[[a,b] for a in 5.+np.arange(3)*2.4 for b in 4.3-.95*np.arange(5)]]

dlab=0.31
labposns=[[]]
labposns+=[[[-0.25,0]]]
labposns+=[[[-0.12-dlab,0],[-0.12,0]]]
labposns+=[[[-0.15+a,.12+b] for a,b in [[0,0],[dlab/2.,dlab*.84],[dlab,0]]]]
labposns+=[[[1.85+a,0.12+b] for a,b in [[0,0],[dlab/2.,dlab*.84],[dlab,0],[3.*dlab/2.,dlab*.84]]]]

def f_setlim(ax,f=False):
    ax.set_xlim(0,13.9)
    ax.set_ylim(-1.5,4.9)
    if not f:
        #ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])

        
# for pid, ana_timestamp in zip(plate_ids, pre_anas):
#     for ana_int in pre_ints:
for pid, ana_timestamp in zip(plate_ids, post_anas):
    for ana_int in post_ints:
        ana_csvpath = glob("%s/%s*/ana__%i*.csv" %(ANALYSIS_DIR, ana_timestamp, ana_int))[0]
        fomdf = pd.read_csv(ana_csvpath, skiprows=8)
        compinds = [i for i,x in enumerate(fomdf.columns) if ".PM.AtFrac" in x]
        colnames = [fomdf.columns[i] for i in compinds]
        ana_elements = [x.replace(".PM.AtFrac", "") for x in colnames]
        fomdf[colnames] = fomdf[colnames].apply(lambda x: pd.Series.round(x, 1) * intervs)
        ana_compsint = fomdf[colnames].astype('int').values.tolist()
        aggdf = fomdf.groupby(colnames).agg({k: 'mean' for k in ana_fomnames}).reset_index()
        aggfomdict = {x: aggdf[x].values for x in ana_fomnames}
        agg_compsint = aggdf[colnames].astype('int').values.tolist()
        
        dupedict = {}
        for k in ana_fomnames:
            dupedict[k] = defaultdict(list)
            for fomv, compv in zip(fomdf[k].values, ana_compsint):
                dupedict[k][tuple(compv)].append(fomv)

        for fomk in ana_fomnames:
            dupedict[fomk] = {k: v for k,v in dupedict[fomk].items() if len(v)>1}
            
        els = ana_elements
        allelinds=range(6)
        allcomp_xy={}
        alllabs_xyel=[]
        for elorder in range(1,5):
            for combcount,elinds in enumerate(comb(allelinds,elorder)):
                xy0=np.array(xystarts_comporder[elorder][combcount])
                alllabs_xyel+=[[xy0[0]+a,xy0[1]+b,els[eli]] for (a,b),eli in zip(labposns[elorder],elinds)]
                dord=d[elorder]
                for c,xy in zip(dord['comps'],dord['xy']):
                    ct=np.zeros(6,'int32')
                    ct[list(elinds)]=c[:elorder]
                    allcomp_xy[tuple(ct)]=xy+xy0

        def addlabels(ax,fontsize=12,**kwargs):
            for x,y,s in alllabs_xyel:
                ax.text(x,y,s,ha='center',va='center',fontsize=fontsize,**kwargs)

        allxy=np.array(allcomp_xy.values())

        for fomname in ana_fomnames:
            plt.figure(figsize=(13,6.5))
            ax=plt.gca()
            ax.set_aspect(1)
            fom = aggfomdict[fomname]
            patches = []
            dupepatches = []
            dupefoms = []
            for c,z in zip(agg_compsint,fom):
                x,y=allcomp_xy[tuple(c)]
                patches.append(make_patch_xy(x, y))
                if tuple(c) in dupedict[fomname].keys():
                    for i,v in enumerate(dupedict[fomname][tuple(c)]):
                        dupepatches.append(make_patch_xy(x, y-0.1*(i+1)))
                        dupefoms.append(v)
            dupefoms = np.array(dupefoms)
            p = PatchCollection(patches, cmap=mpl.cm.viridis_r)
            dp = PatchCollection(dupepatches, cmap=mpl.cm.viridis_r)
            p.set_array(fom)
            dp.set_array(dupefoms)
            allfoms = np.concatenate((fom, dupefoms))
            clim = [np.percentile(allfoms, 2), np.percentile(allfoms, 98)]
            if 'Jmin' in fomname:
                p.set_clim(clim)
                dp.set_clim(clim)
            ax.add_collection(p)
            ax.add_collection(dp)
            f_setlim(ax)
            addlabels(ax)
#             plt.savefig('%s-prePETS-ana__%s-%s_heatmap.eps' %(pid, ana_int, fomname))
#             plt.savefig('%s-prePETS-ana__%s-%s_heatmap.png' %(pid, ana_int, fomname))
            plt.savefig('%s-postPETS-ana__%s-%s_heatmap.eps' %(pid, ana_int, fomname))
            plt.savefig('%s-postPETS-ana__%s-%s_heatmap.png' %(pid, ana_int, fomname))
            plt.close()

        lines=['x,y,A,B,C,D,E,F']
        lines+=['%.5f,%.5f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f' %tuple(list(v)+list(np.float32(k)/intervs)) for k,v in allcomp_xy.items()]
        comp_points_filestr='\n'.join(lines)
        with open('comp_points__%d_intervals.csv' %(intervs),mode='w') as f:
            f.write(comp_points_filestr)

        lines=['x,y,lab']
        lines+=['%.5f,%.5f,%s' %tuple(tup) for tup in alllabs_xyel]
        el_labels_filestr='\n'.join(lines)
        with open('label_posns__%d_intervals.csv' %(intervs),mode='w') as f:
            f.write(el_labels_filestr)

