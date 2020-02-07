import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import CirclePolygon
from matplotlib.collections import PatchCollection
from itertools import combinations as comb
from glob import glob

ANALYSIS_DIR = 'L:/processes/analysis/eche'
# ana_timestamp = '20191126.143112'
plate_ids = ['5052', '5070', '5151', '5137', '5156', '5162', '5176']
pre_anas = ['20191126.143112', '20191126.143939', '20191126.144647', '20191126.145624', '20191126.153131', '20191126.160608', '20191126.164123']
pre_ints = [2, 4]
post_anas = ['20191126.171630', '20191126.173907', '20191126.180144', '20191126.182425', '20191126.184655', '20191203.090307', '20191203.090925']
post_ints = [2, 4, 6]
ana_fomnames = ['Jmin.mAcm2']

# for pid, ana_timestamp in zip(plate_ids, pre_anas):
    # for ana_int in pre_ints:
for pid, ana_timestamp in zip(plate_ids, post_anas):
    for ana_int in post_ints:
        ana_csvpath = glob("%s/%s*/ana__%i*.csv" %(ANALYSIS_DIR, ana_timestamp, ana_int))[0]
        df = pd.read_csv(ana_csvpath, skiprows=8)
        fomdict = {x: df[x].values for x in ana_fomnames}

        compinds = [i for i,x in enumerate(df.columns) if ".PM.AtFrac" in x]
        colnames = [df.columns[i] for i in compinds]
        ana_elements = [x.replace(".PM.AtFrac", "") for x in colnames]
        ana_compsint = (df[colnames].round(1)*10).astype('int').values.tolist()


        intervs=10
        # patch_xyc=lambda ax, x, y, c, **kwargs:ax.add_patch(CirclePolygon((x, y),radius=1./intervs/3.**.5,resolution=6, color=c, **kwargs))
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
        #for nzeros in range(6):
        #    compsint=[[ee, dd, b, c, (intervs-a-b-c-dd-ee), a] for a in np.arange(0,intervs+1)[::-1] for b in np.arange(0,intervs+1-a) for c in np.arange(0,intervs+1-a-b) for dd in np.arange(0,intervs+1-a-b-c) for ee in np.arange(0,intervs+1-a-b-c-dd) if [ee, dd, b, c, (intervs-a-b-c-dd-ee), a].count(0)==nzeros] [::-1]
        #    print len(compsint)
        compsint

        #compact
        xystarts_comporder=[[]]#empty list of el order 0
        xystarts_comporder+=[[[.4,b] for b in 3.9-.2*np.arange(6) ]]
        xystarts_comporder+=[[[.95,b] for b in 3.9-.2*np.arange(15) ]]
        xystarts_comporder+=[[[a,b] for a in 2.1+np.arange(4) for b in 3.7-.8*np.arange(5)]]
        xystarts_comporder+=[[[a,b] for a in 4.3+np.arange(3)*2.3 for b in 3.7-.8*np.arange(5)]]

        dlab=0.18
        labposns=[[]]
        labposns+=[[[-0.2,0]]]
        labposns+=[[[-0.1-dlab,0],[-0.1,0]]]
        labposns+=[[[-0.02+a,.1+b] for a,b in [[0,0],[dlab/2.,dlab/2.*3.**.5],[dlab,0]]]]
        labposns+=[[[1.98+a,0.1+b] for a,b in [[0,0],[dlab/2.,dlab/2.*3.**.5],[dlab,0],[3.*dlab/2.,dlab/2.*3.**.5]]]]

        # els='Ni,Fe,La,Ce,Co,Mn'.split(',')
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

        def addlabels(ax,fontsize=8,**kwargs):
            for x,y,s in alllabs_xyel:
                ax.text(x,y,s,ha='center',va='center',fontsize=fontsize,**kwargs)


        allxy=np.array(allcomp_xy.values())
        #def f_setlim(ax):
        #    ax.set_xlim(min(allxy[:,0])-0.1,max(allxy[:,0])+0.1)
        #    ax.set_ylim(min(allxy[:,1])-0.1,max(allxy[:,1])+0.1)
        def f_setlim(ax):
            ax.set_xlim(0,13)
            ax.set_ylim(0,4.15)


        #larger fontsize
        xystarts_comporder=[[]]#empty list of el order 0
        xystarts_comporder+=[[[.5,b] for b in 4.5-.25*np.arange(6) ]]
        xystarts_comporder+=[[[1.2,b] for b in 4.5-.25*np.arange(15) ]]
        xystarts_comporder+=[[[a,b] for a in 2.5+np.arange(4)*1.1 for b in 4.3-.95*np.arange(5)]]
        xystarts_comporder+=[[[a,b] for a in 5.+np.arange(3)*2.4 for b in 4.3-.95*np.arange(5)]]

        dlab=0.31
        labposns=[[]]
        labposns+=[[[-0.25,0]]]
        labposns+=[[[-0.12-dlab,0],[-0.12,0]]]
        labposns+=[[[-0.15+a,.12+b] for a,b in [[0,0],[dlab/2.,dlab*.84],[dlab,0]]]]
        labposns+=[[[1.85+a,0.12+b] for a,b in [[0,0],[dlab/2.,dlab*.84],[dlab,0],[3.*dlab/2.,dlab*.84]]]]

        # els='Ni,Fe,La,Ce,Co,Mn'.split(',')
        # els='A,B,C,D,E,F'.split(',')
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
        #def f_setlim(ax):
        #    ax.set_xlim(min(allxy[:,0])-0.1,max(allxy[:,0])+0.1)
        #    ax.set_ylim(min(allxy[:,1])-0.1,max(allxy[:,1])+0.1)
        def f_setlim(ax,f=False):
            ax.set_xlim(0,13.9)
            ax.set_ylim(0,4.9)
            if not f:
                #ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])

        compsint
        # for fomind in range(6):
        for fomname in ana_fomnames:
            plt.figure(figsize=(13,5))
            ax=plt.gca()
            ax.set_aspect(1)
            # fom=np.float64([c[fomind] for c in compsint])/intervs
            fom = fomdict[fomname]
            patches = []
            for c,z in zip(ana_compsint,fom):
                x,y=allcomp_xy[tuple(c)]
                patches.append(make_patch_xy(x, y))
                # patch_xyc(ax,x,y,[z,0.1,0.1])
            p = PatchCollection(patches, cmap=mpl.cm.viridis_r)
            p.set_array(fom)
            # clim = p.get_clim()
            if 'Jmin' in fomname:
                p.set_clim([np.percentile(fom, 2), np.percentile(fom, 98)])
            ax.add_collection(p)
            f_setlim(ax)
            addlabels(ax)
            # plt.savefig('%s-prePETS-ana__%s-%s_heatmap.eps' %(pid, ana_int, fomname))
            # plt.savefig('%s-prePETS-ana__%s-%s_heatmap.png' %(pid, ana_int, fomname))
            plt.savefig('%s-postPETS-ana__%s-%s_heatmap.eps' %(pid, ana_int, fomname))
            plt.savefig('%s-postPETS-ana__%s-%s_heatmap.png' %(pid, ana_int, fomname))

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

#plt.figure()
#ax=plt.gca()
#ax.set_aspect(1)
#fom=np.float64([c.count(0) for c in compsint])/5.
#for c,z in zip(compsint,fom):
#    x,y=allcomp_xy[tuple(c)]
#    patch_xyc(ax,x,y,[z,0.1,0.1])
#f_setlim(ax)
#
