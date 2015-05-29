import time
import os, os.path
import sys
import numpy
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import operator
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar



from matplotlib.figure import Figure
import numpy.ma as ma
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import pylab
import pickle
#from fcns_math import *
#from fcns_io import *
#from fcns_ui import *
#PyCodePath=os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]

matplotlib.rcParams['backend.qt4'] = 'PyQt4'



wd=os.getcwd()


#sys.path.append(os.path.join(PyCodePath,'PythonCompositionPlots'))
from myternaryutility import TernaryPlot
from myquaternaryutility import QuaternaryPlot
from quaternary_FOM_stackedtern2 import *
from quaternary_FOM_stackedtern5 import *
from quaternary_FOM_stackedtern20 import *
from quaternary_FOM_stackedtern30 import *
from quaternary_FOM_stackedtern9of100 import *
from quaternary_ternary_faces import *
from quaternary_faces_shells import *
from quaternary_folded_ternaries import *

#        self.plotw_stack.fig.clf()
#        self.plotw_stack_axl, self.plotw_stack_stpl=makefcn(fig=self.plotw_stack.fig, ellabels=self.ellabels)
#        self.cbax_stack=self.plotw_stack.fig.add_axes(self.plotw_stack_cbaxrect)
        
        
class plottypeDialog(QDialog):
    def __init__(self, parent=None, ellabels=['A', 'B', 'C', 'D']):
        super(plottypeDialog, self).__init__(parent)
        mainlayout=QVBoxLayout()
        self.ellabels=ellabels
        self.plotw=plotwidget(self)
        self.plotw.fig.clf()
        mainlayout.addWidget(self.plotw)

        self.plottypeComboBox = QComboBox(self)
        mainlayout.addWidget(self.plottypeComboBox)
        QObject.connect(self.plottypeComboBox,SIGNAL("activated(QString)"),self.plot)
        
        self.setLayout(mainlayout)

        self.ternaryfaceoptions=[\
        ('tern. layered shells',  ternaryfaces_shells), \
        ('folded tern. layers', ternaryfaces_folded), \
        ('only tern. faces', ternaryfaces), \
        ]
        self.ternaryface_uiinds=[1, 2, 3]

        self.stackedternoptions=[\
        ('20% interv ternaries', (make5ternaxes, scatter_5axes)), \
        ('10% interv ternaries', (make10ternaxes, scatter_10axes)), \
        ('5% interv ternaries', (make20ternaxes, scatter_20axes)), \
        ('3.3% interv ternaries', (make30ternaxes, scatter_30axes)), \
        ('9 plots at 1% interv', (make9of100ternaxes, scatter_9of100axes)), \
        ]
        self.stackedtern_uiinds=[4, 5, 6, 7, 8]
        
        self.fillplotoptions()
    
    def fillplotoptions(self):
        self.plottypeComboBox.clear()
        self.plottypeComboBox.insertItem(0, 'none')
        for count, tup in enumerate(self.ternaryfaceoptions):
            self.plottypeComboBox.insertItem(999, tup[0])
        for count, tup in enumerate(self.stackedternoptions):
            self.plottypeComboBox.insertItem(999, tup[0])
        self.plottypeComboBox.setCurrentIndex(1)
    
    def loadplotdata(self, quatcomps, cols, nintervals=None):
        self.cols=cols
        self.quatcomps=quatcomps
        if nintervals is None:
            pairwisediffs=(((quatcomps[1:]-quatcomps[:-1])**2).sum(axis=1))**.5/2.**.5
            mindiff=(pairwisediffs[pairwisediffs>0.005]).min()
            self.nintervals=round(1./mindiff)
        else:
            self.nintervals=nintervals
            
    def plot(self, **kwargs):
        self.plotw.fig.clf()
        i=self.plottypeComboBox.currentIndex()
        if i in self.ternaryface_uiinds:
            selclass=self.ternaryfaceoptions[self.ternaryface_uiinds.index(i)][1]
            self.ternaryfaceplot(selclass, **kwargs)
        if i in self.stackedtern_uiinds:
            makefcn, scatterfcn=self.stackedternoptions[self.stackedtern_uiinds.index(i)][1]
            self.stackedternplot(makefcn, scatterfcn, **kwargs)
        self.plotw.fig.canvas.draw()
        
    def ternaryfaceplot(self, plotclass, **kwargs):
        ax=self.plotw.fig.add_axes((0, 0, 1, 1))
        tf=plotclass(ax, nintervals=self.nintervals)
        tf.label()
        tf.scatter(self.quatcomps, self.cols, **kwargs)
        
    def stackedternplot(self, makefcn, scatterfcn, **kwargs):
        self.axl, self.stpl=makefcn(fig=self.plotw.fig, ellabels=self.ellabels)
        scatterfcn(self.quatcomps, self.cols, self.stpl, edgecolor='none', **kwargs)
        



class plotwidget(FigureCanvas):
    def __init__(self, parent, width=12, height=6, dpi=72, projection3d=False):

        #plotdata can be 2d array for image plot or list of 2 1d arrays for x-y plot or 2d array for image plot or list of lists of 2 1D arrays

        self.fig=Figure(figsize=(width, height), dpi=dpi)
        if projection3d:
            self.axes=self.fig.add_subplot(111, navigate=True, projection='3d')
        else:
            self.axes=self.fig.add_subplot(111, navigate=True)

        self.axes.hold(True)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        #self.parent=parent
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        #NavigationToolbar(self, parent)
        
        self.toolbar=NavigationToolbar(self.figure.canvas, self)
        self.toolbar.setMovable(True)#DOESNT DO ANYTHING
        self.mpl_connect('button_press_event', self.myclick)
        self.clicklist=[]

    def myclick(self, event):
        if not (event.xdata is None or event.ydata is None):
            arrayxy=[event.xdata, event.ydata]
            print 'clicked on image: array indeces ', arrayxy, ' using button', event.button
            self.clicklist+=[arrayxy]
            self.emit(SIGNAL("genericclickonplot"), [event.xdata, event.ydata, event.button])

if __name__ == "__main__":
    class MainMenu(QMainWindow):
        def __init__(self, previousmm, execute=True, **kwargs):
            super(MainMenu, self).__init__(None)
            self.ui=plottypeDialog(self, **kwargs)
            
            intervs=20
            compsint=[[b, c, (intervs-a-b-c), a] for a in numpy.arange(0,intervs+1)[::-1] for b in numpy.arange(0,intervs+1-a) for c in numpy.arange(0,intervs+1-a-b)][::-1]
            print len(compsint)
            comps=numpy.float32(compsint)/intervs
            pylab.figure()
            stpquat=QuaternaryPlot(111)
            cols=stpquat.rgb_comp(comps)

            self.ui.loadplotdata(comps, cols)
            
            if execute:
                self.ui.exec_()
    mainapp=QApplication(sys.argv)
    form=MainMenu(None)
    form.show()
    form.setFocus()
    mainapp.exec_()

