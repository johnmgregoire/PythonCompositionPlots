import matplotlib.cm as cm
import numpy, sys
import pylab
import operator, copy, os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from quaternary_folded_ternaries import ternaryfaces_folded
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

from myquaternaryutility import QuaternaryPlot


class plotwidget(FigureCanvas):
    def __init__(self, parent, width=12, height=6, dpi=72, projection3d=False):

        #plotdata can be 2d array for image plot or list of 2 1d arrays for x-y plot or 2d array for image plot or list of lists of 2 1D arrays
        self.projection3d=projection3d
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
        NavigationToolbar(self, self)

        self.mpl_connect('button_press_event', self.myclick)
        self.clicklist=[]
        self.cbax=None
    

    def myclick(self, event):
        if not (event.xdata is None or event.ydata is None):
            arrayxy=[event.xdata, event.ydata]
            print 'clicked on image: array indeces ', arrayxy, ' using button', event.button
            self.clicklist+=[arrayxy]
            self.emit(SIGNAL("genericclickonplot"), [event.xdata, event.ydata, event.button, event.inaxes])

class dialog(QDialog):
    def __init__(self, parent=None, title='', folderpath=None):
        super(dialog, self).__init__(parent)

        
        plotw=plotwidget(self)
        
        ax=plotw.axes
        
        intervs=20
        compsint=[[b, c, (intervs-a-b-c), a] for a in numpy.arange(0,intervs+1)[::-1] for b in numpy.arange(0,intervs+1-a) for c in numpy.arange(0,intervs+1-a-b)][::-1]
        print len(compsint)
        comps=numpy.float32(compsint)/intervs

        pylab.figure()
        stpquat=QuaternaryPlot(111)
        cols=stpquat.rgb_comp(comps)
        stpquat.scatter(comps, c=cols, s=100, edgecolors='none')
        stpquat.label()

        self.tf=ternaryfaces_folded(ax, nintervals=intervs)
        self.tf.label()
        self.tf.scatter(comps, cols, s='patch')
        
        
        QObject.connect(plotw, SIGNAL("genericclickonplot"), self.plotclick)
        
        mainlayout=QGridLayout()
        mainlayout.addWidget(plotw, 0, 0)

        
        self.setLayout(mainlayout)
    
    def plotclick(self, coords_button_ax):
        xc, yc, button, ax=coords_button_ax
        print self.tf.toComp(xc, yc)
        
class MainMenu(QMainWindow):
    def __init__(self):
        super(MainMenu, self).__init__(None)
        
        x=dialog()
        x.exec_()
        
mainapp=QApplication(sys.argv)
form=MainMenu()
form.show()
form.setFocus()
mainapp.exec_()

