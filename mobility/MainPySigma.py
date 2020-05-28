# -*- coding: utf-8 -*-
"""
Prepared by Nicholas Peraino
Lumigen Instrument Center
Wayne State Univeristy


This is a rewrite of "Sigma"

   T. Wyttenbach, G. von Helden,... M.T. Bowers; JASMS 1997, 8,
   275 
   
via the FORTRAN code provided by Thomas Wyttenbach.

The goal is to make a community accessible module for advancement of ion-mobility research.


From the readme file:
    
Theory
------
The algorithm for the cross section calculations has been
described in the literature and it is sometimes referred to as
the "projection model":
    
   G. von Helden, M.-T. Hsu, N.G. Gotts, M.T. Bowers; JPC 1993,
   97, 8182
   T. Wyttenbach, G. von Helden,... M.T. Bowers; JASMS 1997, 8,
   275 (section "Model", last paragraph):
       
"The molecule is projected onto a randomly chosen plane in space,
and a circle with the corresponding collision radius is drawn
in that plane at the position of each projected atom. Then points
in the plane are randomly picked within a square of area A that
encloses the projected molecule. If a selected point is inside
one or more circles, it is counted as a hit. The ratio of hits to
the number of tries multiplied by A is the collision cross section
of that particular projection. The procedure is repeated for many
different randomly selected projections..."

The radius drawn around each atom is either taken from a parameter
file (hard sphere mode) or is calculated from Lennard-Jones
parameters supplied by a second parameter file (Lennard-Jones
mode). In the latter case the radii are temperature dependent.
"""

from PySigmaFunctions import *
import PySigmaFunctions as psf
import tkinter as tk
from tkinter.filedialog import askopenfilename
import math
import potential
from projectionApproximation import *

#Ask to open the input geometry commented out to save development time
#root = tk.Tk()
#root.withdraw()
#geometryInputFile = askopenfilename()

mol = molecule(getInputGeometry("rna.inp"))

#mode sets the radius correction method as a string: 
#   'solid sphere' constant radius, based on vander waals radius of the atom.  This is simple case, and is made available in case future calculations can be absed off of the uncorrected radius.
#   'temperature'  temperature dependent atom radius.
#   'size'  includes temperature, accounts for superposition effects due to increasing number of atoms.  This is the default mode.

PA = projectionApproximation(mol, temp=300, accuracy=1, buffr=1.4, maxCycles=1)

mol2 = getInputGeometry("rna.inp")
buffr = 1.19
accuracy = 1
conv = accuracy/100
temp=5
correction = 'solid sphere'

a = molecule.getCenterMass(mol2)
             
b = molecule.translateCenter(mol2)
c = molecule.rotateGeometry(b)
d = molecule.projectionArea(c, buffr=1.4)
e = monteCarloIntegration(d, c, buffr, conv, maxit=1000000000)






#This is a data viewing program just for checking.
hitx = []
hity =[]
for i in e[1]:
    hitx.append(i[0])
    hity.append(i[1])

carbonListx=[]   
carbonListy=[]
hydrogenListx=[]
hydrogenListy=[]
nitrogenListx=[]
nitrogenListy=[]
oxygenListx=[]
oxygenListy=[]

for i in c:
    if i[0] == 'C':
        carbonListx.append(float(i[2]))
        carbonListy.append(float(i[3]))
    if i[0] == 'H':
        hydrogenListx.append(float(i[2]))
        hydrogenListy.append(float(i[3]))
    if i[0] == 'N':
        nitrogenListx.append(float(i[2]))
        nitrogenListy.append(float(i[3]))
    if i[0] == 'O':
        oxygenListx.append(float(i[2]))
        oxygenListy.append(float(i[3]))     
        
        
        
        
        
        

from PyQt5 import QtWidgets
import pyqtgraph as pg
import sys  # We need sys so that we can pass argv to QApplication

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.graphWidget = pg.PlotWidget()
        self.setCentralWidget(self.graphWidget)
        # plot data: x, y values
        self.graphWidget.plot(hitx, hity, pen=None, symbol='o', symbolSize=1, symbolBrush='w' )
        self.graphWidget.plot(carbonListx, carbonListy, pen=None, symbol='o', symbolSize=10, symbolBrush='y' )        
        self.graphWidget.plot(oxygenListx, oxygenListy, pen=None, symbol='o', symbolSize=10, symbolBrush='r' )         
        #self.graphWidget.plot(hydrogenListx, hydrogenListy, pen=None, symbol='o', symbolSize=10, symbolBrush='b' ) 
        self.graphWidget.plot(nitrogenListx, nitrogenListy, pen=None, symbol='o', symbolSize=10, symbolBrush='g' )        
        self.graphWidget.setTitle()
        self.graphWidget.setLabel('left', text='Y-axis', units = 'Angstroms')
        self.graphWidget.setLabel('bottom', text='X-axis', units='Angstroms')

def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()