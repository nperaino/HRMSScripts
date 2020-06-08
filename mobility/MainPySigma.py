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

from PySigmaFunctions import getInputGeometry
import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename 
from projectionApproximation import projectionApproximation
from tqdm import tqdm 
from PyQt5 import QtWidgets
import pyqtgraph as pg
import sys
import csv


print('This script will calculate the temperature and size dependent collisional cross setion (CCS)')
print('from an input file generated using Gabedit as a Gaussian input file.  This way, Gabedit or Gaussian')
print('can be used to optimize the gas phase geometry before calculating the CCS.  It parases by line counting')  
print('')
print('')
print('')
input('Please press Enter to select an input file')
root = tk.Tk()
root.withdraw()
geometryInputFile = askopenfilename()
geometry = getInputGeometry(geometryInputFile)
print('Enter 1 for single temperature run')
print('Enter 2 for curve')
tempType = input('Enter Selection:')
print(tempType)
print('')
print('')
print('')
print('Please Enter buffer gas radius.  Use 0.109 if you are not sure.')
bufferRadius = float(input('Enter buffer gas radius:'))
print('')
print('')
print('')
print('Please enter an accuracy.  1 is ok for smaller molecules less than 100 atoms')
print('but larger molecules will take longer so raise up to reduce time')
acc = float(input('Accuracy: '))
print('')
print('')
print('')
print('Please enter the maximum number of rotations.  Very symmetric molecules need fewer')
print('If you are unsure, run a type 2 temperature curve and see if it is very inconsitent.')
print('50 seems to be ok for small molecules.')
maxC = int(input('Enter maximum cycles:'))

if tempType == '1':
    print('Please enter the temperature in Kelvin.')
    tempK = float(input('Temperature: '))
    PA = projectionApproximation(geometry, temp = tempK, accuracy=acc, buffr=bufferRadius , maxCycles=maxC, mode='temperature and size')
    print('The projection area is: ',PA)
    input('press enter to quit')
elif tempType == '2':
    Temp= []
    CCS =[]
    for i in tqdm(range(20,605,20)):
        CCS.append(projectionApproximation(geometry, temp=i, accuracy=acc, buffr=bufferRadius , maxCycles=maxC, mode='temperature and size'))
        Temp.append(i)
    print('Please select a save file location:')

    saveData = asksaveasfilename()      
    CCSdata =[]
    for i,j in zip(Temp, CCS):
        CCSdata.append([i,j])
    with open(saveData+'.csv', 'w', newline='') as csvfile: 
        csvwriter = csv.writer(csvfile)
        headers=['Temperature','CCS']
        # writing the headers 
        csvwriter.writerow(headers) 
        # writing the data rows
        csvwriter.writerows(CCSdata)
    print('A graph will display the results.')
    title = input('please input a title for the graph:')
    class MainWindow(QtWidgets.QMainWindow):

        def __init__(self, *args, **kwargs):
            super(MainWindow, self).__init__(*args, **kwargs)
            self.graphWidget = pg.PlotWidget()
            self.setCentralWidget(self.graphWidget)
            self.graphWidget.setTitle("<span style=\"color:black;font-size:10px\">"+title+"  </span>")
            self.graphWidget.setLabel('left', 'CCS (A<sup>2</sup>)', color='black', size=30)
            self.graphWidget.setLabel('bottom', 'Temperature (K)', color='black', size=30)
            self.graphWidget.setBackground('w')
            pen = pg.mkPen(color=(0, 0, 0))
            self.graphWidget.showGrid(x=True, y=True)
            self.graphWidget.plot(Temp, CCS, pen=pen)


    def main():
        app = QtWidgets.QApplication(sys.argv)
        main = MainWindow()
        main.show()
        sys.exit(app.exec_())


    if __name__ == '__main__':
        main()

