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
import tkinter as tk
from tkinter.filedialog import askopenfilename
from projectionApproximation import *


#Ask to open the input geometry commented out to save development time
#root = tk.Tk()
#root.withdraw()
#geometryInputFile = askopenfilename()

geometry = getInputGeometry("c60.inp")

#mode sets the radius correction method as a string: 
#   There is no real useful reason to use any mode other than 'temperature and size'.
#   'solid sphere' constant radius, based on vander waals radius of the atom.  This is simple case, and is made available in case future calculations can be absed off of the uncorrected radius.
#   'temperature'  temperature dependent atom radius.
#   'temperature and size'  includes temperature, accounts for superposition effects due to increasing number of atoms.  This is the default mode.
#   Temp cannot be 0, MaxCycles seems to do best over 30 but certainly depends on input.  If you dont give it enough cycles it doesn't converge.


#PA = projectionApproximation(geometry, temp=100, accuracy=1, buffr=1.4, maxCycles=100, mode='temperature and size')
#print(PA)


Temp= []
CCSsolid =[]
from tqdm import tqdm 
for i in tqdm(range(50,605,50)):
    CCSsolid.append(projectionApproximation(geometry, temp=i, accuracy=1, buffr=0.109, maxCycles=50, mode='temperature'))
    Temp.append(i)