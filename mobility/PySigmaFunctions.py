# -*- coding: utf-8 -*-
"""
Prepared by Nicholas Peraino
Lumigen Instrument Center
Wayne State Univeristy

This is a separate file of functions to improve readability of MainPySigma
It is all the functions it uses for molecule transformations and list operations.
"""

import csv
import math
import random


"""Create an input geometry file as a Gaussian input using Gabedit.
This will allow for easier molecular mechanics optimization or double use of the input.
It is parsing for the sequence:
    "blank line"
    title or description line
    "blank line"
Then it starts writing to the table
"""
#Function for translating input into list for processing
def getInputGeometry(geometryInputFile):    
    geometryXY=[]
    blankCount=0
    with open(geometryInputFile) as csvfile:
        geometryStack = csv.reader(csvfile, delimiter="\t")
        for row in geometryStack:
            if blankCount <2:            
                if len(row) == 0:
                    blankCount = blankCount+1
                else:
                    pass     
            else:
                if not len(row) == 0:
                    geometryXY.append(row)
                else:
                    pass  
    del geometryXY[0]
    return geometryXY
        

"""Converts an atom label into its Exact mass
"""        
def getMass(atom):
    if "C":
        return float(12.000000)
    if "H":
        return float(1.007825)
    if "O":
        return float(15.994915)
    if "N":
        return float(14.003074)
    if "B":
        return float(11.009305)
    if "F":
        return float(18.998403)
    if "He":
        return float(4.002603)
    if "S":
        return float(31.972072)
    
"""Converts an atom label into its atomic van der waals radius in Angstrom.
"""

def getAtomRadius(atom):
    if "C":
        return float(1.6)
    if "H":
        return float(1.1)
    if "O":
        return float(1.51)
    if "N":
        return float(1.20)
    if "B":
        return float(1.21)
    if "F":
        return float(1.44)
    if "S":
        return float(1.79)
      
def randomEulerAngles():
    a1 = random.uniform(0,math.pi)
    a2 = random.uniform(0,2*math.pi)
    a3 = random.uniform(0,2*math.pi)
    return (a1, a2, a3)


def getWeight(atoms):
    weight = 0
    for row in atoms:
        weight = weight + getMass(row[0])
    return weight

def getCenterMass(atoms):
    xCM = 0
    yCM = 0
    zCM = 0
    weight = 0
    for row in atoms:
        xCM = xCM + getMass(row[0]) * float(row[2])
        yCM = yCM + getMass(row[0]) * float(row[3])
        zCM = zCM + getMass(row[0]) * float(row[4])
        weight = weight + getMass(row[0])
            
    xCM = round(xCM/weight,8)
    yCM = round(yCM/weight,8)
    zCM = round(zCM/weight,8)
        
    return [xCM, yCM, zCM]

    
def translateCenter(atoms):
    translated=[]
    CM = getCenterMass(atoms)
    for row in atoms:         
        try:
            x = round(float(row[2]),8)
        except ValueError:
            x = 0
        try:
            y = round(float(row[3]),8)
        except ValueError:
            y=0
        try:
            z = round(float(row[4]),8)
        except ValueError:
            z = 0
                    
        translatedX = x-float(CM[0])
        translatedY = y-float(CM[1])
        translatedZ = z-float(CM[2])       
        translated.append([row[0],row[1], translatedX, translatedY, translatedZ, row[5]])              
    return translated

def rotateGeometry(atoms):
    #make a tuple of random angles to rotate around.
    randomAngles = randomEulerAngles()
    a1= randomAngles[0]
    a2= randomAngles[1]
    a3= randomAngles[2]
        
    #make a matrix to transfrom the coordinates with using the random angles.      
    eumat11 = math.cos(a2)*math.cos(a3)-math.cos(a1)*math.sin(a3)*math.sin(a2)
    eumat21 = -math.sin(a2)*math.cos(a3)-math.cos(a1)*math.sin(a3)*math.cos(a2)
    eumat31 = math.sin(a1)*math.sin(a3)
    eumat12 = math.cos(a2)*math.sin(a3)+math.cos(a1)*math.cos(a3)*math.sin(a2)
    eumat22 = -math.sin(a2)*math.sin(a3)+math.cos(a1)*math.cos(a3)*math.cos(a2)
    eumat32 = -math.sin(a1)*math.cos(a3)
    eumat13 = math.sin(a1)*math.sin(a2)
    eumat23 = math.sin(a1)*math.cos(a2)
    eumat33 = math.cos(a1)
        
    # rather than calling numpy, the math was pretty straightforward
    # so I left it in this way, but for visual understanding here it is:
    
    #eumat = [[eumat11, eumat12, eumat13],
    #         [eumat21, eumat22, eumat23],
    #         [eumat31, eumat32, eumat33]]
        
    rotated=[]
    for row in atoms:
        try:
            x = round(float(row[2]),8)
        except ValueError:
            x = 0
        try:
            y = round(float(row[3]),8)
        except ValueError:
            y=0
        try:
            z = round(float(row[4]),8)    
        except ValueError:
            z = 0
                
        rotatedX = x*eumat11+y*eumat21+z*eumat31
        rotatedY = x*eumat12+y*eumat22+z*eumat32
        rotatedZ = x*eumat13+y*eumat23+z*eumat33
        rotated.append([row[0], row[1], rotatedX, rotatedY, rotatedZ, row[5]])         
    return rotated
#  def calcCollisionCrossSection(atoms, accuracy=2, maxIterations=1000, gasWeight=4):

def projectionArea(geometry, buffr):
    maxx=-10000000
    maxy=-10000000
    minx=10000000
    miny=10000000
    for row in geometry:
        maxx=max(maxx,(row[2]+(row[5])))
        maxy=max(maxy,(row[3]+(row[5])))
        minx=min(minx,(row[2]-(row[5])))
        miny=min(miny,(row[3]-(row[5])))
    
    maxx=maxx+buffr
    maxy=maxy+buffr
    minx=minx-buffr
    miny=miny-buffr
        
    return (minx,maxx,miny,maxy)
 
#do a projection area on the molecule to get the area,
#prepare the atoms first by doing a getCenterMass, then translateCenter.
        
def monteCarloIntegration(area, geometry, buffr, conv, maxit=1000000):
    minx = area[0]
    maxx = area[1]
    miny = area[2]
    maxy = area[3]
    
    boxsize=(maxx-minx)*(maxy-miny)
    
    nhit=0
    ntry=0
    hits=[]
    miss=[]
    while ntry < maxit:
        ntry = ntry+1
        xr = random.uniform(0,1)*(maxx-minx) + minx
        yr = random.uniform(0,1)*(maxy-miny) + miny
        for row in geometry:
            x = row[2]
            y = row[3]
            rad = (row[5]) + buffr    
            if (abs(x-xr)) <= rad or (abs(y-yr)) <= rad:            
                dist = math.sqrt((x-xr)**2 + (y-yr)**2)
                if dist <= rad:
                    nhit = nhit + 1
                    hits.append([xr, yr])
                    area = boxsize*nhit/ntry               
                    if ntry < 1000:
                        break
                    else:
                        err = math.sqrt(((ntry*nhit-nhit**2)/(ntry**3)))
                        if err <= conv and ntry>200:  
                            #print('Area= ',area,' Estimated error = ',err*100,' %', 'after ', ntry, ' iterations')
                            return float(area)
                        else:
                            break
                else:
                    continue
            else:
                miss.append([xr,yr])

    print('Warning, failed to converge after ', maxit, ' iterations')
    return float(area)

def atomNumbers(geometry):
    numberedGeometry=[]
    for row in geometry:
        if row[0] == 'C':
            numberedGeometry.append([6,row[1],row[2],row[3],row[4]])
        if row[0] == 'H':
            numberedGeometry.append([1,row[1],row[2],row[3],row[4]])
        if row[0] == 'O':
            numberedGeometry.append([8,row[1],row[2],row[3],row[4]])
        if row[0] == 'N':
            numberedGeometry.append([7,row[1],row[2],row[3],row[4]])
    return numberedGeometry

