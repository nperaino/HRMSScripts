# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:56:23 2020

@author: Nick

This is the main function for the projection approximation
"""

from potential import phs4
from PySigmaFunctions import translateCenter, rotateGeometry, getAtomRadius, monteCarloIntegration, projectionArea
import math
from statistics import mean
import adj_r
import adj_heion
from psa import superpositionArea
from tqdm import tqdm

def projectionApproximation(geometry, temp=300, accuracy=1, buffr=0.109, maxCycles=100, mode = 'temperature and size'):
    data=[]
    datacor=[]
    mina=0
    maxa=0
    conv = accuracy/100
    #he = 4
    #weight = getWeight(geometry)    

# This sets the arom radii based on the correction methods.
    newgeometry=[]
    if mode == 'solid sphere':
        for row in geometry:
            newgeometry.append([row[0], row[1], row[2], row[3], row[4], getAtomRadius(row[0])])
    if mode == 'temperature':
        newgeometry = adj_r.adj_r(geometry, temp)      
    if mode == 'temperature and size':
        newgeometry = adj_heion.adj_heion(geometry, temp)

    for i in range(1,maxCycles):
                #Translate center mass coordinates to position "0,0,0"
                #So that we can radnomly rotate around the center mass
        centerGeometry = translateCenter(newgeometry)
                #Select a random rotation around this center mass.
        rotationGeometry = rotateGeometry(centerGeometry)
                #Calculate the miaximum projection Area in the XY plane
                #pramaterizes the box that includes the surface of the molecule.
        projectionBox = projectionArea(rotationGeometry, buffr)
                #Monte Carlo integration of the space to find where the molecule is.
                #Collision is determine by measurign the distance between the buffer gas
                #and any of the atoms in the structure.  If a distance is found to be less
                #than the two radii of the gas and atom, then there is a collision.
        area = monteCarloIntegration(projectionBox, rotationGeometry, buffr, conv)
        
        #First Order Temperature Correction
        tStar = temp/(1.7179E4/(area/math.pi)**2)
        omStar = phs4(tStar)     
        areacor= omStar*area
        
        data.append(area) 
        datacor.append(areacor)
 
        averageArea = mean(data)
        averageCorrectedArea = mean(datacor)
        
#     --- do at least 10 different orientations, then calculate
#     statistical errors (sigma) assuming a gaussian distribution
#     and test for convergence.
        
        chi=0.0       
        if i > 10:
            for k in range(0,i):
                chi=chi+(data[k]-area)**2
            deviation = math.sqrt(chi/(i-1))/math.sqrt(i)
        mina=min(mina,area)
        maxa=max(maxa,area)
        
#
#     --------- Test for convergence :
#

        if (i > 30):
            if ((deviation/averageArea) < (conv)):
                #mob=(1.0/math.sqrt(temp*weight*he/(weight+he)))*1.85e4/averageArea
                #dmob=(1.0/math.sqrt(temp*weight*he/(weight+he)))*1.85e4/(averageArea+deviation)
                
                
#           Do some first order temperature correction, 1.7179D+4 includes 0.205
#           for He polarizability. !!Change for other drift gases!!!
#           Note: it makes little to no difference if the temperature correction
#                   is done for each iteration or just at the end. Since a
#                    correction for each iteration seems better, it is given
#                    back from this subroutine
                tStar = temp/ (1.7179E4/(averageArea/math.pi) **2)
                omStar = phs4(tStar)
                averageCorrectedArea = omStar*averageArea
            
                if mode == 'solid sphere':
                    return averageCorrectedArea
                else:
                    return averageArea
    if mode == 'solid sphere':
        return averageCorrectedArea
    else:
        return averageArea

"""
This calculates the geometry for molecules larger than 1000 atoms.  It takes much longer,
but still works for small molecules but it can be a waste of time.
"""

def PSA(geometry, temp=300, accuracy=1, buffr=1.09, maxCycles=100):
    data=[]
    mina=0
    maxa=0
    conv = accuracy/100
    newgeometry = adj_heion.adj_heion(geometry, temp)
    for i in range(1,maxCycles):
                #Translate center mass coordinates to position "0,0,0"
                #So that we can radnomly rotate around the center mass
        centerGeometry = translateCenter(newgeometry)
                #Select a random rotation around this center mass.
        rotationGeometry = rotateGeometry(centerGeometry)
                #Calculate the miaximum projection Area in the XY plane
                #pramaterizes the box that includes the surface of the molecule.
        projectionBox = projectionArea(rotationGeometry, buffr)
                #Monte Carlo integration of the space to find where the molecule is.
                #Collision is determine by measuring the distance between the buffer gas
                #and any of the atoms in the structure.  If a distance is found to be less
                #than the two radii of the gas and atom, then there is a collision.
        
        area = superpositionArea(projectionBox, rotationGeometry, buffr, conv, temp, maxit=200000)
        data.append(area) 
        averageArea = mean(data)
        
#     --- do at least 10 different orientations, then calculate
#     statistical errors (sigma) assuming a gaussian distribution
#     and test for convergence.
        
        chi=0.0       
        if i > 10:
            for k in range(0,i):
                chi=chi+(data[k]-area)**2
            deviation = math.sqrt(chi/(i-1))/math.sqrt(i)
        mina=min(mina,area)
        maxa=max(maxa,area)
        
#
#     --------- Test for convergence :
#

        if (i > 30):
            if ((deviation/averageArea) < (conv)):          
                return averageArea
    else:
        print('Reached maximum cycles.')
        return averageArea