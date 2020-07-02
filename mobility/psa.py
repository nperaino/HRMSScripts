# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:25:31 2020
Nick Peraino

The following are the neccesary functions for the size Projection Superpositions Approximation.
The method and parameters were published:
    
"A novel projection approximation algorithm for the fast and 
accurate computation of molecular collision cross sections.""
(I) Method
Christian Bleiholder, Thomas Wyttenbach, Michael T. Bowers
https://doi.org/10.1016/j.ijms.2011.06.014
(II) Model parameterization and definition of empirical shape factors for proteins
Christian Bleiholder, Stephanie Contreras, Thanh D. Do, Michael T. Bowers
https://doi.org/10.1016/j.ijms.2012.08.027


The supporting information in the Method describes everything well,
and the neccesary corrections to the probability based collision model,
and size factor are calculated here.
"""

from PySigmaFunctions import  threeDmontecarloIntegration, translateCenter, getAtomRadius, projectionVolume
import math
import random
import csv

def superpositionArea(area, geometry, buffr, conv, temp, maxit=200000):  
    minx = area[0]
    maxx = area[1]
    miny = area[2]
    maxy = area[3]
    
    boxsize=(maxx-minx)*(maxy-miny)
    
    nhit=0
    ntry=0
    while ntry < maxit:
        ntry = ntry+1
        xr = random.uniform(0,1)*(maxx-minx) + minx
        yr = random.uniform(0,1)*(maxy-miny) + miny
        for row in geometry:
            x = row[2]
            y = row[3]
            atom = row[0]
            rad = row[5]   
            #This bit deviates from the other methods by addition of the probability for tailing.
            hit=probHit(x,y, xr, yr, buffr, rad, atom, temp)     
            if hit == True:            
                dist = math.sqrt((x-xr)**2 + (y-yr)**2)
                if dist <= rad:
                    nhit = nhit + 1
                    area = boxsize*nhit/ntry               
                    if ntry < 1000:
                        break
                    else:
                        err = math.sqrt(((ntry*nhit-nhit**2)/(ntry**3)))
                        if err <= conv and ntry>200: 
                            return float(area)
                        else:
                            break
                else:
                    continue
            else:
                continue
    print('Warning, failed to converge after ', maxit, ' iterations')
    return float(area)


"""
This accounts for the tailing collision radius as a probablilty density.
there is a slight increase in size with increasing atom number - the superposition part.
This is accounted by addition of the probability that if the selected montecarlo radius of collision
is less than or equal to the hard sphere radius probability of collision is guaranteed or 1.
Otherwise, there is a reducing probability of collision by a function of the molecule size, not an
instantaneous zero probability.
it is equation 5 from the first paper using parameters file from the second paper.   

    ra      εb      k1      l1      k2      l2      alpha
H   2.67146 0.01059 0.22267 0.94441 0.21167 0.74496 0.64682
C   3.49512 0.01765 0.15835 0.62905 0.15179 0.43198 0.85604
N   3.49512 0.01765 0.13021 0.79003 0.12483 0.34075 0.95252
O   3.49512 0.01765 0.17505 0.62290 0.19355 0.52855 1.10393
Na  3.97000 0.00154 0.12817 0.57818 0.28347 0.40100 0.65908
K   4.47000 0.00154 0.04470 0.69859 0.09170 0.55420 1.10628

a Values given in Å.
b Values given in kcal/mol.

"""

def probHit(x, y, xr, yr, buffr, rad, atom, temp):
    def tailingFunction(x, y, rad, temp, r, e, k1, l1, k2, l2, alpha):
        #This is the probablility tailing function. I broke up the exponents into constants
        #for readability.  The coordinate can be the x or y coordinate for the projection.
        constantX1 = -(abs(x)-rad)*k1*temp**l1
        constantX2 = -(abs(x)-rad)*k2*temp**l2
        constantY1 = -(abs(y)-rad)*k1*temp**l1
        constantY2 = -(abs(y)-rad)*k2*temp**l2
        
        #probability of hit in x axis
        ProbTx = alpha*(math.e**constantX1)+(alpha-1)*(math.e**constantX2)
        #probability of hit in y axis
        ProbTy = alpha*(math.e**constantY1)+(alpha-1)*(math.e**constantY2)
      
        #non mutually exclusive probability of hit in X or Y axis
        return ProbTx+ProbTy-(ProbTx*ProbTy)
    
    
    #first read the parameters file to get the coefficients for the probabibility function.
    probabilityParameters=[]
    
    with open('psa.par') as csvfile:
        probParameters = csv.reader(csvfile, delimiter=',')
        for row in probParameters:
            probabilityParameters.append(row)      
    
    for row in probabilityParameters:
        #I am defining r as none here so that if it an atom does not have parameters,
        #we can just assign it as an Oxygen atom.  I am expecting that for most unknowns,
        #It is probably Sulfur or some other smallish, diffuse atom... we have to put something or skip it I guess.
        r = None
        if str(atom) == str(row[0]):
                r = (float(row[1]))
                e = (float(row[2]))
                k1 = (float(row[3]))
                l1 = (float(row[4]))
                k2 = (float(row[5]))
                l2 = (float(row[6]))
                alpha = (float(row[7]))  
        # if r is still not assigned, then 
    if r == None:
        r = 3.49512
        e = 0.01765
        k1 = 0.17505
        l1 = 0.62290
        k2 = 0.19355
        l2 = 0.52855
        alpha = 1.10393              
    if (abs(x-xr)) <= rad or (abs(y-yr)) <= rad:            
        dist = math.sqrt((x-xr)**2 + (y-yr)**2)
        if dist <= rad:
            return True
    else:
        probabilityHit = tailingFunction(x, y, rad, temp, r, e, k1, l1, k2, l2, alpha)
        return random.random() < probabilityHit
        
"""
The shape factor is defined by the ratio of the convex envelope of the projections divided by the total
molecular surface area.  A convexhull is applied to the list of hits from the montecarlo to get the convex shape,
and the molecular surf
    
"""
def shapeFactor(geometry, acc):
    conv = acc/100
    vanderwaalsGeometry=[]
    for row in geometry:
        vanderwaalsGeometry.append([row[0], row[1], row[2], row[3], row[4], getAtomRadius(row[0])])   
    centerGeometry = translateCenter(vanderwaalsGeometry)
    projectionCube = projectionVolume(centerGeometry, buffr=0)
    integration = threeDmontecarloIntegration(projectionCube, vanderwaalsGeometry, conv, maxit=1000000)
    
    Ace = integration[1]
    Amol = integration[2]
    shapeFactor = Ace/Amol
    return shapeFactor
