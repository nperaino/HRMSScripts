# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:56:23 2020

@author: Nick

This is the main function for the projection approximation
"""
import PySigmaFunctions
from potential import *
from PySigmaFunctions import *
import math
from statistics import mean

def projectionApproximation(geometry, temp=300, accuracy=1, buffr=1.09, maxCycles=20):
    data=[]
    sumarcor=0
    suma=0
    mina=0
    maxa=0
    conv = accuracy/100
    for i in range(1,maxCycles):
                #Translate center mass coordinates to position "0,0,0"
                #So that we can radnomly rotate around the center mass
        centerGeometry = molecule.translateCenter(geometry.atoms)
                #Select a random rotation around this center mass.
        rotationGeometry = molecule.rotateGeometry(centerGeometry)
                #Calculate the miaximum projection Area in the XY plane
                #pramaterizes the box that includes the surface of the molecule.
        projectionBox = molecule.projectionArea(rotationGeometry, buffr)
                #Monte Carlo integration of the space to find where the molecule is.
                #Collision is determine by measurign the distance between the buffer gas
                #and any of the atoms in the structure.  If a distance is found to be less
                #than the two radii of the gas and atom, then there is a collision.
        area = monteCarloIntegration(projectionBox, rotationGeometry, buffr, conv)
                #Temperature correction.
        
        tStar = temp/(1.7179E4/(area[0]/math.pi)**2)
        omStar = phs4(tStar)     
        reducedArea = phs4(tStar)
        data.append(reducedArea)
        areacor= omStar*area[0]
        sumarcor=sumarcor+areacor
        aareacor=sumarcor/i
        suma=suma+area[0]
        aarea=suma/i

#     --- do at least 10 different orientations, then calculate
#     statistical errors (sigma) assuming a gaussian distribution
#     and test for convergence.
        
        chi=0.0       
        if i > 10:
            for k in range(0,i):
                chi=chi+(data[k]-aarea)**2
            deviation = math.sqrt(chi/(i-1))/math.sqrt(i)

        mina=min(mina,area[0])
        maxa=max(maxa,area[0])
#
#     --------- Test for convergence :
#
        if (i > 30):
#           if ((100.*deviation/aarea).lt.(acc)) then
            if ((deviation/aarea) < (conv)):
                mob=(1.0/math.sqrt(temp*weight*he/(weight+he)))*1.85e4/aarea
                dmob=(1.0/math.sqrt(temp*weight*he/(weight+he)))*1.85e4/(aarea+deviation)
                
                
#           Do some first order temperature correction, 1.7179D+4 includes 0.205
#           for He polarizability. !!Change for other drift gases!!!
#           Note: it makes little to no difference if the temperature correction
#                   is done for each iteration or just at the end. Since a
#                    correction for each iteration seems better, it is given
#                    back from this subroutine
            ave = mean(data)
            tStar = temp/(1.7179E4/(ave/math.pi)**2)

            omStar = phs4(tStar)
            

#
#            write(iout,*) ' Iterations            =',i
#            write(iout,*) ' Note: the deviations below have nothing'
#            write(iout,*) ' to do with the uncertanty of the average'
#            write(iout,*) ' Maximum Cross-section =',maxa
#            write(iout,*) ' Minimum Cross-section =',mina
#            write(iout,*) ' Average Deviation     =',adev
#            write(iout,*) ' Standard Deviation    =',sdev
#            write(iout,*) ' Variance              =',var
#            write(iout,*) ' Skewness              =',SKEW
#            write(iout,*)
#            write(iout,*) ' Average Cross-section =',
#     1           ave,'+/-',deviation
#            write(iout,*) ' Mobility              =',mob,'+/-',mob-dmob
#            write(iout,*) ' Temperature           =',degK
#            write(iout,*) ' The following values make only sense'
#            write(iout,*) ' for a hard sphere calculation (no -s)!!'
#            write(iout,*) ' Reduced Temperature    =',tstar
#            write(iout,*) ' Omegastar              =',omstar
#            write(iout,*) ' 1st ord. temp cor. Cross-sec. =',areacor
#            write(iout,*) ' 2nd ord. temp cor. Cross-sec. =',aareacor
#
            cross=ave
            crosst=aareacor
            print('\n deveiation = ' + str(ave), 
                  '\n 1st ord. temp cor. Cross-sec = ' + str(aareacor), 
                  '\n 2nd ord. temp cor. Cross-sec = ' + str(areacor))
            return