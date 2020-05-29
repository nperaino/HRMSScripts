# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:44:40 2020

@author: Nick
"""
from PySigmaFunctions import getAtomRadius, translateCenter, rotateGeometry, projectionArea, monteCarloIntegration, atomNumbers
import csv
import math
from scipy.optimize import minimize_scalar
from potential import *



def adj_heion(geometry, temp): 

# ****************************************************************************
#     read the 6-n params file. Format: atnum, rep, rmin,emin
#     rmin is location of He-X potential minimum (X is atom in polyatomic ion)
#     emin is depth of He-X potential minimum
#*****************************************************************************
    sixnfile=[]
    with open('adj.pot') as csvfile:
        LJParameters = csv.reader(csvfile, delimiter=' ')
        for row in LJParameters:
            sixnfile.append(row)      
    numberedGeometry = atomNumbers(geometry)
    #sets the LJ parameters for the geometry.
    rep=[]
    rmin=[]
    emin=[]
    for atomrow in numberedGeometry:
        for row in sixnfile:
            if str(atomrow[0]) == str(row[0]):
                #row[1] is rep, row[2] is rmin, row[3] is emin
                rep.append(float(row[1]))
                emin.append(float(row[3]))
                
                # ******************************
                # Size scaling takes place here
                # ******************************
                natoms = len(geometry)
                rvar1 =   .86882
                rvar2 =  -.99427
                rvar3 =   .99913
                if natoms > 222:
                    scalfr = 1.1
                else:
                    scalfr = rvar1 + rvar2*natoms + natoms**rvar3
                scaledrmin = float(row[2])*scalfr
                rmin.append(scaledrmin)              
                continue          
    #sets the LJ parameters for the buffer gas.
    #called as atom 0 in the LJ input parameters.
    #declared the variables as 0 just to initialize them.
    rminbuf=0
    polbuf=0
    eminbuf=0   
    for row in sixnfile:
        if str(row[0]) == str(0):
            polbuf=float(row[1])
            rminbuf=float(row[2])
            eminbuf=float(row[3])
            
# ******************************
#This is all very redundant since it is just a copy from adj_r.        

#iterate over the geometry and correct each minimum radius by finding the radius that gives the minimum
#energy of the energy well calculated from the Lennard Jones parameters defining the n-6-4 potential.
            
    cfour=[]
    rnew=[]
    enew=[]
    gam=[]
    atomr=[]

    for i in range(0,len(geometry)):
# Units: kcal/mol for potential, Angstroms for radii
#     calculate constant C4 of r^-4 term
        cfour.append(0.5*polbuf*(3.3205E2/(len(geometry)**2)))

#     find new position of potential minimum [dV(n_6_4)/dr = 0]
#     it should reasonably be between 1 and 10 angstrom or something is wrong.
#     pass the equation coefficients to the function using current iterations rep, emin, rmin, and cfour values
#     so the bounded Brent optimization function can accept the function.
        #minimum = minimize_scalar(dvdr,args=(rep[i], emin[i], rmin[i], cfour[i]), bounds=(1, 10), tol=3E-8, method='bounded')        
        minimum = minimize_scalar(vnsixfour,args=(rep[i], emin[i], rmin[i], cfour[i]),bracket=(1,10), options={'xtol':3E-8}, method='golden')
        rnew.append(minimum.x)

#     get new well depth
        enew.append(-1*vnsixfour(rnew[i],rep[i],emin[i],rmin[i],cfour[i]))
        tstar=temp/(enew[i]*5.032E2)

#     calculate gamma
        xnumratr = 4 * cfour[i] * (6-rep[i])
        dnomnatr = 3 * rep[i] * enew[i] * rnew[i]**4 + cfour[i] * (rep[i]-12)
        gam.append(1 + xnumratr/dnomnatr)

#     calculate omegastar
        qN = QRED(tstar,gam[i],rep[i])

#     calculate the new radius
        atomr.append(rnew[i]*math.sqrt(qN))
            
#        print('old:', numberedGeometry[i],rmin[i],emin[i])
#        print('new:', numberedGeometry[i],rnew[i],enew[i],gam[i])
#D       write(*,*) 'omega', qN,atomr(i),tstar 
        buffr=0.0
    newgeometry=[]
    for row, radius in zip(geometry, atomr):
        newgeometry.append([row[0], row[1], row[2], row[3], row[4], radius])
    return newgeometry

    
# The first derivative of a n-6 potential with a r-4 tagged on
def dvdr(x,rep,emin,rmin,cfour):
    gam=1.0
    pre=rep*emin/(rep*(3.0+gam)-(12.0*(1+gam)))
    vrep=-pre*(12.0*(1+gam)/x *(rmin/x)**rep)
    vsix=pre*4.0*gam*6.0/x *(rmin/x)**6
    vfour=4.0*cfour/(x**5)
    dvdr= vrep+vsix+vfour
    return dvdr.real
          
# The n -6 -4 potential
def vnsixfour(x,rep,emin,rmin,cfour):
    gam=1.0
    pre=rep*emin/(rep*(3.0+gam)-(12.0*(1+gam)))
    vrep=pre*12.0/rep *(1+gam)*(rmin/x)**12
    vsix=-pre*4.0*gam*(rmin/x)**6
    vfour=-cfour/x**4
    vnsixfour=vrep+vsix+vfour
    return vnsixfour.real

def QRED(t,g,r):
    Q3 = P16_6_4(t,g)
    Q2 = P12_6_4(t,g)
    Q1 = P8_6_4(t,g)   
    Q = [round(Q1,5), round(Q2,5), round(Q3, 5)]
    X = [8., 12., 16.]
    QN = interpolate.pchip_interpolate(X, Q, r)   
    return QN
 

        
#some settings for testing       
#temp = 500            
#mol = getInputGeometry("C60.inp")       
#m = adj_r(mol, temp)
#n = atomNumbers(mol)