# -*- coding: utf-8 -*-
"""
Created on Wed May 20 12:10:29 2020

@author: Nick

# ***********************************************************************
# read the 6-n params file. Format: atnum, rep, rmin,emin
# here is an example file for MM3 params
# 0 0.205 1.53 0.026
# 1 12.0 1.62 0.02
# 6 12.0 2.04 0.027
# 8 12.0 1.82 0.059
# 11 12.0 2.70 0.057
#
# references: MM3 params: Allinger et al. Theochem, Vol312, 69 (1994)
# and 
# Halgren, JACS, 114, 7827 (1992) https://doi.org/10.1021/ja00046a032
#
# ***********************************************************************

temperature corrected radii was set by using -s flag in Sigma.

"""
from PySigmaFunctions import atomNumbers
import csv
import math
from scipy.optimize import minimize_scalar
from potential import P16_6_4, P12_6_4, P8_6_4
from scipy import interpolate



#def adj_r(natoms,maxatoms,atomr,buffr,atnum, degK,sixnfile,rep,emin,rmin,cfour,rnew,enew,gam)
def adj_r(geometry, temp):
    
    #reads through the LJ parameters file and sets up
    #a list with the parameters for each molecule in the 
    #geometry.   
    sixnfile=[]
    rep=[]
    rmin=[]
    emin=[]    
    with open('adj.pot') as csvfile:
        LJParameters = csv.reader(csvfile, delimiter=" ")
        for row in LJParameters:
            sixnfile.append(row)      
    numberedGeometry = atomNumbers(geometry)
    
    #sets the LJ parameters for the geometry.
    for atomrow in numberedGeometry:
        for row in sixnfile:
            if str(atomrow[0]) == str(row[0]):
                #row[1] is rep, row[2] is rmin, row[3] is emin
                rep.append(float(row[1]))
                rmin.append(float(row[2]))
                emin.append(float(row[3]))
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

#iterate over the geometry and correct each minimum radius by finding the radius that gives the minimum
#energy of the energy well calculated from the Lennard Jones parameters defining the n-6-4 potential.
            
    cfour=[]
    rnew=[]
    enew=[]
    gam=[]
    atomr=[]

    for i in range(0,len(rmin)):
        radd=rminbuf
        rmin[i]=2.0*rmin[i]
        radd=2.0*rminbuf
        
# Combination rules from the halgren paper "HHG" for radii and (14) for emin
        reg0=(rmin[i]**3 + radd**3)/(rmin[i]**2 + radd**2)
        rmin[i]=reg0
        emin[i]=(4.0*eminbuf*emin[i])/(math.sqrt(eminbuf)+math.sqrt(emin[i]))**2

# Units: kcal/mol for potential, Angstroms for radii
#     calculate constant C4 of r^-4 term
        cfour.append(0.5*polbuf*(3.3205E2/(len(rmin)**2)))

#     find new position of potential minimum [dV(n_6_4)/dr = 0]
#     it should reasonably be between 1 and 10 angstrom or something is wrong.
#     pass the equation coefficients to the function using current iterations rep, emin, rmin, and cfour values
#     so the bounded Brent optimization function can accept the function.
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
    newgeometry=[]
    for row, radius in zip(geometry, atomr):
        newgeometry.append([row[0], row[1], row[2], row[3], row[4], radius])
    return newgeometry
         
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
