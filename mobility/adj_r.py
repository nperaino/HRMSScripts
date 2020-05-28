# -*- coding: utf-8 -*-
"""
Created on Wed May 20 12:10:29 2020

@author: Nick

# ***********************************************************************
#     read the 6-n params file. Format: atnum, rep, rmin,emin
# here is an example file for MM3 params
# 0 0.205 1.53 0.026
# 1 12.0 1.62 0.02
# 6 12.0 2.04 0.027
# 8 12.0 1.82 0.059
# 11 12.0 2.70 0.057
#
# references: MM3 params: Allinger et al. Theochem, Vol312, 69 (1994)
# and Halgren, JACS, 114, 7872 (1992)
#
# ***********************************************************************

temperature corrected radii was set by using -s flag in Sigma.

"""
from PySigmaFunctions import *
import csv
import math
import graph
from potential import *


#def adj_r(natoms,maxatoms,atomr,buffr,atnum, degK,sixnfile,rep,emin,rmin,cfour,rnew,enew,gam)
def adj_r(geometry, temp):
    
    #reads through the LJ parameters file and sets up
    #a list with the parameters for each molecule in the 
    #geometry.
    
    sixnfile=[]
    with open('mm3.pot') as csvfile:
        LJParameters = csv.reader(csvfile, delimiter=" ")
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
            
    cfour=[]
    rnew=[]
    enew=[]
    gam=[]
    atomr=[]
    for i in range(0,len(geometry)):
        radd=rminbuf
        #rmin[i]=rmin[i]+radd

        rmin[i]=2.0*rmin[i]
        radd=2.0*rminbuf
    # Combination rules from the halgren paper "HHG" for radii and (14) for emin
        reg0=(rmin[i]**3 + radd**3)/(rmin[i]**2 + radd**2)
        rmin[i]=reg0
    # emin[i]=math.sqrt(eminbuf*emin[i])
        emin[i]=(4.0*eminbuf*emin[i])/(math.sqrt(eminbuf)+math.sqrt(emin[i]))**2
# Units: kcal/mol for potential, Angstroms for radii
#     calculate constant C4 of r^-4 term
        cfour.append(0.5*polbuf*(3.3205E2/(len(geometry)**2)))
#     find new position of potential minimum [dV(n_6_4)/dr = 0]
#       write (*,*) 'Alternate zbrent entry, shouldnt'
        rnew.append(zbrent(rep[i],emin[i],rmin[i],cfour[i],float(1.0),float(10.0),float(0.0001)))
#
#     get new well depth
        enew.append(-vnsixfour(rep[i],emin[i],rmin[i],rnew[i],cfour[i]))
        tstar=temp/(enew[i]*5.032E2)
#     calculate gamma
#         gam(i)=1-((3.3205D+2/(natoms**2))*polbuf)/(6.0*rnew(i)**4)
#         This was Gert's formula. It is not correct from my point of view.
#         Thomas 12/12/95
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
    print(newgeometry)
    return newgeometry
            

 
    
#this is from Numerical Recipes "Brent's Method" for root finding.
def zbrent(rep,emin,rmin,cfour,x1,x2,tol):
    itmax=100
    eps=3.e-8
    a=x1
    b=x2
    
    fa = dvdr(a,rep,emin,rmin,cfour).real
    fb = dvdr(b,rep,emin,rmin,cfour).real

    if(fb*fa > 0.):
        print('zbrent : root must be bracketed')
    fc=fb
    for i in range(0,itmax):
        if ((fb*fc).real > 0):
            c=a
            fc=fa
            d=b-a
            e=d
    
        if (abs(fc) < abs(fb)):
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        tol1=2.*eps*abs(b)+0.5*tol
        xm=.5*(c-b)
        if (abs(xm) <= tol1 or fb == 0):
            zbrent=b
            return zbrent
         
        if (abs(e) >= tol1 and abs(fa) > abs(fb)):
            s=fb/fa
            if(a == c):
                p=2.*xm*s
                q=1.-s
            else:
                q=fa/fc
                r=fb/fc
                p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                q=(q-1.)*(r-1.)*(s-1.)

            if (p.real > 0.):
                q=-q
            p=abs(p)
            
            if((2.*p.real) < (min(3.*xm.real*q.real-abs(tol1.real*q.real),abs(e.real*q.real)))):
               e=d
               d=p/q
            else:
               d=xm
               e=d
        else:
            d=xm
            e=d
        a=b
        fa=fb
        if(abs(d) > tol1):
            b=b+d
        else:
            b=b+math.copysign(tol1.real,xm.real)
        fb=dvdr(b,rep,emin,rmin,cfour)
    zbrent=b
    return zbrent
    
# The first derivative of a n-6 potential with a r-4 tagged on
def dvdr(x,rep,emin,rmin,cfour):
    gam=1.0
    pre=rep*emin/(rep*(3.0+gam)-(12.0*(1+gam)))
    vrep=-pre*(12.0*(1+gam)/x *(rmin/x)**rep)
    vsix=pre*4.0*gam*6.0/x *(rmin/x)**6
    vfour=4.0*cfour/(x**5)
    dvdr= vrep+vsix+vfour
    return dvdr 
          
# The n -6 -4 potential
def vnsixfour(rep,emin,rmin,r,cfour):
    gam=1.0
    pre=rep*emin/(rep*(3.0+gam)-(12.0*(1+gam)))
    vrep=pre*12.0/rep *(1+gam)*(rmin/r)**12
    vsix=-pre*4.0*gam*(rmin/r)**6
    vfour=-cfour/r**4
    vnsixfour=vrep+vsix+vfour
    return vnsixfour

def QRED(t,g,r):
    Q3 = P16_6_4(t,g)
    Q2 = P12_6_4(t,g)
    Q1 = P8_6_4(t,g)   
    Q = [Q1, Q2, Q3]
    X = [8., 12., 16.]
    QN = interpolate.pchip_interpolate(X, Q, r)   
    return QN
 

        
        
temp = 500      
        
mol = getInputGeometry("C60.inp")       
m = adj_r(mol, temp)
n = atomNumbers(mol)