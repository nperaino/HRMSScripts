#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#prepared by Nicholas Peraino, Lumigen Instrument Center, Wayne State University
"""
Created on Mon Sep 23 12:12:19 2019

@author: Nick Peraino
"""


from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
import itertools
import re
import csv
import tkinter as tk

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Fragmentor 1.2')
print('')
print('This script calculates the Mass Spec fragmentations for an expected compound based on RDkit')
print('implementation of substructure matching based on SMARTS.')
print('Input files must be entered in as SMILES')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('')
print('')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')

#the r must be included to avoid unicode errors wherein the SMILES code starts a unicodeescape.
#In the future a line of code must add this r to files which are read in.

mcy = input("Enter a SMILES value for the structure:")
# example smiles to use:   mcy-DASP3-LR = r'O=C([C@@H](CC(N[C@@H](CCC/N=C(N)\N)C(N[C@@H](/C=C/C(C)=C/[C@H](C)[C@@H](OC)CC1=CC=CC=C1)[C@H](C)C(N[C@@H](C(O)=O)CC2)=O)=O)=O)NC([C@H](CC(C)C)NC([C@@H](C)NC(C(N(C)C2=O)=C)=O)=O)=O)O'
m = Chem.MolFromSmiles(mcy)
print('')
print('')
print('Please select an output folder for the data.')
input('Press Enter to select a file')
from tkinter.filedialog import asksaveasfilename
root = tk.Tk()
root.withdraw()
fragments = asksaveasfilename()
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('')
print('')
print('SMARTS manual is available at Daylight Chemical Information Systems https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html')
print('')
print('')
print('Examples for peptide analysis:')
print('Methoxy = $([CHX4][O][CH3])')
print('benzylic = $([CH2X4][c])')
print('methyl-methoxy-ethane = $([CH3][CHX4][CHX4][O][CH3])')
print('decarboxylation  = $([CX3](=[OX1])O)')
print('peptide fragmentation = $([C,N][C,N][C]=[O])')
print('Recursive inclusion of all above: $([C,N][C,N][C]=[O]),$([CHX4][O][CH3]),$([CH2X4][c]),$([CH3][CHX4][CHX4][O][CH3]),$([CX3](=[OX1])O)')
print('Good peptide set: $([C,N]-[C,N]),$([C,N]-[C,N]-[C,N](-[C,N]-[C,N])=[O]),$([C,N]-[C,N]-[C,N](-[C,N]-[C,N])=[O]),$([C,N]-[C,N]-[C,N]-[C]=[O]),$([C,N]-[C,N]-[C]=[O]),$([C,N]-[C]=[O]),$([CHX4][O][CH3]),$([CH2X4][c]),$([CH3][CHX4][CHX4][O][CH3]),$([CX3](=[OX1])O)')
print('Separate the codes which are included in $(___) by commas.')  
print('')
print('')  
fr = input("Enter a SMARTS code for desired fragment.  Type [peptides] without the brackets, case sensitive, for shortcut to peptide set.")
if fr=='peptides':
    fr='$([SX2]),$([S]-[S]),$([S]-[C]),$([C]-[S]),$([C,N]-[C,N]),$([C,N]-[C,N]-[C,N](-[C,N]-[C,N])=[O]),$([C,N]-[C,N]-[C,N](-[C,N]-[C,N])=[O]),$([C,N]-[C,N]-[C,N]-[C]=[O]),$([C,N]-[C,N]-[C]=[O]),$([C,N]-[C]=[O]),$([CHX4][O][CH3]),$([CH2X4][c]),$([CH3][CHX4][CHX4][O][CH3]),$([CX3](=[OX1])O)'


#--------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------

#This sections determines the fragments resulting from breaking a carbonyl bond.
#--------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------

#amide $([N][CX3]=[O]),


#Finds Atom pairs exhibiting the SMARTS code for the bond. Converts from tuple to list, then removes duplicates.
bis = m.GetSubstructMatches(Chem.MolFromSmarts('['+fr+']'))
bonds = [item for t in bis for item in t]
#bonds = list(OrderedDict.fromkeys(bonds))

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('')
print('')
print('The compound has '+str(len(bonds))+' bonds matching the entered SMARTS code.')
print('How many bonds should be broken in any combination?')
print('Note: breaking more bonds will add significantly more time to the calculation and is going to produce a lot of redundant fragments.  3 or 4 is usually fine')
      
bnds = input('Input the number of bonds to break.')

if int(bnds)<=0:
    bnds=1
if int(bnds)>=int(len(bonds)+1):
    bnds=len(bonds)+1

#Generates list of combonications of fragments from 1 to all bonds matching the SMART code above.
combo=[]
from tqdm import tqdm
for i in tqdm(range(1,int(bnds)+1)):
    temp=(list(itertools.combinations(bonds,i)))
    mb = [list(elem) for elem in temp]
    combo.append(mb)
    
count = 0
for i in combo:
    count = count + len(i)
print(count)
print('#################################################################################')
#Generates string of all Smiles generated from every bond break combination created above.  The exception allow for the "iteration" to occur if only one bond is found.
ok=''
from tqdm import tqdm

cur=0
for frags in combo:
    for group in frags:
        percent = (cur/count * 100)
        print(percent)
        try:
            nm = Chem.FragmentOnBonds(m,group)
            ok= Chem.MolToSmiles(nm,True)+'.'+ok+'.'
            cur = cur +1
        except AttributeError:
            nm = Chem.FragmentOnBonds(m,combo)
            ok= Chem.MolToSmiles(nm,True)+'.'+ok+'.'
            cur = cur +1
        except ValueError:
            ok = ok+'.'
            cur = cur +1

#Data cleanup and duplicate removal.  Cleans out the atoms on break points and changes them from * to proton.
print('Cleaning Data...')    
ok2 = re.sub(r'\[(?:[^\]|]*\|)?([^\]|]*)\*\]','[H]', ok)
ok3 = re.sub(r'\*','[H]',ok2)    
Stack = [str(x) for x in ok3.split('.')]
Stack2 = [x for x in Stack if x]
Stack2.sort()
Stack.clear()
Stack3=list(Stack2 for Stack2,_ in itertools.groupby(Stack2))
Stack3.append(mcy)
Stack2.clear()

#Calculate Exact mass and prepare output csv.
mass=[]
for fragment in Stack3:
    mass.append([Descriptors.ExactMolWt(Chem.MolFromSmiles(fragment)),fragment])

    
with open(fragments+'.csv', 'w', newline='') as csvfile: 
    csvwriter = csv.writer(csvfile)
    headers=['Calculated Mass','Name']
    # writing the headers 
    csvwriter.writerow(headers) 
    # writing the data rows 
    csvwriter.writerows(mass) 

#Fragment Drawing Output, leave hashed out for now.  some valueerror or something.

#from rdkit.Chem import Draw#
#
#for name in Stack3:
#    try:
#        m = Chem.MolFromSmiles(name,sanitize=False)
#        Chem.SanitizeMol(m)
#        fig = Chem.Draw.MolToMPL(m, size=(200, 200), kekulize=True, wedgeBonds=True, firImage=True)
#    except AttributeError:
#        print(name)
#       
input("Press enter to exit.")