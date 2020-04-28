#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#This script calculates the ESI adducts for a list of target masses and finds
#them in the spectrum csv.

#prepared by Nicholas Peraino, Lumigen Instrument Center, Wayne State University


import csv
import tkinter as tk

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Spectrum Search 1.10')
print('')
print('This script searches a spectrum for files on a mass list accounting for various ESI adducts')
print('Input files must be properly formatted to parse correctly')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
print('')
print('')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Target List should be formatted as a csv file containing Target and Name as a header, followed by the unadducted exact mass and a name')
print('Example:')
print('')
print('        Target, Name')
print('        178.0477, Biotin')
print('        441.1396, Folic Acid')
print('        etc.')
print('')
print('')
print('')
print('Please select the list of targets csv file.')
input('Press Enter to select a file')
from tkinter.filedialog import askopenfilename
root = tk.Tk()
root.withdraw()
masslistfile = askopenfilename()
print(masslistfile)
print('')
print('')
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Spectrum file should be formatted as csv file containing Mass and Intensity as a header, followed by the observed mass and intensities')
print('Example:')
print('')
print('        Mass, Intensity')
print('        320.1150, 100')
print('        321.1188, 11.9')
print('        322.1120, 4.5')
print('')
print('')
print('Please select the spectrum csv file.')
input('Press Enter to select a file.')
from tkinter.filedialog import askopenfilename
spectrumfile = askopenfilename()
print(spectrumfile)
print('')
print('')
print('Please select an output folder for the data.')
input('Press Enter to select a file')
from tkinter.filedialog import asksaveasfilename
hits = asksaveasfilename()
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('--------------------------------------------------------------------------------------------------------------------------------------')
masserror = abs(float(input("Enter a value for mass error (ppm) accepted: ")))
minintensity = abs(float(input("Enter a value for minimum intensity accepted: ")))
chargeerror = abs(float(input("Enter a value for acceptable charge variance for confimation by isotope distribution: ")))
#This will ask for user to open file
#from tkinter.filedialog import askopenfilename + tkinter.Tk().withdraw()
#filename = askopenfilename()



rawtarget=[]
with open(masslistfile) as csvfile:
        mass_list = csv.reader(csvfile, delimiter=",")
        for row in mass_list:
            mlh=len(row)                
            for i in row:
                try:
                    rawtarget.append(float(i))
                except ValueError:
                    rawtarget.append(str(i))

#print(rawtarget[2::2]) list of targeted masses
#print(rawtarget[3::2]) list of targeted mass names
                    
                    
                    
                    
#-----------------------------------------------------------------------#
#load the adducts and convert to floats
#-----------------------------------------------------------------------#



rawadductlist=[]
with open('adducts.csv') as csvfile:
        adducts = csv.reader(csvfile, delimiter=",")
        for row in adducts:
            alh=len(row) 
            for i in row:
                try:
                    rawadductlist.append(float(i))
                except ValueError:
                    rawadductlist.append(str(i))
                    
#print(rawadductlist[4::4]) list of adduct names                    
#print(rawadductlist[5::4]) list of adduct charges
#print(rawadductlist[6::4]) list of adduct mass multipliers
#print(rawadductlist[7::4]) list of adduct masses


#-----------------------------------------------------------------------#
#load the spectrum data and convert to floats
#-----------------------------------------------------------------------#



rawspectrum=[]
with open(spectrumfile) as csvfile:
        spec = csv.reader(csvfile, delimiter=",")
        for row in spec:
            slh=len(row)
            for i in row:
                try:
                    rawspectrum.append(float(i))
                except ValueError:
                    rawspectrum.append(str(i))
                except TypeError:
                    rawspectrum.append(str(i))


#------------------------------------------------------------------------#
#adduct the target masses and format files for named output.
#------------------------------------------------------------------------#                   
print('Preparing List of Adducts...')                    
                    
adducted2=[]
from tqdm import tqdm             
for i in tqdm(rawtarget[2::2]):
    for adduct, charge in zip(rawadductlist[7::4], rawadductlist[6::4]):
        adducted2.append(round((i+adduct)*charge,4))

numadducts=len(rawadductlist[4::4])
adducted = [adducted2[i:i+numadducts] for i in range(0, len(adducted2), numadducts)]


spectrum=[]
largespectrum=[]
del rawspectrum[0:2]
for mass,intensity in zip(rawspectrum[::2], rawspectrum[1::2]):
    largespectrum.append([(round(mass,4)),(round(intensity,2))])
  

for row in largespectrum:
    if row[1]>minintensity:
        spectrum.append(row)
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')

#--------------------------------------------------------------------------#
### can probably be improved but the the only way I could get the names to work was to set
### separate lists for the names and charges because of the range requirement.
#--------------------------------------------------------------------------#


ions=rawadductlist[4::4]
charges=rawadductlist[5::4]

print('Naming Targets... Step 1 of 4')

named=[]
from tqdm import tqdm 
for i in tqdm(rawtarget[3::2]):
    for j in range(numadducts):
        named.append(i)
        
## makes a list of target names with copies long enough to merge with the adducts
names=[named[i:i+numadducts] for i in range(0, len(named), numadducts)]
print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')
#merges adducted mass values with the ion names and charges    

nameprep1=[]
for j in range(len(rawtarget[3::2])):
    for i in range(0, numadducts):
        nameprep1.append([adducted[j][i],ions[i],round(charges[i])])

#merges the target names with the ion information.
print('Merging Target Identifiers...Step 2 of 4')

namedadduct=[]
from tqdm import tqdm 
for i,j in tqdm(zip(nameprep1,named)):
    namedadduct.append([j,i[0],i[1],i[2]])
    


#--------------------------------------------------------------------------#
#Search the spectrum against the mass list to find hits.
#--------------------------------------------------------------------------#


print('Searching spectrum within', masserror, 'ppm mass error and', minintensity, 'minimum intensity... Step 3 of 4')
hitlist=[]
from tqdm import tqdm 
for row in tqdm(namedadduct):
    for i in range(len(spectrum)):
       if abs((spectrum[i][0]-row[1])/row[1]*1000000)<masserror:
           hitlist.append(row + spectrum[i] + [round(((spectrum[i][0]-row[1])/row[1]*1000000),2)])

print('')
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')

print('Matching charge distribution... Step 4 of 4')

isotopelist=[]
from tqdm import tqdm
for row in tqdm(hitlist):
    for i in range(len(spectrum)):
        try:
            if (1/(spectrum[i][0]-row[4]))>0 and abs(abs(1/(spectrum[i][0]-row[4]))-row[3])<chargeerror and (row[5]/spectrum[i][1])<20:
                isotopelist.append(row + spectrum[i] + [round(abs(1/(spectrum[i][0]-row[4])),2)])
        except ZeroDivisionError:
            continue
#--------------------------------------------------------------------------#
#Sorting by intensity.  Change val to change sorting as 1=calculated mass
#--------------------------------------------------------------------------#
           
         
            
def sortintensity(val): 
    return val[5]  
hitlist.sort(key = sortintensity, reverse = True)  

def sortintensity(val): 
    return val[5]  
isotopelist.sort(key = sortintensity, reverse = True)  

with open(hits+'.csv', 'w', newline='') as csvfile: 
    # creating a csv writer object 
    csvwriter = csv.writer(csvfile)
    headers=['Name','Calculated Mass', 'Ion Identity','Charge', 'Found Mass', 'Intensity', 'PPM Error']
    isohead=['Name','Calculated Mass', 'Ion Identity','Charge', 'Found Mass', 'Intensity', 'PPM Error', 'Isotope Mass', 'Isotope Intensity', 'Apparent Charge']
    # writing the headers 
    csvwriter.writerow(headers) 
    # writing the data rows 
    csvwriter.writerows(hitlist) 
    blank1=[]
    csvwriter.writerow(isohead)
    csvwriter.writerows(isotopelist)


#--------------------------------------------------------------------------#
#Prepare output as string and format column widths.
#--------------------------------------------------------------------------#
print('') 
print('')   
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')    
print('This is the list of hits on the list, unfiltered by isotope charge matching.')
print('')  
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('') 
print('') 
strlistprep=[]
for row in hitlist:
    for item in row:
        strlistprep.append(str(item))
headers.extend(strlistprep)
try:
    strlist= [headers[i:i+len(hitlist[0])] for i in range(0, len(headers), len(hitlist[0]))]
except IndexError:
    strlist=[['No Targets Found']]
widths = [max(map(len, col)) for col in zip(*strlist)]

for row in strlist:
    print("  ".join((val.ljust(width) for val, width in zip(row, widths))))
print('') 
print('') 
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('')    
print('This is the list of hits matching the expected charge from isotope distribution.')
print('')  
print('--------------------------------------------------------------------------------------------------------------------------------------')
print('') 
print('') 
isoprep=[]
for row in isotopelist:
    for item in row:
        isoprep.append(str(item))

isohead.extend(isoprep)

try:
    isolist= [isohead[i:i+len(isotopelist[0])] for i in range(0, len(isohead), len(isotopelist[0]))]
except IndexError:
    isolist=[['No Target Isotopes Found']]
    
widths = [max(map(len, col)) for col in zip(*isolist)]

for row in isolist:
    print("  ".join((val.ljust(width) for val, width in zip(row, widths))))
input("Press enter to exit.")