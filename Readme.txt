How to run:
Download everything to a file and load up the MainPySigma script.  This will ask for whether you want to run a single temperature
or generate a curve.  There are the various methods and parameters developed over time included, but I have it set to only run the
temperature and size corrected method

Open an input file and select the settings you want to run at as prompted.
If a temperature curve is selected, it will save the data as a csv for graphing elsewhere.
A summary graph will output at the end if you just want to check quickly.
A run at 1% accuracy and 50 iterations takes only a minute or so for a temperature curve of C60.inp example.

There are two scripts for analyzing compounds HRMSMS and finding them in a spectrum.  The first script is the fragmentor which takes a SMILES input and outputs a list of fragments based on SMARTS codes. This list is formatted to enter into the Spectrum Search script, assuming you have exported the exact mass from an MSMS spectrum that might be from the first SMILES.  It will iterate over the lists to find which fragment masses are in the spectrum.  You need to check the structures to see that they make sense, because it does not decided whether or not something like an Oxonium or protonation event is possible, it will just find the mass and list it.  

Current work is incorporating these scripts into a GUI for analysis of mzXML files and will be made public once it is more completeed.  Currently you can view the mzXML chromatogram and mass spectra at retention times, and the fragmentor generates some structures to test your SMARTS code on before committing a 5-10 minute search function.

Notes:

Data must be formatted as in the included test data.
Be carfeul to make sure that if you are using excel, there are no "blank" columns by selecting the empty columns of cells and deleting them if you are having an issue.

Export centroided data and paste into Excel from Xcalibur and delete the header information.

A file labeled as adducts.csv must be in the Fragmentor folder for the spectrum search to function.  It must be formatted as in the example or it wont parse correctly.  I have included two, one for a lot of adducts and one with only two.  The one with two should be used for search MSMS spectra.  You must change the file name to change adducts at this time.  You can add stuff to the list too if you want to do conjugations or reasonable reactions like a methylation event or dehydration, it just generates all combinations of what is on the list with your targets.

The fragmentor generates MSMS data while spectrum search generates hits on a spectrum based on the adduct list.

The input files should be CSV for the SpectrumSearch and can be made by clicking and dragging accross the whole chromatogram and then right clicking and exporting exact mass on the mass spectrum. 

Then paste into Excel and save as CSV with the header removed so that it only has two columns: Mass and Intensity.

The other input for Spectrum Search is the Target list, designed for parent masses.  This input requires a CSV that contains the Name and Exact mass (unadducted) of the compound of interest.  In the case of MSMS searching you can directly use the output from fragmentor and the output from spectrum search will give a list of SMILES hits that can be converted in Excel using Chemdraw to give structures.  You must use 32 bit excel to use this unfortunately, which may require a re-install.

you can easily generate SMILES using chemdraw and then highlight the molecule and press Ctrl+alt+C, then go to the window and paste it.

You often have to click on the cmd window again after saving a file or something to allow for entering data into the window.   It is especially frustrating when using spectrum search, but just click on the cmd window again and press enter to open each file.

You shouldn't really need to break more than 4 combination of bonds, and more than 6 probably wont work because this is something like 10MM combinations if you have more than 50 bonds, and a lot of the resultant products will be redundant (ie as in breaking a peptide chain, if you break the first bond you get the first amino acid say, then in a combination of the first and 5th bond, you would also get the first amino acid).
