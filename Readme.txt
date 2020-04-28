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