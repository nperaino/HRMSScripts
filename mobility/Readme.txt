Input file is generated as an XYZ Gaussian input.  I use Gabedit to do geometry optimization and then create a gaussian input.
There is no specific reason I am doing it this way other than it is quick.

This will mimic the graph for C60 fullerene temperature dependent size data shown in:

Effect of the Long-Range Potential on Ion
Mobility Measurements
Thomas Wyttenbach, Gert von Helden, Joseph J. Batka, Jr., Douglas
Carlat, and Michael T. Bowers
J Am Soc Mass Spectrom
1997, 8, 275-282)

But significant subroutine editing has been done regarding removal of interpolation parameterization since SciPy can handel a lot of this
more simply than in Fortran.  I have not yet figured out why the buffer gas radius is off by a factor of 10.  So if buffr variable is set
to 1.09, it is wrong, but 0.109 gives the correct predictions.

This is a rewrite of "Sigma" in python.

   T. Wyttenbach, G. von Helden,... M.T. Bowers; JASMS 1997, 8,
   275 
   
via the FORTRAN code provided by Thomas Wyttenbach.

The goal is to make a community accessible module for advancement of ion-mobility research.


There are three modes built up from a simple van der waals radius approximation to molecule size
to a temperature and size dependent super position correction.  The "best" model is to apply the
'size and temperature' mode which will acount for size deviations as a result of... size and temperature.
There is not a mode for just size, as development first started with the simple model then added temperature,
then size, and then in the projections superposition approximation, a shape function for concavity was added 
(not here yet).
