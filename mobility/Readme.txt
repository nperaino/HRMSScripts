This script will allow for the PA and PSA calculations.

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
more simply than in Fortran.

This is a rewrite of "Sigma" in python.

   T. Wyttenbach, G. von Helden,... M.T. Bowers; JASMS 1997, 8,
   275 
   
via the FORTRAN code provided by Thomas Wyttenbach.

The goal is to make a community accessible module for advancement of ion-mobility research.


There are three modes built up from a simple van der waals radius approximation to molecule size
to a temperature and size dependent super position correction.  The "best" model is to apply the
'size and temperature' mode which will acount for size deviations as a result of... size and temperature.

The PSA algorithm is implemented here but is much slower than standard PA. You will get the same answer for small
molecules as you would for large ones, so it is a waste of time for molecules under ~1000 atoms.
It is implemented here using the shape factor derived from the ratio of the convexhull/concavehull area
and the collision probability tailing function is from:

"A novel projection approximation algorithm for the fast and 
accurate computation of molecular collision cross sections.""
(I) Method
Christian Bleiholder, Thomas Wyttenbach, Michael T. Bowers
https://doi.org/10.1016/j.ijms.2011.06.014
(II) Model parameterization and definition of empirical shape factors for proteins
Christian Bleiholder, Stephanie Contreras, Thanh D. Do, Michael T. Bowers
https://doi.org/10.1016/j.ijms.2012.08.027

I am not sure how accurate it is yet, and the alpha value for the alpha-shape function (concave hull) is a bit uncertain to me.
It seemed like 6 made sense because it maximized the calculated area for c60 and insulin, but this could be completely wrong.
Still need to round out the script a bit but the main functions are working well.

