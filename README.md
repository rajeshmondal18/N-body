# N-body
_____________________________________
A simple set of codes to simulate (semi-numerically) HI maps during reionization.

Use the makefile for compilation in the following manner:

make ionz_main

It will create the executable 'ionz_main'

You need to install FFTW 2.x.x with following flags: '--enable-float' and '--enable-type-prefix' to compile this set of codes.

Please acknowledge these papers 
1. Bharadwaj and Srikant 2004 (http://adsabs.harvard.edu/abs/2004JApA...25...67B), if you are using serial version of the code
2. Mondal et al. 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.449L..41M), if you are using parallel version of the code.
