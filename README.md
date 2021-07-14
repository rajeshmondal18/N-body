# N-body

This is an open-source code to compute the non-linear evolution of the cosmological matter density contrast. It is based on the Particle Mesh (PM) technique. It is written in C. The code has been parallelized (Mondal et al. 2015) for shared-memory machines using Open Multi-Processing (OpenMP).

Read the user's guide 'nbody_doc.pdf' for a detailed description.
_____________________________________

Download the code by cloning the git repository using

$ git clone https://github.com/rajeshmondal18/N-body.git

You need to install FFTW-3.x.x with the following flags: '--enable-float',  '--enable-threads' and '--enable-openmp' to compile this set of codes.
Look at the installation instruction http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix
_____________________________________
Use the makefile for compilation in the following manner:

$ make nbody_comp

It will create the executable 'nbody_comp'.

To Run:

$ ./nbody_comp
_____________________________________
Please acknowledge these papers 
1. Bharadwaj and Srikant 2004 (http://adsabs.harvard.edu/abs/2004JApA...25...67B), if you are using the serial version of the code
2. Mondal et al. 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.449L..41M), if you are using the current parallelized version of the code.
