BF-rvplotter
============

I assume you have some spectral observations of an SB2 binary star and would really like RV curves.

* BF_python.py is a radial velocity curve extractor using the Broadening Function technique. It is based on Rucinski's BFall_IDL.pro and uses the PyAstronomy tools, among other fun packages.

See more about BFs here: http://www.astro.utoronto.ca/~rucinski/BFdescription.html
and here: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

<b>IMPORTANT UPDATE</b>: The RV axis was being defined incorrectly in earlier versions of BF_python.py. This has been fixed as of January 2015.

* BF_functions.py is used by BF_python.py
* rvplotmaker.py makes a pretty plot of your final SB2 RV curve
* splot.py makes a simple plot of FITS spectra
* telfit_aband.py creates a model telluric spectrum of the O2 A-band and/or uses the BF technique to calculate velocity shifts. Useful if you forgot to take standard spectra or have a fussy spectrograph with a changing wavelength solution.
* telfit_plot.py plots a model telluric spectrum on top of an actual spectrum
* specshift.py applies velocity shifts to a set of FITS spectra
* spec_convolve.py convolves a spectrum to a resolution of your choice using FORTRAN convspec (bonus: the convspec directory includes a compiled FORTRAN binary which works on my OS X 10.10.5!)
