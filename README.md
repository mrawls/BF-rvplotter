BF-rvplotter
============

I assume you have some spectral observations of a hopefully double-lined binary star and would really like RV curves.

* `BF_python.py` is a radial velocity curve extractor using the Broadening Function technique. It is based on Rucinski's BFall_IDL.pro and uses the PyAstronomy tools, among other fun packages.

See more about BFs here: http://www.astro.utoronto.ca/~rucinski/BFdescription.html
and here: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

* `BF_functions.py` is used by `BF_python.py`
* `rvplotmaker.py` makes a pretty plot of your final RV curve
* `splot.py` makes a simple plot of FITS spectra
* `telfit_aband.py` creates a model telluric spectrum of the O2 A-band and/or uses the BF technique to calculate velocity shifts. Useful if you forgot to take standard spectra or have a fussy spectrograph with a changing wavelength solution.
* `telfit_plot.py` plots a model telluric spectrum on top of an actual spectrum
* `specshift.py` applies velocity shifts to a set of FITS spectra
* `spec_convolve.py` convolves a spectrum to a resolution of your choice using FORTRAN convspec (bonus: the convspec directory includes a compiled FORTRAN binary which works on my OS X 10.10.5!)

Walkthrough
============

So you have a bunch of spectra of a presumably double-lined eclipsing binary, and you want a nice radial velocity (RV) curve via the broadening function (BF) technique. Let's further assume you don't have RV standard star observations and you want to use the telluric A-band to nudge all the spectra into something closer to an RV standard alignment.

This process is written specifically for a set of optical ARCES echelle spectra which have been reduced and combined into a single order, but it can work for any kind of 1D spectra, such as APOGEE IR spectra (skip the telluric line alignment part in that case).

**(1) Get some appropriate template/model files.**

**telluricspec.txt** (or .fits) - probably created with the telfit MakeModel

**templateBF.fits** (or .txt) - either a model spectrum with stellar parameters vaguely close to the target, or a real star spectrum also vaguely similar to the target. Either way, you need similar (or higher) resolution and the same wavelength coverage as the target. If you use a real star, you need to know its barycentric velocity (due to the Earth/Sun orbit) as well as its systemic velocity (relative to our Solar System) and remove these first (`specshift.py` can help!).

**(2) Make some text files.**

**infiles_telfit.txt** - a list of full paths to all the observed spectra (easy!). No comments or shenanigans in this file.

**bjdinfile.txt** - three columns: image (filename does not have to be correct), observation time (BJD of exposure midpoint), barycentric velocity (hopefully from the header or something, figure it out). The first row MUST correspond to the templateBF.fits file you plan to use (e.g., template 0 0).

**gausspars_telfit.txt** - a list of guessed shift velocities so the gaussian fitter can find them. Must be the same length as bjdinfile.txt with six columns (it fits two gaussians with three parameters each by default): amp1, shift1, width1, amp2, shift2, width2. Set all shift1 = 0 and shift2 = something far from 0, like 20, so the second gaussian is out of the way since we don't need it here. The first line should just be six zeros since it corresponds to the template.

**shifts_telfit.txt** - a list of observed spectra followed by small velocity shifts that need to be removed. This is the only input needed for `specshift.py`, and is automatically created when you run `telfit_aband.py` (requires `telluricspec.txt`, `bjdinfile.txt`, `infiles_telfit.txt`, and `gausspars_telfit.txt`). The shift velocities are also spit out to the screen. Running specshift.py will then apply the shifts and create new FITS files in the working directory.

**infiles_BF.txt** - similar to `infiles_telfit.txt`, but now it's a list of the filenames spat out by `specshift.py`, and the top row must be `templateBF.fits`. No comments or shenanigans in this file. (If you want to skip all the telfit/specshift business because your target spectra are magically on the same RV standard already, this can just be `templateBF.fits` followed by the original spectrum filenames.)

**gaussparsBF.txt** - similar to `gausspars_telfit.txt`, only now we're guessing actual RVs for shift1 and shift2 instead of "0 and get out of the way." Feel free to copy `gausspars_telfit.txt` to get started, edit this file after you run `BF_python.py` and are looking at the BF peaks, then run `BF_python.py` again until all the BF peaks are being fit nicely.

**(3) Prepare some python scripts.**

**telfit_aband.py** - calculate small RV offsets based on telluric A-band model.
INPUT: `infiles_telfit.txt`, `bjdinfile.txt`, `gausspars_telfit.txt`.
OUTPUT: shows you some plots, prints shifts to screen, and also saves shifts in the new `shifts_telfit.txt` file.

**specshift.py** - applies velocity shifts to spectra.
INPUT: `shifts_telfit.txt`.
OUTPUT: new shifted FITS spectra in working directory.

**BF_python.py** - does broadening function magic. You must specify BF parameters and the binary's period and zeropoint in the code itself.
INPUT: `infiles_BF.txt`, `bjdinfile.txt`, `gaussparsBF.txt`, **BF_functions.py**.
OUTPUT: text file with radial velocities! The columns are BJD_mid, phase, Kep_BJD, RV1, err, RV2, err, source.

**rvplotmaker.py** - plots the radial velocity curve as a function of time AND orbital phase. Put the star name in the code itself if you want it printed on the plot correctly.
INPUT: output from `BF_python.py`.
OUTPUT: plot only.

**(4) Stand back and attempt science.**

Run the programs in the order listed above. You'll need to run `BF_python.py` more than once as you adjust the values in `gaussparsBF.txt`.
