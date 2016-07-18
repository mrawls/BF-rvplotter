BF-rvplotter
============

I assume you have some spectra of a hopefully double-lined binary star, saved either as 1D FITS files or simple two-column text files, and would really like RV curves.

**Step 0**: Read more about BFs [here](http://www.astro.utoronto.ca/~rucinski/BFdescription.html) and [here](http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html).

Everything here _should_ work with Python 3. All the packages used come with [Astroconda](http://astroconda.readthedocs.io/en/latest/) except for [PyAstronomy](https://pyastronomy.readthedocs.io/en/latest/). If you're starting from scratch, that's how I recommend installing Python and getting everything up and running.

Note this is **not** a standalone python package, but rather a collection of handy programs.

###The main programs:
* `BF_python.py` is a radial velocity curve extractor using the Broadening Function technique based on Rucinski's `BFall_IDL.pro`
* `BF_functions.py` is used by `BF_python.py`
* `rvplotmaker.py` uses the output file from `BF_python.py` to make a pretty plot of your final RV curve

In practice, you will edit the input/output filenames and various parameters at the top of `BF_python.py` for each system you use it on. If you fork this repo and sync it with the master branch (recommended!) instead of downloading individual programs, you will want to save a record of your own personal "header" for use with `BF_python.py`.

###Bonus other stuff in `auxiliary` directory:
* `splot.py` makes a simple plot of FITS spectra
* `telfit_aband.py` creates a model telluric spectrum of the O2 A-band and/or uses the BF technique to calculate velocity shifts. Useful if you forgot to take standard spectra or have a fussy spectrograph with a changing wavelength solution.
* `telfit_plot.py` plots a model telluric spectrum on top of an actual spectrum
* `specshift.py` applies velocity shifts to a set of FITS spectra
* `spec_convolve.py` convolves a spectrum to a resolution of your choice using FORTRAN convspec (bonus: the convspec subdirectory includes a compiled FORTRAN binary which works on my OS X 10.10.5!)
* `spec_rotbroad.py` simulates rotational broadening given a vsini
* `spectra_rename.py` is a specialized script to move observations that are sorted by date into directories sorted by object
* `rv_latextable.py` turns your RV curve data into something resembling a LaTeX table
* `BF_functions_better.py` is a work in progress to make `BF_functions.py` more object-oriented, and `specload_test.py` is a short program to test its functionality

Walkthrough
============

So you have a bunch of spectra of a presumably double-lined eclipsing binary, and you want a nice radial velocity (RV) curve via the broadening function (BF) technique. This process is written specifically for a set of optical ARCES echelle spectra which have been reduced and combined into a single order, but it can work for any kind of 1D spectra, such as APOGEE IR spectra.

####1. Get an appropriate template/model file.

  * **template_spectrum.fits** (or .txt) - either a model spectrum with stellar parameters vaguely close to the target, or a real star spectrum also vaguely similar to the target. Either way, you need similar (or higher) resolution and the same wavelength coverage as the target. If you use a real star, you need to know its barycentric velocity (due to the Earth/Sun orbit) as well as its systemic velocity (relative to our Solar System) and remove these first (`specshift.py` can help!).

####2. Make some text files.

  * **infilesBF.txt** - a list of full paths to all the observed spectra. No comments or shenanigans in this file. The first one is the template.

  * **bjdinfile.txt** - three columns: filename (this string is not used, it is just a reference column), observation time (I recommend the BJD of exposure midpoint), barycentric velocity (hopefully from the header? that's on you). The zeroth row of this file MUST correspond to the `template_spectrum` file you plan to use (e.g., the first line could read template 0 0).

  * **gaussparsBF.txt** - a list of starting guesses for velocities so the gaussian fitter can find them. Must be the same length as bjdinfile.txt with six columns (it fits two gaussians with three parameters each by default): amp1, shift1, width1, amp2, shift2, width2. The first line can just be six zeros since it corresponds to the template. You can optionally fit a third gaussian by adding amp3, shift3, width3 to any line. You will edit this file after you run `BF_python.py` and are looking at the BF peaks, then run `BF_python.py` again until all the BF peaks are being fit nicely.

####3. Prepare some python scripts.

  * **BF_python.py** - does broadening function magic. You must specify BF parameters and the binary's period and zeropoint in the code itself.
INPUT: `infiles_BF.txt`, `bjdinfile.txt`, `gaussparsBF.txt`, `BF_functions.py`.
OUTPUT: text file with radial velocities! Now with a useful commented-out header, too.

  * **rvplotmaker.py** - plots the radial velocity curve as a function of time AND orbital phase. Put the star name in the code itself if you want it printed on the plot correctly.
INPUT: output from `BF_python.py`.
OUTPUT: plot only.

####4. Stand back and attempt science.

  * Run the programs in the order listed above. You'll need to run `BF_python.py` more than once as you adjust the values in `gaussparsBF.txt`.


The option exists to use telluric lines to line everything up first in case you don't have RV standard stars. This would apply if you don't have RV standard star observations and you want to use the telluric A-band to nudge all the spectra into something closer to an RV standard alignment. In this case, you would first run `telfit_aband.py` followed by `specshift.py`. Then you would proceed with the main steps above.

**telfit_aband.py** - OPTIONAL (calculates small RV offsets based on a telluric A-band model)
INPUT: `infiles_telfit.txt`, `bjdinfile.txt`, `gausspars_telfit.txt`.
OUTPUT: shows you some plots, prints shifts to screen, and also saves shifts in the new `shifts_telfit.txt` file.

**specshift.py** - OPTIONAL (applies velocity shifts to spectra)
INPUT: `shifts_telfit.txt`.
OUTPUT: new shifted FITS spectra in working directory.

**shifts_telfit.txt** - OPTIONAL (only if you are doing a telluric line velocity correction). This file should be a list of observed spectra (full filenames) followed by small velocity shifts that need to be removed. This is the only input for `specshift.py`, and is automatically created when you run `telfit_aband.py` (which requires `infiles_telfit.txt`, `bjdinfile.txt`, and `gausspars_telfit.txt`; these are analogous to the infiles used for a regular BF run).
