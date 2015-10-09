from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl
from scipy import ndimage
import pandas as pd
import gaussfitter as gf
import BF_functions as bff
'''
Program to extract radial velocities from a double-lined binary star spectrum.
Uses the Broadening Function technique.

Meredith Rawls
2014-2015

Based loosely on Rucinski's BFall_IDL.pro, and uses the PyAstronomy tools.
http://www.astro.utoronto.ca/~rucinski/BFdescription.html
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

In practice, you will run this twice: once to do the initial BF, and then again
to properly fit the peaks of each BF with a Gaussian.

INPUT
infiles:	single-column file with one FITS or TXT filename (w/ full path) per line
			1st entry must be for the template star (e.g., arcturus or phoenix model)
			(the same template is used to find RVs for both stars)
			NO comments are allowed in this file
			FUN FACT: unless APOGEE, these should be continuum-normalized to 1 !!!
bjdinfile: 	columns 0,1,2 must be filename, BJD, BCV (e.g., from IRAF bcvcorr)
			top row must be for the template star (e.g., arcturus)
			(the 0th column is never used, but typically looks like infiles_BF.txt)
			one line per observation
			comments are allowed in this file using #
gausspars:	your best initial guesses for fitting gaussians to the BF peaks
			the parameters are [amp1, offset1, width1, amp2, offset2, width2]
			the top line is ignored (template), but must have six values
			one line per observation
			comments are allowed in this file using #

OUTPUT
outfile:	a file that will be created with 8 columns: BJD midpoint, orbital phase,
			Kepler BJD, RV1, RV1 error, RV2, RV2 error, and source (a string)

IMMEDIATELY BELOW, IN THE CODE
You need to specify whether you have APOGEE (near-IR) or "regular" (e.g., ARCES)
spectra with the 'isAPOGEE' flag. You also need to set the binary's PERIOD and BJD0,
both in days, and the constant RV and BCV of whatever template you are using.
'''

##########
# YOU NEED TO HAVE THESE INPUT FILES
# THE OUTPUT FILE WILL BE CREATED FOR YOU

#infiles =   '../../RG_spectra/8430105/infiles_arcesBF.txt'
#bjdinfile = '../../RG_spectra/8430105/bjdinfile_arcesBF.txt'
#gausspars = '../../RG_spectra/8430105/gaussfit_arcesBF.txt'
#outfile =   '../../RG_spectra/8430105/rvoutfile2_arcesBF.txt'

# (the original infiles)
infiles = '../../TelFit/9246715_telfit/infiles_BF_shift.txt'
bjdinfile = '../../RG_spectra/9246715/bjds_baryvels.txt'
gausspars = '../../RG_spectra/9246715/gaussfit_pars.txt'
outfile = '../../RG_spectra/9246715/redo_plot_BFoutput.txt'

# STUFF YOU NEED TO DEFINE CORRECTLY !!!
isAPOGEE = False

# ORBITAL PERIOD AND ZEROPOINT !!!
period = 171.277697; BJD0 = 2455170.514777 # 9246715
#period = 63.327106; BJD0 = 2454976.635546 # 8430105
#period = 120.390971; BJD0 = 2454957.586519 # 10001167

# RADIAL VELOCITY AND BCV INFO FOR TEMPLATE (km/s; set both to 0 if using a model)
#rvstd = -64.422; bcvstd = 10.747 # HD168009 (fullspec.0026)
#rvstd = -21.619; bcvstd = 16.571 # HD182488 (fullspec.0028)
#rvstd = -21.123; bcvstd = 12.499 # HD196850 (fullspec.0032)
rvstd = 0; bcvstd = 0 # model template

# PARAMETERS FOR THE BROADENING FUNCTION (you can adjust w00, n, and stepV below)
amp = 6.0		    # arbitrary amplitude to stretch the smoothed BFs by in y, for clarity
smoothstd = 1.5     # stdev of Gaussian to smooth BFs by (function of instrument resolution)
m = 211             # length of the BF (must be longer if RVs are far from 0)
#m = 171

# STUFF TO MAKE PLOTS LOOK NICE
rvneg = -69; rvpos = 69; ymin = -0.05; ymax = 0.45 # 9246715
#rvneg = -69; rvpos = 69; ymin = -0.05; ymax = 0.35 # 8430105
#rvneg = -170; rvpos = 5; ymin = -0.05; ymax = 0.15 # 10001167

# some previously set values for posterity ...
# ARCES ARCTURUS OBSERVATION
#rvstd = 20.71053 # this is the TOTAL RV OFFSET FROM REST of the ARCES Arcturus observation
#bcvstd = 0 # this goes with the above rvstd
#rvstd = -5.19 # from SIMBAD, for Arcturus ... WRONG(ish)
#bcvstd = -0.155355148339 # APOGEE Arcturus bcv
#bcvstd = 18.4574 # ARCES Arcturus bcv
##########

print('Welcome to BF_python.py! It\'s broadening function time.')
print(' ')
print('You have set Porb = {0} days, BJD0 = {1}'.format(period, BJD0))
print(' ')

# CREATE NEW SPECTRUM IN LOG SPACE
# The wavelengths are w1.
# The BF will be evenly spaced in velocity with length m.
# The velocity steps are r (km/s/pix).
# Guidelines for choosing w00, n, stepV, and m are in the functions file.
if isAPOGEE == True:
	w00 = 15145
	n = 15000 #20000 #15000 ... testing apogee?
else:
	w00 = 5400          # 9246715 accommodates TRES
	n = 42500 #38750    # 9246715 accommodates TRES
	#w00 = 4000
	#n = 90000
stepV = 1.5
w1, m, r = bff.logify_spec(isAPOGEE, w00, n, stepV, m)

# READ IN ALL THE THINGS
nspec, filenamelist, datetimelist, wavelist, speclist, source = bff.read_specfiles(infiles, bjdinfile, isAPOGEE)

# INTERPOLATE THE TEMPLATE AND OBJECT SPECTRA ONTO THE NEW LOG-WAVELENGTH GRID
# OPTION TO PLOT THIS (commented out for now)
##plt.figure(1)
newspeclist = []
yoffset = 1
#plt.axis([w1[0], w1[-1], 0, 27])
#plt.xlabel('wavelength')
for i in range (0, nspec):
	newspec = np.interp(w1, wavelist[i], speclist[i])
	newspeclist.append(newspec)
#	plt.plot(w1, newspeclist[i]+yoffset, label=datetimelist[i].iso[0:10])
#	yoffset = yoffset + 1
##plt.legend()
#plt.show()

# BROADENING FUNCTION TIME
svd = pyasl.SVD()
# Single Value Decomposition
svd.decompose(newspeclist[0], m)
bflist = []
bfsmoothlist = []
for i in range (0, nspec):
	# Obtain the broadening function
	bf = svd.getBroadeningFunction(newspeclist[i]) # this one is like a matrix
	bfarray = svd.getBroadeningFunction(newspeclist[i], asarray=True)
	# Smooth the array-like broadening function
	bfsmooth = amp*pd.rolling_window(bfarray, window=5, win_type='gaussian', std=smoothstd, center=True)
	# The rolling window makes nans at the start because it's a punk.
	for j in range(0,len(bfsmooth)):
		if np.isnan(bfsmooth[j]) == True:
			bfsmooth[j] = 0
	bflist.append(bf)
	bfsmoothlist.append(bfsmooth)
# Obtain the indices in RV space that correspond to the BF
bf_ind = svd.getRVAxis(r, 1) + rvstd - bcvstd

# PLOT THE RESULT OF THE SINGLE VALUE DECOMPOSITION TO SEE IF IT IS TERRIBLE OR NOT
# NOT READY YET
#plt.figure(2)
#plt.plot(m, svd)
#plt.show()

# OPTION TO PLOT THE SMOOTHED BFs (commented out for now)
##plt.figure(3)
#plt.axis([rvneg, rvpos, -0.2, 12])
#plt.xlabel('Radial Velocity (km s$^{-1}$)')
#plt.ylabel('Broadening Function (arbitrary amplitude)')
#yoffset = 0.0
#for i in range(1, nspec):
#	plt.plot(bf_ind, bfsmoothlist[i]+yoffset, color='b')
#	yoffset = yoffset + 0.4
#plt.show()

# FIT THE SMOOTHED BF PEAKS WITH TWO GAUSSIANS
# you have to have pretty decent guesses in the gausspars file for this to work.
bffitlist = bff.gaussparty(gausspars, nspec, filenamelist, bfsmoothlist, bf_ind)
rvraw1 = []; rvraw2 = []; rvraw1_err = []; rvraw2_err = []
rvraw1.append(0), rvraw2.append(0), rvraw1_err.append(0), rvraw2_err.append(0)
for i in range(1, len(bffitlist)):
    rvraw1.append(bffitlist[i][0][1])
    rvraw2.append(bffitlist[i][0][4])
    rvraw1_err.append(bffitlist[i][2][1])
    rvraw2_err.append(bffitlist[i][2][4])

# CALCULATE ORBITAL PHASES AND FINAL RV CURVE
phase, bjdfunny, rv1, rv2, rv1_err, rv2_err = bff.rvphasecalc(bjdinfile, outfile, nspec, period, BJD0, rvraw1, rvraw1_err, rvraw2, rvraw2_err, rvstd, bcvstd, source)

# PLOT THE FINAL SMOOTHED BFS + GAUSSIAN FITS IN INDIVIDUAL PANELS
# manually adjust this multi-panel plot based on how many spectra you have
#plt.figure(4)
windowcols = 4		                        # how many window columns there should be
windowrows = np.rint(nspec/windowcols)+1	# how many window rows there should be
xmin = rvneg
xmax = rvpos
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Uncorrected Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=26)
fig.text(0.07, 0.5, 'Broadening Function', ha='center', va='center', size=26, rotation='vertical')
for i in range (1,nspec):
	ax = fig.add_subplot(windowrows, windowcols,i) # out of range if windowcols x windowrows < nspec
	ax.yaxis.set_major_locator(MultipleLocator(0.2))
	if i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21 and i!=25:
		ax.set_yticklabels(())
	#if i!=20 and i!=21 and i!=22 and i!=23 and i!=24 and i!=25:
	if i < nspec-windowrows:
	#if i!=13 and i!=14 and i!=15 and i!=16:
		ax.set_xticklabels(())
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.axis([xmin, xmax, ymin, ymax])
	plt.tick_params(axis='both', which='major', labelsize=14)
	plt.text(xmax - 0.25*(np.abs(xmax-xmin)), 0.8*ymax, '%.3f $\phi$' % (phase[i]), size=12)
	plt.text(xmax - 0.35*(np.abs(xmax-xmin)), 0.6*ymax, '%s' % (datetimelist[i].iso[0:10]), size=12)
	## comment out below if you don't want to print source on plot windows ##
	if source[i] == 'arces': plt.text(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, 'ARCES', color='#0571b0', size=12)
	elif source[i] == 'tres': plt.text(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, 'TRES', color = '#008837', size=12)
	elif source[i] == 'arces': plt.txt(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, 'APOGEE', color = 'k', size=12)
	#else: plt.text(0.9*xmin, 0.4*ymax, 'SOURCE?', color = 'k', size=12)
	## comment out above if you don't want to print source on plot windows ##
	plt.plot(bf_ind, bfsmoothlist[i], color='k', lw=1.5, ls='-', label='Smoothed BF')
	plt.plot(bf_ind, bffitlist[i][1], color='#e34a33', lw=2, ls='--', label='Two-Gaussian fit')
	#plt.axvline(x=0, color='0.75')	
	if i==23: ax.legend(bbox_to_anchor=(2.1,0.7), loc=1, borderaxespad=0., 
						frameon=False, handlelength=3, prop={'size':20})
plt.show()