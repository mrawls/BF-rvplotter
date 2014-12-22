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
'''
Program to extract radial velocities from a double-lined spectrum.
Uses the Broadening Function technique.

Meredith Rawls
Summer 2014

Based loosely on Rucinski's BFall_IDL.pro and uses the PyAstronomy tools.
See more here:
http://www.astro.utoronto.ca/~rucinski/BFdescription.html
and
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

There are lots of packages you need. Have fun with that.

In practice, you will run this twice: once to do the initial BF, and then again
to properly fit the peaks of each BF with a Gaussian.

There are lots of optional plots you can un-comment for sanity checks.

REQUIRED INFILES
infiles:		single-column file with one FITS filename (or path+filename) per line
			1st entry must be for the template star (e.g., arcturus)
			(note that the same template is used to find RVs for both stars)
			no comments are allowed in this file
			FUN FACT: unless APOGEE, these should be continuum-normalized to 1 !!!
bjdinfile: 	columns 0,1,2 must be FITS filename, BJD, BCV (from IRAF bcvcorr)
			top row must be for the template star (e.g., arcturus)
			(the 0th column is never used, but typically looks like infiles_BF.txt)
			one line per observation
			comments are allowed in this file using #
			subsequent columns after 0,1,2 are initially ignored**
gausspars:	your best initial guesses for fitting gaussians to the BF peaks
			the parameters are [amp1, offset1, width1, amp2, offset2, width2]
			the top line is ignored (template), but must have six values
			one line per observation
			comments are allowed in this file using #

Finally, you'll need to specify various parameters near the top of the code.

There are TWO True/False FLAGS you can set early in the code. One is explained below.
The other is for whether you are working with APOGEE near-IR spectra and/or "regular"
(i.e., ARCES or other optical) spectra.

You will need to tweak the number of panels being plotted at the end of the code,
depending on how many observations you have.

**OPTION TO MAKE FINAL PLOT NICE WITH STAR 1 & STAR 2 PROPERLY ASSIGNED
By default, the left peak is assigned 'star 1' and the right peak is assigned 'star 2.'
This is obviously not always correct.
This is how to fix it once you have good RV values, but some assigned to the wrong star.
-->	edit bjds_baryvels.txt columns 3,4,5,6 to be rvraw1, rvrawerr1, rvraw2, rvrawerr2
	once you look at the BF amplitudes and figure out which star is which in each case.
-->	set variable bjdfilehasrvs = True below (default is False)
-->	run this program again
'''

print('Welcome to BF_python.py! It\'s broadening function time.')
print(' ')

# YOU NEED TO HAVE THESE FILES
infiles = 'infiles2_BF.txt'
bjdinfile = 'bjds_baryvels2.txt'
gausspars = 'gaussfit_pars2.txt'

# FILE THAT WILL BE WRITTEN TO
outfile = 'rvs_out_apogee.txt'

# STUFF YOU NEED TO DEFINE CORRECTLY !!!
bjdfilehasrvs = False
isAPOGEE = True
period = 171.277967
BJD0 = 2455170.514777
rvstd = -5.19 # from SIMBAD, for Arcturus
#bcvstd = -0.155355148339 # APOGEE Arcturus observation TIME DEPENDENT!!!
bcvstd = 18.4574 # ARCES Arcturus observation

# The new log-wavelength array will be w1. it will have equal spacing in velocity.
# please specify reasonable values below or else bad things will happen.
# note log = log base 10 because SERIOUSLY, come on now.
if isAPOGEE == True:
	w00 = 15145		# starting wavelength of the log-wave array in Angstroms
	n = 20000			# desired length of the log-wave vector in pixels (must be EVEN)
else:
	w00 = 5400 		# starting wavelength of the log-wave array in Angstroms
	n = 38750 		# desired length of the log-wave vector in pixels (must be EVEN) 
stepV = 1.7 			# step in velocities in the wavelength vector w1
m = 171 				# length of BF (must be ODD)
r = stepV/2.997924e5 	# put stepV in km/s/pix
w1 = w00 * np.power((1+r), np.arange(float(n)))
amp = 5.0 			# amplitude to stretch the smoothed BFs by (y-direction) for clarity
print ('The new log-wavelength scale will span %d - %d A with stepsize %f km/s.' % (w1[0], w1[-1], stepV))
print(' ')

# READ IN DATA FROM FITS FILES
# Read in a text file containing a list of fits files
# The first one listed will be the template spectrum
f1 = open(infiles)
print('Reading FITS files from 1st column in %s' % infiles)
print('The first one had better be your template spectrum.')
print(' ')
speclist = []; wavelist = []
filenamelist = []; datetimelist = []
source = [] # option: keep track of which spectrograph was used (ARCES vs. TRES)
i = 0
for line in f1: # This loop happens once for each spectrum (FITS file)
	infile = line.rstrip()
	# Read in the FITS file with all the data in the primary HDU
	hdu = fits.open(line)
	if isAPOGEE == True: # APOGEE: the data is in a funny place, backwards, and not normalized
		spec = hdu[1].data ### APOGEE
		spec = spec.flatten() ### APOGEE
		spec = spec[::-1] ### APOGEE
		spec = spec / np.median(spec)
	else: # non-APOGEE (regular) option
		spec = hdu[0].data
	head = hdu[0].header
#	# *** begin SPECIAL FOR MORE THAN ONE SPECTROGRAPH (ARCES + TRES) ONLY ***
#	# Jean says... could do "if '.tres.' in line:" / "if '.ec.' in line". meh.
#	if head['imagetyp'] == 'object': source.append('arces')
#	if head['imagetyp'] == 'OBJECT': source.append('tres')
#	# *** end SPECIAL FOR MORE THAN ONE SPECTROGRAPH (ARCES + TRES) ONLY ***
	filenamelist.append(infile)
	datetime = head['date-obs']
	datetimelist.append(Time(datetime, scale='utc', format='isot'))
	# Define the original wavelength scale
	if isAPOGEE == True: # APOGEE: read wavelength values straight from FITS file
		wave = hdu[4].data ### APOGEE
		wave = wave.flatten() ### APOGEE
		wave = wave[::-1] ### APOGEE
	else: # non-APOGEE (linear) option: create wavelength values from header data
		headerdwave = head['cdelt1']
		headerwavestart = head['crval1']
		headerwavestop = headerwavestart + headerdwave*len(spec)
		wave = np.arange(headerwavestart, headerwavestop, headerdwave)
	if len(wave) != len(spec): # The wave array is sometimes 1 longer than it should be?
		minlength = min(len(wave), len(spec))
		wave = wave[0:minlength]
		spec = spec[0:minlength]
	wavelist.append(wave)
	speclist.append(spec)
	i = i + 1	
# save the total number of spectra
nspec = i
f1.close()

# INTERPOLATE THE TEMPLATE AND OBJECT SPECTRA ONTO THE NEW LOG-WAVELENGTH GRID
# option to make a plot for sanity
newspeclist = []
yoffset = 1
#plt.axis([w1[0], w1[-1], 0, 27])
#plt.xlabel('wavelength')
for i in range (0, nspec):
	newspec = np.interp(w1, wavelist[i], speclist[i])
	newspeclist.append(newspec)
#	plt.plot(w1, newspeclist[i]+yoffset, label=datetimelist[i].iso[0:10])
	yoffset = yoffset + 1
#plt.legend()
#plt.show()

# DO A SINGLE VALUE DECOMPOSITION DANCE
svd = pyasl.SVD()
svd.decompose(newspeclist[0], m)

# BROADENING FUNCTIONS WOOOOOOOOO!!!
bflist = []
bfsmoothlist = []
for i in range (0, nspec):
	# Obtain the broadening function
	bf = svd.getBroadeningFunction(newspeclist[i]) # this one is like a matrix
	bfarray = svd.getBroadeningFunction(newspeclist[i], asarray=True)
	# Smooth the array-like broadening function
	bfsmooth = amp*pd.rolling_window(bfarray, window=5, win_type='gaussian', std=1.5)
	# The rolling window makes nans at the start because it's a punk.
	for j in range(0,len(bfsmooth)):
		if np.isnan(bfsmooth[j]) == True:
			bfsmooth[j] = 0
	bflist.append(bf)
	bfsmoothlist.append(bfsmooth)
# Obtain the indices in RV space that correspond to the BF
bf_ind = stepV*(np.arange(-m/2, m/2))

# OPTION TO PLOT THE SMOOTHED BFs
#plt.axis([-70, 70, -0.2, 9.5])
#plt.xlabel('Radial Velocity (km s$^{-1}$)')
#plt.ylabel('Broadening Function (arbitrary amplitude)')
#yoffset = 0.0
#for i in range(1, nspec):
#	plt.plot(bf_ind, bfsmoothlist[i]+yoffset, color='b')
#	#plt.plot(bf_ind, bflist[i]+yoffset) # unsmoothed BFs
#	yoffset = yoffset + 0.4
#plt.show()

# FIT THE SMOOTHED BF PEAKS WITH TWO GAUSSIANS
# you have to have pretty decent guesses for it to fit the main peaks.
# edit the file gaussfit_pars.txt to adjust these values.
f1 = open(gausspars)
param = np.loadtxt(f1, comments='#')
f1.close()
bffitlist = []
bffitlist.append(0)
gauss1 = [[] for i in range(nspec)]
gauss2 = [[] for i in range(nspec)]
gauss1[0] = [0,0]
gauss2[0] = [0,0]
rvraw1 = []; rvraw2 = []
rvraw1_err = []; rvraw2_err = []
rvraw1.append(0); rvraw2.append(0); rvraw1_err.append(0); rvraw2_err.append(0)
error_array = np.ones(len(bfsmoothlist[0]))*0.01 # dummy array with 0.01 error values
print('Two-Gaussian fit results: peak amplitude, rvraw, rvraw_err (for each star*)')
print ('---------------------------------------------------------------------------')
for i in range(1, nspec):
	bffit = gf.multigaussfit(bf_ind, bfsmoothlist[i], ngauss=2, 
			params=param[i], err=error_array,
			limitedmin=[True,True,True], limitedmax=[True,True,True], 
			minpars=[0.2,-70,1], maxpars=[0.6,70,5], quiet=True, shh=True)
	bffitlist.append(bffit)
	# NOTE: to get the gaussian fit corresponding to bfsmoothlist[i], use bffitlist[i][1].
	gauss1[i] = [bffit[0][0], bffit[0][1], bffit[2][1]] # these are [amp1, rvraw1, rvraw1_err]
	gauss2[i] = [bffit[0][3], bffit[0][4], bffit[2][4]] # these are [amp2, rvraw2, rvraw2_err]
	print ('%s \t %.5f %.5f %.5f \t %.5f %.5f %.5f' % (filenamelist[i][-30:], 
		gauss1[i][0], gauss1[i][1], gauss1[i][2], gauss2[i][0], gauss2[i][1], gauss2[i][2]))
	# NOTE: the following rvraw values will be overwritten below if bjdfilehasrvs = True
	rvraw1.append(bffit[0][1])
	rvraw2.append(bffit[0][4])
	rvraw1_err.append(bffit[2][1])
	rvraw2_err.append(bffit[2][4])
	
# OPTION TO PLOT DETAILED VIEW OF EACH SMOOTHED BF
# (you will need this for guesstimating the location of each Gaussian's peak)
#yoffset = 0.0
#for i in range(0, nspec-1):
#	plt.plot(bf_ind, bffitlist[i][1]+yoffset, color='r')
#	yoffset = yoffset + 0.4
#plt.show()
print(' ')
print('You MUST manually guesstimate the location of each Gaussian\'s peak in %s!' % gausspars)
print('Until you do, the above values will be WRONG and the plot will look TERRIBLE.')
print(' ')

# CALCULATE ORBITAL PHASES AND TRUE RV CURVE
rv1 = []; rv2 = []
rv1.append(0); rv2.append(0)
rv1_err = []; rv2_err = []
rv1_err.append(0); rv2_err.append(0)
print('*******')
print('*Note that Star 1 and Star 2 MAY NOT BE properly assigned above.')
print('For now, you must do this manually; that is why the peak amplitude is printed.')
print('You should update %s manually, set bjdfilehasrvs = True, and run again.' % bjdinfile)
g1 = open(bjdinfile)
g2 = open(outfile, 'w')
# Decide if we should calculate the RVs in place or read them from an updated bjdinfile
if bjdfilehasrvs == True:
### USE THIS if columns 3,4,5,6 of bjdinfile contain RVraw1, RVraw1err, RVraw2, RVraw2err
	print('NOT calculating RVs; reading them from bjd infile...')
	rvraw1 = []; rvraw2 = []
	rvraw1_err = []; rvraw2_err = []
	bjdmid, bcv, rvraw1, rvraw1_err, rvraw2, rvraw2_err = np.loadtxt(g1, comments='#', 
		dtype=np.float64, usecols=(1,2,3,4,5,6), unpack=True)
else:
### USE THIS otherwise (default option)
	print('Calculating RVs...')
	bjdmid, bcv = np.loadtxt(g1, comments='#', dtype=np.float64, usecols=(1,2), unpack=True)
bjdfunny = bjdmid - 2454833.
phase = []
phase.append(0)
for i in range(1, nspec):
	fracP = (bjdmid[i] - BJD0) / period
	if fracP < 0:
		phase.append(1 + (fracP % 1))
		cycle = int(fracP) - 1
	else:
		phase.append((fracP % 1))
		cycle = int(fracP)
	rv1.append(rvraw1[i] + bcvstd - bcv[i] + rvstd) # DON'T MESS UP THE +/- SIGNS
	rv2.append(rvraw2[i] + bcvstd - bcv[i] + rvstd)
	rv1_err.append(rvraw1_err[i])
	rv2_err.append(rvraw2_err[i])
	print ('%.9f %.9f %.9f %.5f %.5f %.5f %.5f' % (bjdmid[i], phase[i], bjdfunny[i], 
			rv1[i], rv1_err[i], rv2[i], rv2_err[i]), file=g2)
g1.close()
g2.close()
print('*******')
print(' ')
print('DID YOU FOLLOW THE INSTRUCTIONS ABOVE??? if so, then...')
print('BJD, orbital phase, and RVs with errors written to %s.' % outfile)
print('Use rvplotmaker.py to plot the RV curve.')

# PLOT THE FINAL SMOOTHED BFS + GAUSSIAN FITS IN INDIVIDUAL PANELS
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Uncorrected Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=26)
fig.text(0.07, 0.5, 'Broadening Function', ha='center', va='center', size=26, rotation='vertical')
for i in range (1,nspec):
	#ax = fig.add_subplot(6,4,i)
	ax = fig.add_subplot(2,2,i)
	ax.yaxis.set_major_locator(MultipleLocator(0.2))
	#if i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21:
	#	ax.set_yticklabels(())
	#if i!=20 and i!=21 and i!=22 and i!=23:
	#	ax.set_xticklabels(())
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.axis([-125, 125, -0.1, 0.55])
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.text(-115, 0.45, '%s' % (datetimelist[i].iso[0:10]), size=12)
	plt.plot(bf_ind, bfsmoothlist[i], color='k', lw=1.5, label='Smoothed BF')
	plt.plot(bf_ind, bffitlist[i][1], color='#e34a33', lw=1.5, label='Two-Gaussian fit')
#	if source[i] == 'arces': plt.text(-115, 0.25, 'ARCES', color='#0571b0', size=12)
#	if source[i] == 'tres': plt.text(-115, 0.25, 'TRES', color = '#008837', size=12)
	plt.text(-115, 0.35, '%.3f $\phi$' % (phase[i]), size=12)
	if i==23: ax.legend(bbox_to_anchor=(2.1,0.7), loc=1, borderaxespad=0., 
						frameon=False, prop={'size':20})
plt.show()