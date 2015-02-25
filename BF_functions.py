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
Functions used in BF_python.py
Read the damn comments
'''

def logify_spec(isAPOGEE = False):
	# The new log-wavelength array will be w1. it will have equal spacing in velocity.
	# please specify reasonable values below or else bad things will happen.
	# note log = log base 10 because SERIOUSLY, come on now.
	if isAPOGEE == True:
		w00 = 15145		# starting wavelength of the log-wave array in Angstroms
		n = 20000			# desired length of the log-wave vector in pixels (must be EVEN)
	else:
		w00 = 5400 #5055 		# starting wavelength of the log-wave array in Angstroms
		n = 38750 #7000 		# desired length of the log-wave vector in pixels (must be EVEN) 
	stepV = 1.7 #2.7 			# step in velocities in the wavelength vector w1
	m = 171 				# length of BF (must be ODD)
	r = stepV/299792.458 	# put stepV in km/s/pix
	w1 = w00 * np.power((1+r), np.arange(float(n)))
	print ('The new log-wavelength scale will span %d - %d A with stepsize %f km/s.' % (w1[0], w1[-1], stepV))
	print(' ')
	return w1, m, r

def read_specfiles(infiles = 'infiles_BF.txt', bjdinfile = 'bjds_baryvels.txt', isAPOGEE = False):
	f1 = open(infiles)
	print('Reading files from 1st column in %s' % infiles)
	print('The first one had better be your template spectrum.')
	print(' ')
	speclist = []; wavelist = []
	filenamelist = []; datetimelist = []
	source = [] # option: keep track of which spectrograph was used (ARCES vs. TRES)
	i = 0
	for line in f1: # This loop happens once for each spectrum
		infile = line.rstrip()
		if infile[-3:] == 'txt':
			print('You have a text file. Reading BJD date from bjdinfile, not FITS header.')
			# treat it like a text file
			filenamelist.append(infile)
			datetimelist.append( Time(np.loadtxt(open(bjdinfile), comments='#', dtype=np.float64, usecols=(1,), unpack=True)[i], scale='utc', format='jd') )
			#datetimelist.append(Time('2015-01-01 00:00:00', scale='utc', format='iso')) # placeholder
			wave, spec = np.loadtxt(open(infile), comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
			if isAPOGEE == True: # we need to normalize it and sort by wavelength
				source.append('apogee')
				spec = spec / np.median(spec)
				spec = spec[np.argsort(wave)]
				wave = wave[np.argsort(wave)]
#				# VACUUM CORRECTION
#				if i != 0:
#					wave = wave / (1 + 5.792105e-2/(238.0185 - 1/np.power(wave,2)) + 1.67917e-3/(57.362 - 1/np.power(wave,2)) )
			else: # unknown source
				source.append('unknown')
			wavelist.append(wave)
			speclist.append(spec)
		else:
			# assume it's a FITS file
		# Read in the FITS file with all the data in the primary HDU
			hdu = fits.open(line)
			if isAPOGEE == True: # APOGEE: the data is in a funny place, backwards, not normalized, and in VACUUM WAVELENGTHS !!
				source.append('apogee')
				spec = hdu[1].data ### APOGEE
				spec = spec.flatten() ### APOGEE
				spec = spec[::-1] ### APOGEE
				spec = spec / np.median(spec)
#				# VACUUM CORRECTION
#				if i != 0:
#					wave = wave / (1 + 5.792105e-2/(238.0185 - 1/np.power(wave,2)) + 1.67917e-3/(57.362 - 1/np.power(wave,2)) )
			else: # non-APOGEE (regular) option
				spec = hdu[0].data
			head = hdu[0].header
		#	# *** begin SPECIAL FOR MORE THAN ONE SPECTROGRAPH (ARCES + TRES) ONLY ***
			if '.tres.' in line: source.append('tres')
			elif '.ec.' in line: source.append('arces')
			else: source.append('unknown')
#			if head['imagetyp'] == 'object': source.append('arces')
#			if head['imagetyp'] == 'OBJECT': source.append('tres')
		#	# *** end SPECIAL FOR MORE THAN ONE SPECTROGRAPH (ARCES + TRES) ONLY ***
			filenamelist.append(infile)
			try:
				datetime = head['date-obs']
			except:
				datetime = head['date']
			datetimelist.append(Time(datetime, scale='utc', format='isot'))
			# Define the original wavelength scale
			if isAPOGEE == True: # APOGEE: read wavelength values straight from FITS file
				wave = hdu[4].data ### APOGEE
				wave = wave.flatten() ### APOGEE
				wave = wave[::-1] ### APOGEE
			else: # non-APOGEE (linear): create wavelength values from header data
				headerdwave = head['cdelt1']
				headerwavestart = head['crval1']
				headerwavestop = headerwavestart + headerdwave*len(spec)
				wave = np.arange(headerwavestart, headerwavestop, headerdwave)
			if len(wave) != len(spec): # The wave array is sometimes 1 longer than it should be?
				minlength = min(len(wave), len(spec))
				wave = wave[0:minlength]
				spec = spec[0:minlength]
			try: # check to see if we have a file with log angstroms
				logcheck = head['dispunit'] 
			except:
				logcheck = 'linear' # hopefully, at least
			if logcheck == 'log angstroms':
				wave = np.power(10,wave) # make it linear
				spec = spec / np.median(spec) # also normalize it to 1
			#print(wave, spec)
			wavelist.append(wave)
			speclist.append(spec)
		i = i + 1	
	# save the total number of spectra
	nspec = i
	f1.close()
	return nspec, filenamelist, datetimelist, wavelist, speclist, source

def gaussparty(gausspars, nspec, filenamelist, bfsmoothlist, bf_ind):
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
	print(' ')
	print('Two-Gaussian fit results: peak amplitude, width, rvraw, rvraw_err (for each star)')
	print ('---------------------------------------------------------------------------')
	for i in range(1, nspec):
		bffit = gf.multigaussfit(bf_ind, bfsmoothlist[i], ngauss=2, 
				params=param[i], err=error_array,
				limitedmin=[True,True,True], limitedmax=[True,True,True], 
				minpars=[0.2,-100,1], maxpars=[0.6,100,5], quiet=True, shh=True)
		bffitlist.append(bffit)
		# NOTE: to get the gaussian fit corresponding to bfsmoothlist[i], use bffitlist[i][1].
		gauss1[i] = [bffit[0][0], bffit[0][2], bffit[0][1], bffit[2][1]] # these are [amp1, width1, rvraw1, rvraw1_err]
		gauss2[i] = [bffit[0][3], bffit[0][5], bffit[0][4], bffit[2][4]] # these are [amp2, width2, rvraw2, rvraw2_err]
		print ('%s \t %.5f %.5f %.5f %.5f \t %.5f %.5f %.5f %.5f' % (filenamelist[i][-15:], 
			gauss1[i][0], gauss1[i][1], gauss1[i][2], gauss1[i][3], gauss2[i][0], gauss2[i][1], gauss2[i][2], gauss2[i][3]))
		rvraw1.append(bffit[0][1])
		rvraw2.append(bffit[0][4])
		rvraw1_err.append(bffit[2][1])
		rvraw2_err.append(bffit[2][4])
	print(' ')
	print('You MUST manually guesstimate the location of each Gaussian\'s peak in %s!' % gausspars)
	print('Until you do, the above values will be WRONG and the plot will look TERRIBLE.')
	print(' ')
	return bffitlist, rvraw1, rvraw1_err, rvraw2, rvraw2_err

def rvphasecalc(bjdinfile, outfile, nspec, period, BJD0, rvraw1, rvraw1_err, rvraw2, rvraw2_err, rvstd, bcvstd, source):
	rv1 = []; rv2 = []
	rv1.append(0); rv2.append(0)
	rv1_err = []; rv2_err = []
	rv1_err.append(0); rv2_err.append(0)
	g1 = open(bjdinfile)
	g2 = open(outfile, 'w')
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
		rv1.append(rvraw1[i] + bcv[i] - rvstd + bcvstd) # DON'T MESS UP THE +/- SIGNS
		rv2.append(rvraw2[i] + bcv[i] - rvstd + bcvstd)
		rv1_err.append(rvraw1_err[i])
		rv2_err.append(rvraw2_err[i])
		print ('%.9f %.9f %.9f %.5f %.5f %.5f %.5f %s' % (bjdmid[i], phase[i], bjdfunny[i], 
				rv1[i], rv1_err[i], rv2[i], rv2_err[i], source[i]), file=g2)
	g1.close()
	g2.close()
	print(' ')
	print('BJD, phase, and RVs written to %s.' % outfile)
	print('Use rvplotmaker.py to plot the RV curve.')
	return phase, bjdfunny, rv1, rv2, rv1_err, rv2_err