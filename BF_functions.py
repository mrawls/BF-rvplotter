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
'''

def read_specfiles(infiles = 'infiles_BF.txt', isAPOGEE = False):
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
			print('You have a text file. The date will be WRONG without a FITS header.')
			# treat it like a text file
			filenamelist.append(infile)
			datetimelist.append(Time('2015-01-01 00:00:00', scale='utc', format='iso')) # placeholder
			wave, spec = np.loadtxt(open(infile), comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
			if isAPOGEE == True: # we need to normalize it
				spec = spec / np.median(spec)
			wavelist.append(wave)
			speclist.append(spec)
		else:
			# assume it's a FITS file
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
			else: # non-APOGEE (linear): create wavelength values from header data
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
	return nspec, filenamelist, datetimelist, wavelist, speclist, source