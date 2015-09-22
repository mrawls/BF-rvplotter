import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
#from astropy.io.fits import getdata
#from astropy.io.fits import getheader
from BF_functions import read_specfiles
'''
Simple python script for plotting spectra from FITS files
Meredith Rawls, originally written May 2014

To run: python splot.py infiles.txt
'''

# CHANGE AXIS RANGE AS DESIRED
#plt.axis([7590, 7720, -0.5, 25])
#plt.axis([6100, 6300, -0.5, 25])
#plt.axis([5000, 8000, -0.5, 10])

import sys
infiles = sys.argv[1]
#bjdfile = '../../FDBinary/9246715/bjds_baryvels.txt'
bjdfile = '../../RG_spectra/7037405/bjdinfile_arces.txt'
	# note: bcvin is never really used here, but you have to give it some file
	# with the correct length because of how read_specfiles is set up
#linelist = '../../MOOG/ares_v1.0/yong_fe2_jeanversion.txt'
#refspec = '../../../Dropbox/KIC9246715/phoenix_air_jean.txt'
#refspec = '../../RG_spectra/APOGEE/model_rg_apogee.txt' OPTION to plot this too, see end
isAPOGEE = False # read_specfiles will check to see if apogee is actually true
#if 'apogee' in infiles:
#	isAPOGEE = True
#	print('Setting isAPOGEE to True')
#else:
#	isAPOGEE = False

# Option: draw vertical lines for comparison of feature location by eye
# (the vertical lines will be every 2 Angstroms between xstart and xend)
#xstart = 7610
#xend = 7680
#for i in range(xstart,xend):
#	plt.axvline(x=i, ymin=-1, ymax=10, color='0.75')

plt.xlabel("Wavelength")
plt.ylabel("Arbitrary Flux")

# Read in a text file containing a list of fits files
#infiles = "infiles.txt"

i = 0.0
yoffset = 0.5

### This only worked for FITS infiles in a standard format
### I replaced it with a more generic call to read_specfiles
# for line in open(infiles):
# 	infile = line.rstrip()
# 	# Read in the FITS file with all the data in the primary HDU
# 	hdu = fits.open(infile)
# 	spec = hdu[0].data + 2*i
# 	head = hdu[0].header
# 	try:
# 		datetime = head['date-obs']
# 	except:
# 		datetime = head['date']
# 	exptime = head['exptime']
# 
# 	# Define the wavelength scale (we need an x-axis for the data)
# 	dwave = head['cdelt1']
# 	wavestart = head['crval1']
# 	wavestop = wavestart + dwave*len(spec)
# 	wave = np.arange(wavestart, wavestop, dwave)
# 	# This is here because the wave array was sometimes 1 longer than it should be?!
# 	if len(wave) != len(spec):
# 		minlength = min(len(wave), len(spec))
# 		wave = wave[0:minlength]
# 		spec = spec[0:minlength]	
	
nspec, filenamelist, datetimelist, wavelist, speclist, source = read_specfiles(infiles, bjdfile, isAPOGEE)
for i, (wave, spec, datetime) in enumerate(zip(wavelist, speclist, datetimelist)):
	# Add this data to the plot
	plt.plot(wave, spec + 2*i, color='b', marker='.', mfc='k')#
	i = i + yoffset

# OPTION: plot text file reference spectrum here
#waveref, fluxref = np.loadtxt(refspec, usecols=(0,1), unpack=True)
#fluxref = fluxref + 2*i #apply yoffset
#plt.plot(waveref, fluxref, color='r', marker='.', mfc='k')

# Plot random vertical lines here
#xvals = np.loadtxt(linelist, comments='#', usecols=(0,), unpack=True)
#for xval in xvals:
#	plt.axvline(x=xval, color='0.75')

plt.show()