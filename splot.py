import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
#from astropy.io.fits import getdata
#from astropy.io.fits import getheader
'''
Simple python script for plotting spectra from FITS files
Meredith Rawls, originally written May 2014
'''

# CHANGE AXIS RANGE AS DESIRED
#plt.axis([7590, 7720, -0.5, 25])
#plt.axis([6100, 6300, -0.5, 25])
plt.axis([5000, 8000, -0.5, 10])

# Option: draw vertical lines for comparison of feature location by eye
# (the vertical lines will be every 2 Angstroms between xstart and xend)
#xstart = 7610
#xend = 7680
#for i in range(xstart,xend):
#	plt.axvline(x=i, ymin=-1, ymax=10, color='0.75')

plt.xlabel("Wavelength")
plt.ylabel("Arbitrary Flux")

# Read in a text file containing a list of fits files
infiles = "infiles.txt"

i = 0.0
yoffset = 0.5
for line in open(infiles):
	infile = line.rstrip()
	# Read in the FITS file with all the data in the primary HDU
	hdu = fits.open(infile)
	spec = hdu[0].data + 2*i
	head = hdu[0].header
	try:
		datetime = head['date-obs']
	except:
		datetime = head['date']
	exptime = head['exptime']

	# Define the wavelength scale (we need an x-axis for the data)
	dwave = head['cdelt1']
	wavestart = head['crval1']
	wavestop = wavestart + dwave*len(spec)
	wave = np.arange(wavestart, wavestop, dwave)
	# This is here because the wave array was sometimes 1 longer than it should be?!
	if len(wave) != len(spec):
		minlength = min(len(wave), len(spec))
		wave = wave[0:minlength]
		spec = spec[0:minlength]
	
	#print infile
	#print len(wave)
	#print len(pixel)
	#print len(spec)	
	
	# Add this data to the plot
	plt.plot(wave, spec, color='b', marker='.', mfc='k', label=datetime[0:10])
	i = i + yoffset
	
	#print (infile, datetime)

# Plot random vertical lines here
#plt.axvline(x=7676.55, color='0.75')

#plt.legend()
plt.show()