import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
'''
Plot a model telluric O2 A-band on top of actual spectra.
'''

# Read in and plot the telluric model
telfit = 'transmission-1013.25-273.15-50.0-45.0-368.50-3.90-1.80-1.40.txt'
wavetel, modeltel = np.loadtxt(telfit, unpack=True) # wavetel is in nm by default
plt.xlabel('Wavelength')
plt.ylabel('Arbitrary Flux')
plt.axis([7580, 7700, -1, 25])
plt.plot(wavetel*10, modeltel, color='r')

# Read in a text file containing a list of fits files
infiles = "infiles_telfit.txt"
f1 = open(infiles)

i = 0.5
yoffset = 0.5
for line in f1:
	infile = line.rstrip()

	# Read in the FITS file with all the data in the primary HDU
	hdu = fits.open(infile)
	spec = hdu[0].data + 2*i
	head = hdu[0].header
	# this is less efficient than the above
	#spec = getdata(infile) + i
	#head = getheader(infile)
	datetime = head['date-obs']
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

	# Define the pixel scale (if you want to plot vs. pixel instead of vs. wavelength)
	#dpixel = 1
	#pixelstart = head['crpix1']
	#pixelstop = pixelstart + dpixel*len(spec)
	#pixel = np.arange(pixelstart, pixelstop, dpixel)
	
	#print infile
	#print len(wave)
	#print len(pixel)
	#print len(spec)	
	
	# Add this data to the plot
	plt.plot(wave, spec, color='b', label=datetime[0:10])
	i = i + yoffset
	
	#print (infile, datetime, exptime)

# Plot random lines here
plt.axvline(x=7676.55, color='0.75')

plt.show()
