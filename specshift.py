import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
'''
Meredith Rawls, 2015

Shift a set of spectrum by a list of constant velocities.
Based partly on 'telluric.py', also by Meredith Rawls, but that one was crappy.
(I wrote that a couple years ago, which is why some of the loop syntax is special.)

Input: list of FITS spectra and shifts to be applied
Output: creates new FITS files with shifted spectra in the working directory
'''

shiftfile = 'telfit_RVshifts_really.txt'
#shiftfile = '../../../Dropbox/KIC9246715/rvstandards/telfit_shifttest.txt'
infiles, velshiftlist = np.loadtxt(shiftfile, dtype={'names': ('infiles', 'velshiftlist'),
	'formats': ('|S77', np.float64)}, usecols=(0,1), unpack=True)

# Manually define a 1D wavelength array that will be used for interpolation
wavestart_new = 3850.
dwave_new = 0.0455
wavelen = 107000
waveref = np.arange(wavelen)*dwave_new + wavestart_new
print('Wavelength range: {0} - {1}'.format(waveref[0], waveref[-1]))
c = 2.99792e18 #angstroms per sec

# Loop over each FITS spectrum, read things in, define wavelength & pixel scales
nspec = len(infiles)
speclist = [] # these will be lists of arrays
wavelist = []
#pixellist = []
datetimelist = []
print ' '
print 'Filename, original wave start, original wave step, velshift being applied:'
for infile, velshift in zip(infiles, velshiftlist):
	hdu = fits.open(infile)
	spec = hdu[0].data
	head = hdu[0].header
	datetime = head['date-obs']
	datetimelist.append(datetime)
	# Define the wavelength scale (we need an x-axis for the data)
	dwave = head['cdelt1']
	wavestart = head['crval1']
	wave = np.arange(len(spec))*dwave + wavestart
	# This is here because the wave array was sometimes 1 longer than it should be
	if len(wave) != len(spec):
		minlength = min(len(wave), len(spec))
		wave = wave[0:minlength]
		spec = spec[0:minlength]
	wavelist.append(wave)
	speclist.append(spec)
	#print infile[-24:], wavestart, dwave, velshift # USE THIS TO CHECK THE ORIGINAL WAVELENGTH SCALE!

realnewwave = np.empty((nspec,wavelen*2), dtype=np.float64)
for i in range (0, nspec):
	for j in range (0, len(wavelist[i])):
		realnewwave[i,j] = wavelist[i][j] * velshiftlist[i] / c + wavelist[i][j]
		# the shifted wavelengths! (a 2D array) ... but it's not linear yet

# Make new spectra that goes with a linear interpolation of realnewwave back onto waveref.
realnewwavesinglelist = []
brandnewspeclist = []
for i in range (0, nspec):
	realnewwavesingle = np.empty(len(wavelist[i]), dtype=np.float64)
	for j in range (0, len(wavelist[i])):
		realnewwavesingle[j] = realnewwave[i,j]
	realnewwavesinglelist.append(realnewwavesingle)
	brandnewspec = np.interp(waveref, realnewwavesingle, speclist[i])
	brandnewspeclist.append(brandnewspec) #goes with waveref

print 'Writing new FITS files... '
print ' '

# Write the final wavelengths and associated spectral data to a new FITS file.
# Create an appropriate header.
headernote = 'Modified w/ telluric shift'
for i in range (0, nspec):
	outfile = 's_' + infiles[i][-24:]
	# Read in the original FITS header
	hdu = fits.open(infiles[i])
	head = hdu[0].header
	# Make a new FITS file with the new data and old header
	hdu = brandnewspeclist[i]
	head['cdelt1'] = (dwave_new, headernote)
	head['crval1'] = (wavestart_new, headernote)
	head['cd1_1'] = (dwave_new, headernote)
	fits.writeto(outfile, hdu, header=head, clobber=True, output_verify='fix')#, output_verify='warn')
print ' '
print 'New FITS files written in working directory. Filenames begin with \'s_\'.'
print ' '