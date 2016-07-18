import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
'''
Meredith Rawls, 2015

Shift a set of spectrum by a list of constant velocities.
Result is New_Wavelength = Old_Wavelength - Shift_Value.
(Think of it like... the spectra are off by X, so you need to remove the stupid X)

Input: text file containing a list of FITS spectra and shifts to be applied
Output: creates new FITS files with shifted spectra **in the working directory**

NOTE: this is how Jean shifts spectra instead. It's, uh, easier, but uses IRAF.

## Shift the spectra by rv1 amount... using iraf dopcor ##
for i in range(1,len(filenamelist)):
    iraf.dopcor(filenamelist[i],'shifted_'+filenamelist[i][-27:],rvraw1[i],
		isvelocity=yes,add=no,dispersion=yes,flux=no,factor=3.0,verbose=yes)
## Shift the spectra by rv1 amount... using iraf dopcor ##
'''

#shiftfile = '../../RG_spectra/10001167/shifts_arcesBF_telfit.txt'
shiftfile = '../../RG_spectra/8702921/shifts_zeroRV.txt'

infiles, velshiftlist = np.loadtxt(shiftfile, dtype={'names': ('infiles', 'velshiftlist'),
	'formats': ('|S77', np.float64)}, usecols=(0,1), unpack=True)

# Loop over each FITS spectrum, read things in, define wavelength & pixel scales
nspec = len(infiles)
speclist = [] # these will be lists of arrays
wavelist = []
#pixellist = []
datetimelist = []
print ' '
print 'Filename, original wave start, original wave step, velshift being applied:'
dwave_list = []
for infile, velshift in zip(infiles, velshiftlist):
	hdu = fits.open(infile)
	spec = hdu[0].data
	head = hdu[0].header
	datetime = head['date-obs']
	datetimelist.append(datetime)
	# Define the wavelength scale (we need an x-axis for the data)
	dwave = head['cdelt1']
	wavestart = head['crval1']
	dwave_list.append(dwave)
	wave = np.arange(len(spec))*dwave + wavestart
	# This is here because the wave array was sometimes 1 longer than it should be
	if len(wave) != len(spec):
		minlength = min(len(wave), len(spec))
		wave = wave[0:minlength]
		spec = spec[0:minlength]
	wavelist.append(wave)
	speclist.append(spec)
	print infile[-24:], wavestart, dwave, velshift # USE THIS TO CHECK THE ORIGINAL WAVELENGTH SCALE!

# Manually define a 1D wavelength array that will be used for interpolation
wavestart_new = 3850.
#dwave_new = 0.0440 #0.0455
dwave_new = np.max(dwave_list)
wavelen = 120000
waveref = np.arange(wavelen)*dwave_new + wavestart_new
print('Wavelength range: {0} - {1}'.format(waveref[0], waveref[-1]))
print('Wavelength step size (CDELT1): {0}'.format(dwave_new))
c = 2.99792e5 #km per sec OMG NOT ANGSTROMS YOU MORON

# Shift the wavelengths by the required amount
# Make new spectra to go with said wavelengths
# Interpolate these onto waveref
brandnewspeclist = []
for idx, (wavelengths, specdata) in enumerate(zip(wavelist, speclist)):
    realnewwave = []
    for item in wavelengths:
        realnewwave.append(item - (item * (velshiftlist[idx] / c)))
    #print(wavelengths - realnewwave) #definitely shouldn't be zero
    brandnewspeclist.append(np.interp(waveref, realnewwave, specdata))

print 'Writing new FITS files... '
print ' '

# Write the final wavelengths and associated spectral data to a new FITS file.
# Create an appropriate header.
#headernote = 'Modified w/ telluric shift'
headernote = 'Shifted to zero RV'
for idx, specdata in enumerate(brandnewspeclist):
	#outfile = 's_' + infiles[idx][-24:]
	outfile = 'rest_' + infiles[idx][-26:]
	# Read in the original FITS header
	hdu = fits.open(infiles[idx])
	head = hdu[0].header
	# Make a new FITS file with the new data and old header
	head['cdelt1'] = (dwave_new, headernote)
	head['crval1'] = (wavestart_new, headernote)
	head['cd1_1'] = (dwave_new, headernote)
	fits.writeto(outfile, specdata, header=head, clobber=True, output_verify='ignore')#output_verify='fix')#, output_verify='warn')
print ' '
print 'New FITS files written in working directory.'
print ' '