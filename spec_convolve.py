import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from astropy.io import fits as pyfits
import sys
from subprocess import call
from pyraf import iraf
from pyraf.iraf import onedspec
'''
Convolve a Phoenix model stellar spectrum to a new resolution using the fortran CONVSPEC.
Originally by Jean McKeever
Modified by Meredith Rawls 7/2015

TO RUN:
python thisprogram.py mymodel.dat

OR, if you want near-IR APOGEE wavelengths returned instead:
python thisprogram.py mymodel.dat apogee

You need to have pyraf/iraf up and running with a legit login.cl in the working directory
You need to have the fortran convspec program compiled in a directory specified below

OUTPUT:
hopefully a model spectrum text file that ends in '.bf.arces.txt' or '.bf.apogee.txt'
'''

# thanks Jean
no='no'
yes='yes'

# Tell python where to look for convspec and the associated parameter file
#convspecexe = '/home/tequila/jeanm12/convspec/convspec'
convspecexe = '../jeanshare/convspec/convspec'
convspecpar = '../jeanshare/convspec.par'

def read_fits(fits):
   '''
   Simple function to read in a FITS file
   '''
   hdu = pyfits.open(fits)
   data = hdu[0].data
   hdr = hdu[0].header
   ll = np.arange(len(data))*hdr['cdelt1']+hdr['crval1']     
   return ll, data


"""
def nrefrac(wavelength, density=1.0):
   #Calculate refractive index of air from Cauchy formula.
   #nput: wavelength in Angstrom, density of air in amagat (relative to STP,
   #e.g. ~10% decrease per 1000m above sea level).
   #Returns N = (n-1) * 1.e6. 
   # The IAU standard for conversion from air to vacuum wavelengths is given
   # in Morton (1991, ApJS, 77, 119). For vacuum wavelengths (VAC) in
   # Angstroms, convert to air wavelength (AIR) via: 
   #  AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
   wl=wavelength
   wl2inv = (1.e4/wl)**2
   refracstp = 272.643 + 1.2288 * wl2inv  + 3.555e-2 * wl2inv**2
   return density * refracstp
"""


"""
only need two columns from lte files
awk '{print $1, $2, $3}' lte045-2.5-0.0a+0.0.BT-Settl.7 | sed 's/\D/e/g' |sort -g -k 1 > model_rg1.dat
awk '{print $1, $2, $3}' lte045-3.0-0.0a+0.0.BT-Settl.7 | sed 's/\D/e/g' |sort -g -k 1 > model_rg2.dat
awk '{print $1, $2, $3}' lte065-4.0-0.0a+0.0.BT-Settl.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms1.dat
awk '{print $1, $2, $3}' lte065-4.5-0.0a+0.0.BT-Settl.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms2.dat


awk '{print $1, $2, $3}' lte048-2.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_rg48.dat

awk '{print $1, $2, $3}' lte063-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms63.dat
awk '{print $1, $2, $3}' lte066-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms66.dat
awk '{print $1, $2, $3}' lte070-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms70.dat
awk '{print $1, $2, $3}' lte067-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms67.dat
awk '{print $1, $2, $3}' lte065-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms65.dat
awk '{print $1, $2, $3}' lte058-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms58.dat
awk '{print $1, $2, $3}' lte055-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms55.dat
awk '{print $1, $2, $3}' lte042-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms42.dat
awk '{print $1, $2, $3}' lte062-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms62.dat
awk '{print $1, $2, $3}' lte028-4.5-0.0a+0.0.BT-Settl.spec.7 | sed 's/\D/e/g' |sort -g -k 1 > model_ms28.dat
"""


# Original spectrum model is entered when you run the program
try:
   model = sys.argv[1]
except:
   print('You need to specify the original model at runtime... python THISPROGRAM.PY mymodel.dat')
   sys.exit(1)

# Make sure we select the appropriate wavelength range   
try:
   if sys.argv[2]=='apogee':
      arces=False
except:
   arces=True

# Read in the model file, and set output wavelength ranges
wl, f, b = np.loadtxt(model, unpack=True)
if arces:
   ind2 = np.where((wl>4000)&(wl<9000))[0]
   print('arces')
else:
   ind2 = np.where((wl>15000)&(wl<17000))[0]
   print('apogee')

# Define stuff that makes sense for Phoenix model spectra
wl = wl[ind2]
f = f[ind2]
DF = -8.0
f = 10**(f+DF)

if arces:
# need to put things in air wavelengths and not leave them in vacuum
   fit = interp.interp1d(wl, f)
   l = np.arange(wl[0], wl[-1], .02)
   newf = fit(l)
   VAC = l
   AIR = VAC / (1.0 +  5.792105e-2/(238.0185 - (1.e4/VAC)**2) + 1.67917e-3/( 57.362 - (1.e4/VAC)**2)) #FROM APOGEE
   bfmod = model[:-4]+'.bf.arces.txt'
else:
# leave things in vacuum wavelengths
   AIR = wl
   newf = f
   bfmod = model[:-4]+'.bf.apogee.txt'

# Various outfiles, some of which will be deleted at the end
textcont = model[:-4]+'_short.txt'
fitscont = model[:-4]+'_short.fits'
fitsnocont = model[:-4]+'_nocont.fits'
convtext = model[:-4]+'.cnv'
convout = model[:-4]+'.cnv.out'

# Continuum fit the spectrum so we can use IRAF's continuum function to remove it
np.savetxt(textcont, np.array([AIR,newf]).transpose())
iraf.rspectext(input=textcont, output=fitscont, flux='No', dtype='interp')
print('iraf rspectext complete')
iraf.continuum(input=fitscont, output=fitsnocont, lines='*',
               type='ratio', replace=no, wavescale=yes, logscale=no,
               listonly=no, interactive=no, sample='*',
               naverage=-25, function='spline3', order=5, niterate=10,
               low_reject=1.5, high_reject=5.0, markrej=no, grow=0,
               override=yes)
print('iraf continuum complete')
print('friendly untouched by convspec FITS saved {0}'.format(fitsnocont))

# Create a text file correctly formatted for convspec
wave, data = read_fits(fitsnocont)
np.savetxt(convtext, np.array([wave, np.zeros(len(wave)), data]).transpose())

# Run convspec, the fortran code
call([convspecexe, convtext, convspecpar])
print('completed convspec')

# Read the output from convspec
a, b = np.loadtxt(convout, unpack=True, usecols=[0,2])
np.savetxt(bfmod, np.array([a,b]).transpose())
print('convspec-modified file saved {0}'.format(bfmod))
print('deleting intermediate files')


call(['rm',textcont])
call(['rm',fitscont])
#call(['rm',fitsnocont])
call(['rm',convtext])
call(['rm',convout])

