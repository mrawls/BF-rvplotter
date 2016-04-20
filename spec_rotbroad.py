from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
from BF_functions import read_one_specfile
'''
Simulate Doppler broadening on a model input spectrum.

Inputs: spectrum as FITS or txt, value of vsini, value of linear LD coefficient
Outputs: broadened spectrum (with a new evenly-spaced wavelength grid)
'''

# adjust these values as desired
dir = '../../RG_spectra/coelho2014_spectra/'

# 9291629
#infile = dir+'norm_t04750_g+3.0_p00p00_hrplc.fits'
#vsini = 19.394
#linearLD = 0.80
#outfile = dir+'norm_t04750_g+3.0_p00p00_rotbroad_v'+str(vsini)+'l'+str(linearLD)+'.txt'

# 3955867
infile = dir+'norm_t05000_g+3.0_m05p00_hrplc.fits'
vsini = 12.36
linearLD = 0.91
outfile = dir+'norm_t05000_g+3.0_m05p00_rotbroad_v'+str(vsini)+'l'+str(linearLD)+'.txt'

# read in spectrum
wave, spec = read_one_specfile(infile)

# check if spectrum is already evenly spaced in wavelength
# if it's not, resample it so it is
def all_same(items):
    return all(x == items[0] for x in items)
difflist = []
for idx, point in enumerate(wave[0:-1]):
    difflist.append(wave[idx+1] - point)
if all_same(difflist) == True:
    print('Wavelength grid already evenly spaced')
else:
    print('Resampling spectrum so it is evenly spaced in wavelength') 
    dwave = np.max(difflist)
    newwave = np.arange(wave[0], wave[-1], dwave)
    newspec = np.interp(newwave, wave, spec)
    wave = newwave
    spec = newspec

# apply rotational broadening
print('Applying rotational broadening, please be patient...')
bspec = pyasl.rotBroad(wave, spec, linearLD, vsini)

# save the result
fout = open(outfile, 'w')
for w, s in zip(wave, bspec):
    print(w, s, file=fout)

# plot the result
plt.title("Rotational broadening")
plt.xlabel("Wavelength")
plt.ylabel("Arbitrary flux")
plt.plot(wave, spec, 'k-')
plt.plot(wave, bspec, 'r-')
plt.show()