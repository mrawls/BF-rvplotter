import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import h5py
from echelle_funcs import *
'''
Reads in an ARCES spectrum that has been reduced but is still 2D (multiorder) FITS
Uses info in WAT2_* FITS headers to compute a wavelength solution for each order: *raw spectrum*
Attempts to measure the continuum and divide by it: *flattened spectrum*
--> this result is very sensitive to the parameters you feed continuumFlattenSpec!!
Saves the *raw* and *flattened* spectra (wavelength, flux) arrays to a new HDF5 file
--> this file is formatted so Starfish can read in the raw spectra! (still 2D/multiorder)
(wtf is Starfish: http://iancze.github.io/Starfish/current/overview.html)
'''


fitsfile = '../../../RG_spectra/tequila_201604/kplr7037405/140422.0024.ec.fits' # works
#fitsfile = '../../../KIC_8848288/tyc3559/TYC3559-02080-1.0024.ec.fits' # doesn't work
hdf5out = '../../../Starfish/7037405_Starfish/test.hdf5'
#hdf5out = '../../../KIC_8848288/test.hdf5'

fitsspec = fits.open(fitsfile)
fitsfluxes = fitsspec[0].data
norders = len(fitsspec[0].data)
nfluxpts = len(fitsspec[0].data[0])
header = fitsspec[0].header
fitswaves = getOrderWavelengths(fitsspec)
TruncateOrders = False
SaveFlatFluxes = False

# Sort each order so it goes from lowest to highest wavelength
# Note that this only sorts values within each order (the order of the orders is unchanged...)
for idx, (waves, fluxes) in enumerate(zip(fitswaves, fitsfluxes)):
    fitsfluxes[idx] = fluxes[np.argsort(waves)]
    fitswaves[idx] = waves[np.argsort(waves)]

# Continuum flatten and plot the orders we care about
flatwavelist = []
flatfluxlist = []
orderstart = 11    # red end
orderstop = 85 #92 # blue end
allrawwaves = fitswaves[orderstart:orderstop]
allrawfluxes = fitsfluxes[orderstart:orderstop]
for waves, fluxes in zip(allrawwaves, allrawfluxes):
    flatwaves, flatfluxes = continuumFlattenSpec(waves, fluxes, window=50, fitplot=False)
    # option to cut the ends off each order
    if TruncateOrders == True:
        trunc = int(len(flatwaves)/5)
        flatwavelist.append(flatwaves[trunc:-trunc])
        flatfluxlist.append(flatfluxes[trunc:-trunc])
        plt.plot(flatwaves[trunc:-trunc], flatfluxes[trunc:-trunc])
    else:
        flatwavelist.append(flatwaves)
        flatfluxlist.append(flatfluxes)
        plt.plot(flatwaves, flatfluxes)
    plt.axhline(y=1, ls=':', color='k')
    plt.axhline(y=0, ls=':', color='k')
plt.show()

# Save the raw and flattened spectral data in an HDF5 file
allflatwaves = np.array([np.array(item) for item in flatwavelist])
allflatfluxes = np.array([np.array(item) for item in flatfluxlist])

# Create 'sigmas' arrays for the spectrum (a wild guess at error bars)
allrawsigmas = []
allflatsigmas = []
for fluxes in allrawfluxes:
    rawsigmas = np.sqrt(np.abs(fluxes))
    allrawsigmas.append(rawsigmas)
for fluxes in allflatfluxes:
    flatsigmas = np.sqrt(np.abs(fluxes))
    allflatsigmas.append(flatsigmas)

hdf5spec = h5py.File(hdf5out, 'w')

rawflshdf5 = hdf5spec.create_dataset('fls', data=allrawfluxes)
rawsigmashdf5 = hdf5spec.create_dataset('sigmas', data=allrawsigmas)
rawwlshdf5 = hdf5spec.create_dataset('wls', data=allrawwaves)

# currently encountering a bug sometimes when trying to save flattened fluxes...
# "TypeError: Object dtype dtype('O') has no native HDF5 equivalent" ??? cute.
if SaveFlatFluxes == True:
    flatflshdf5 = hdf5spec.create_dataset('flatfls', data=allflatfluxes)
    flatsigmashdf5 = hdf5spec.create_dataset('flatsigmas', data=allflatsigmas)
    flatwlshdf5 = hdf5spec.create_dataset('flatwls', data=allflatwaves)
    print('Raw and flattened spectrum written to {0}'.format(hdf5out))
else:
    print('Raw spectrum only written to {0}'.format(hdf5out))

hdf5spec.close()

## Maybe in the future... interpolate onto a full 1D spectrum grid?
##
#newwaves = np.arange(4400, 5900, 0.15)
#newwaves = np.arange(3920, 9000, 0.15)
#allflatwaves = np.array(flatwavelist[::-1]).flatten()
#allflatfluxes = np.array(flatfluxlist[::-1]).flatten()
#allflatwaves = [item for sublist in allwaves for item in sublist]
#allflatfluxes = [item for sublist in allfluxes for item in sublist]
#print(len(allwaves), len(allfluxes), allwaves[0:10], allwaves[-10:-1])
#newfluxes = np.interp(newwaves, allwaves, allfluxes)
#plt.plot(newwaves, newfluxes)
#plt.show()