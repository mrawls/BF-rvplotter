#from __future__ import print_function # good idea for python < 3
#from __future__ import division # good idea for python < 3
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
'''
reads in a PHOENIX model spectrum (FITS format)
makes a continuum flattened (normalized) spectrum over a specified wavelength range
saves the result to a text file

'specfile' is the FITS HiRes PHOENIX spectrum file which lives in directory 'dir'
'outtxt' is the new text file that will be written
'idxstart' and 'idxend' set the index of the wavelengths that are written out
(you can set the indices to 0 and -1, but that's a much bigger wavelength range than APOGEE)

you can get PHOENIX spectra from here: http://phoenix.astro.physik.uni-goettingen.de
(click "version 2.0 of the spectral library" and login as a guest when prompted)
'''

def moving_average(series, window=100, sigma=50):
    '''
    Calculates the moving average of a series by using a Gaussian kernel to smooth it.
    This function is used by continuumFlattenSpec.
    Borrowed from http://www.nehalemlabs.net/prototype/blog/2014/04/12/how-to-fix-scipys-interpolating-spline-default-behavior/

    Input: series or array, window size, and Gaussian width (sigma).
    Returns: the smoothed data and the smoothed variance.
    '''
    from scipy.signal import gaussian
    from scipy.ndimage import filters
    b = gaussian(window, sigma)
    average = filters.convolve1d(series, b/b.sum())
    var = filters.convolve1d(np.power(series-average,2), b/b.sum())    
    return average, var

def continuumFlattenSpec(waves, fluxes, win=50, varwin=60, varsig=100, weight=1.0, fitplot=True):
    '''
    Fits a spline to a spectrum and divides to continuum flatten it.
    Weighting inspired by http://www.nehalemlabs.net/prototype/blog/2014/04/12/how-to-fix-scipys-interpolating-spline-default-behavior/
    
    !!WARNING!! this function is not one-size-fits-all and the input parameters will
    probably need to be significantly tweaked for your specific case! It is VERY EASY
    to accidentally overfit or underfit your spectrum and get a poor continuum fit!
    
    More specifically, this function traces the shape of the spectrum in two ways:
    (1) the maximum values of a rolling window function with size 'win'
        --> this traces the TOP of the spectrum and may be an overestimate if it's spikey
    (2) a Gaussian-smoothed version of the spectrum with window size 'varwin' and Gaussian width 'varsig'
        --> this traces the spectrum's shape but is systematically LOWER than the continuum
        --> the variance of the spectrum is also saved during this step
    
    Then, scipy's UnivariateSpline function is used to fit a 5th order spline to either
    (1) OR a combination of (1) and (2). In both cases, the fit is weighted by
    by 'weight' divided by the square root of /sqrt(variance).
    This kind of weighting helps balance overfitting vs. underfitting, as described 
    at the link above, but it is not perfect.
    
    Finally, a plot of the result is shown (assuming 'fitplot' is True) so the user can
    decide if the continuum fit is acceptable or not.
    
    Input: wavelength array, flux array, and the parameters described above
    Returns: a pair of numpy arrays corresponding to wavelengths and flattened fluxes.
    Each array has length = original_length - 2*window.
    '''
    from scipy.interpolate import UnivariateSpline
    import pandas as pd
    from PyAstronomy import pyasl    
    # Fit the continuum
    window = win
    specsmooth_top = pd.Series(fluxes).rolling(window=window, center=True).max()
    specsmooth_bot, variance = np.array(moving_average(fluxes, window=varwin, sigma=varsig))
    variance = variance[window/2:-window/2]
    w = weight/np.sqrt(variance)
    specsmooth = []
    for top in np.array(specsmooth_top)[window/2:-window/2]:
        specsmooth.append(top) # BETTER CHOICE IF SPECTRUM IS NOT SPIKEY
        #specsmooth.append((top + bot) / 2.) # BETTER CHOICE IF SPECTRUM IS SPIKEY
    waves = waves[window/2:-window/2]
    fluxes = fluxes[window/2:-window/2]
    continuum = UnivariateSpline(waves, specsmooth, k=5, w=w)(waves)
    # Plot the result
    if fitplot == True:
        fig = plt.figure()
        fig.add_subplot(211)
        plt.plot(waves, fluxes, 'b-', label='spectrum')
        plt.plot(waves, np.array(specsmooth_top)[window/2:-window/2], 'k:')
        plt.plot(waves, np.array(specsmooth_bot)[window/2:-window/2], 'k:')
        plt.plot(waves, specsmooth, 'k-', label='smoothed')
        plt.plot(waves, continuum, 'r-', label='continuum')
        plt.legend()
        fig.add_subplot(212)
        plt.plot(waves, fluxes/continuum, 'b-', label='flattened')
        plt.axhline(y=1, ls=':', color='k')
        plt.legend()
        plt.show()    
    return waves, fluxes/continuum # numpy arrays


### MAIN PROGRAM BEGINS HERE ###

#######################
## edit values below ##
#######################
# directory where wavelength file lives:
wavedir = '../../../PHOENIX/PHOENIX-ACES-AGSS-COND-2011/'
# directory where spectrum file lives:
specdir = wavedir + 'Z-0.0/'
# spectrum file you want to run the program on:
specfile = 'lte05500-2.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
# wavelength region that will be fit with a spline:
wavestart = 14000; waveend = 18000 #wavestart = 3700; waveend = 9000
# wavelength region that will be saved (must be no longer than the above):
truncstart = 15000; truncend = 17000 #truncstart = 4400; truncend = 5800
# spline fitting parameters (final result is VERY SENSITIVE to these!!):
splinewindow = 100
variancewindow = 60
variancesigma = 150
splineweight = 0.9
#######################
## edit values above ##
#######################

# define full input and output file paths
wavefits = wavedir + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
infits = specdir + specfile
outtxt = specfile[:-5] + '_' + str(truncstart) + '-' + str(truncend) + '-norm.txt'

# read in PHOENIX FITS spectrum and corresponding wavelength FITS file
hdu = fits.open(infits)
spec = hdu[0].data
hduwave = fits.open(wavefits)
wave = hduwave[0].data

# put the wavelength array in goddamn air units
# reference: Huswer et al. 2013
wave = wave / (1.0 + 0.05792105/(238.0185-(1e4/wave)**2) + 0.00167917/(57.362-(1e4/wave)**2))

# truncate spectrum to the range you want to fit
idxstart = np.where(wave > wavestart)[0][0]
idxend = np.where(wave > waveend)[0][0]
wavechunk = wave[idxstart:idxend]
specchunk = spec[idxstart:idxend]

# measure the continuum level and divide by it to flatten the spectrum
print('Flattening spectrum...')
wavenorm, specnorm = continuumFlattenSpec(wavechunk, specchunk, win=splinewindow,
                                            varwin=variancewindow, varsig=variancesigma, 
                                            weight=splineweight, fitplot=True)

# truncate the flattened spectrum to final desired wavelength range
idxtruncstart = np.where(wavenorm > truncstart-100)[0][0]
idxtruncend = np.where(wavenorm > truncend+100)[0][0]
wavenorm = wavenorm[idxtruncstart:idxtruncend]
specnorm = specnorm[idxtruncstart:idxtruncend]

# plot the final spectrum
plt.plot(wavenorm, specnorm, ls='None', marker='.', label='final flattened spectrum')
plt.legend()
plt.show()

# write the result to a text file
f = open(outtxt, 'w')
for wentry, sentry in zip(wavenorm, specnorm):
    print(wentry, sentry, file=f)
f.close()
print('New spectrum written to {0}'.format(outtxt))