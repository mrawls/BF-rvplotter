def moving_average(series, window=100, sigma=50):
    '''
    Calculate the moving average of a series using scipy. Used by continuumFlattenSpec.
    Borrowed from http://www.nehalemlabs.net/prototype/blog/2014/04/12/how-to-fix-scipys-interpolating-spline-default-behavior/
    '''
    from scipy.signal import gaussian
    from scipy.ndimage import filters
    b = gaussian(window, sigma)
    average = filters.convolve1d(series, b/b.sum())
    var = filters.convolve1d(np.power(series-average,2), b/b.sum())
    return average, var

def getOrderWavelengths(hdu):
    '''
    Retrieve wavelength arrays for a multiorder echelle spectra.
    ONLY LINEAR, LOG, OR CHEBYSHEV POLYNOMIAL DISPERSIONS ARE SUPPORTED FOR NOW
    
    input: FITS hdu object 
    output: list of wavelength arrays, one per order
    
    some code is from https://github.com/kgullikson88/General/blob/master/readmultispec.py
    '''
    import numpy as np
    from astropy.io import fits
    header = hdu[0].header
    #print(list(header.keys()))

    # before things get interesting, try a normal linear dispersion
    try:
        dwave = header['cdelt1']
        wavestart = header['crval1']
        wavestop = headerwavestart + headerdwave*len(spec)
        wave = np.arange(wavestart, wavestop, dwave)
    except KeyError:
        pass
    
    # OK then, each order has a set of N=lenwat useful parameters hidden in the WAT2 header entries
    WAT2string = ''
    for key in header.keys():
        if 'WAT2_' in key: # build a giant string
            #print(key, header[key])
            if len(header[key]) < 68: # catch the entries that need a trailing space added (length 67)
                WAT2string += header[key] + ' '
            elif len(header[key]) > 68:
                raise ValueError('Length of header entry {0} is {1}; do you even FITS?!'.format(key, len(header[key])))
            else:
                WAT2string += header[key]
    #print(WAT2string)
    WAT2list = [] # use the giant string to create a list we can slice into orders
    lenwatlist = []
    for item in WAT2string.split("\""):
        if 'spec' not in item:
            lenwatlist.append(len(item.split(' ')))
            WAT2list.extend(item.split(' '))
    norders = len(hdu[0].data)
    lenwat = lenwatlist[1]
    WAT2orders = [] # placeholder list for info in each order
    for idx in np.arange(0, lenwat*norders, lenwat):
        WAT2orders.append(WAT2list[idx:idx+lenwat])
        #print(WAT2list[idx:idx+lenwat])
    fitswaves = []
    for idx, order in enumerate(WAT2orders):
        #print(order)
        dtype = order[2]
        w1 = np.float64(order[3])   # dispersion coordinate of 1st physical pixel
        dw = np.float64(order[4])   # average dispersion interval per physical pixel
        nwave = int(order[5])       # number of wavelength points in the order
        z = np.float64(order[6])    # Doppler factor
        apmin, apmax = float(order[7]), float(order[8])  # original pixel limits along the spatial axis, not used
        if dtype == '0':   # linear
            wavelengths = (w1 + dw * np.arange(nwave, dtype=np.float64)) / (1. + z)
        elif dtype == '1': # log
            wavelengths = (w1 + dw * np.arange(nwave, dtype=np.float64)) / (1. + z)
            wavelengths = np.power(10., wavelengths)
        elif dtype == '2': # nonlinear
            if np.float64(order[6]) != 0:
                print(np.float64(order[6]))
                raise Warning('Nonzero Doppler factor in order {0}, not accounting for this.'.format(idx))
            wt_i, w0_i, ftype_i = np.float64(order[9]), np.float64(order[10]), int(order[11])
            cheb_order = int(order[12])
            # ftype_i means 1 for a Chebyshev polynomial; 2-6 for other polynomials
            # ONLY CONSIDERING CHEBYSHEV FOR NOW!
            if ftype_i != 1:
                raise ValueError('Sorry, the nonlinear dispersion is not a Chebyshev polynomial.')           
            pmin = np.float64(order[13])
            pmax = np.float64(order[14])
            pmiddle = (pmax + pmin) / 2
            prange = pmax - pmin
            coeffs = []
            for cidx in range(0, cheb_order):
                coeffs.append(np.float64(order[-5+cidx])) 
            xs = (np.arange(nwave, dtype=np.float64) + 1 - pmiddle) / ((prange) / 2)        
            p0 = np.ones(nwave, dtype=np.float64)
            p1 = xs
            wavelengths = p0 * coeffs[0] + p1 * coeffs[1]
            for i in range(2, cheb_order):
                p2 = 2 * xs * p1 - p0
                wavelengths = wavelengths + p2 * coeffs[i]
                p0 = p1
                p1 = p2
            #print(wavelengths) # it works!
        else:
            raise ValueError('Spectrum type not recognized.')
        fitswaves.append(list(wavelengths))
    fitswaves = np.array(fitswaves, dtype=np.float64)
    return fitswaves # numpy array

def continuumFlattenSpec(waves, fluxes, window=50, fitplot=True):
    '''
    Fits a spline to a spectrum and divides to continuum flatten it.
    Returns waves and fluxes/continuum with length = original_length - 2*window.
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import UnivariateSpline
    import pandas as pd
    from PyAstronomy import pyasl
    
    # Decide how long the final arrays will be
    newlength = len(waves) - 2*window
    
    # Identify outlier points
    iin, iout = pyasl.slidingPolyResOutlier(waves, fluxes, points=window*2, deg=1, stdlim=3, mode='above', controlPlot=False)
    #print("Number of outliers: ", len(iout))
    #print("Indices of outliers: ", iout)
    if len(iout) > 3*window:
        raise ValueError('More than {0} outliers found, adjust stdlim'.format(3*window))
    # Remove outliers
    waves, fluxes = np.array(waves[iin]), np.array(fluxes[iin])
    
    # Fit the continuum
    specsmooth_top = pd.Series(fluxes).rolling(window=window, center=True).max()
    #specsmooth_bot = pd.Series(fluxes).rolling(window=window, win_type='gaussian', center=True).mean(std=window)
    specsmooth_bot, variance = np.array(moving_average(fluxes, window=2*window, sigma=window))
    variance = variance[window/2:-window/2]
    specsmooth = []
    for top, bot, flux in zip(np.array(specsmooth_top)[window/2:-window/2], np.array(specsmooth_bot)[window/2:-window/2], fluxes[window/2:-window/2]):
        specsmooth.append((top + bot) / 2.)
    waves = waves[window/2:-window/2]
    fluxes = fluxes[window/2:-window/2]
    # inspired by http://www.nehalemlabs.net/prototype/blog/2014/04/12/how-to-fix-scipys-interpolating-spline-default-behavior/
    continuum = UnivariateSpline(waves, specsmooth, k=5, w=1./np.sqrt(variance))(waves)

    # Plot the result
    if fitplot == True:
        fig = plt.figure()
        fig.add_subplot(211)
        plt.plot(waves, fluxes, 'b-', label='spectrum')
        plt.plot(waves, np.array(specsmooth_top)[window/2:-window/2], 'k:')
        plt.plot(waves, np.array(specsmooth_bot)[window/2:-window/2], 'k:')
        plt.plot(waves, specsmooth, 'k-', label='smoothed')
        plt.plot(waves, continuum, 'r-', label='continuum')
        fig.add_subplot(212)
        plt.plot(waves, fluxes/continuum, 'b-', label='flattened')
        plt.axhline(y=1, ls=':', color='k')
        plt.show()
    
    return waves[0:newlength], fluxes[0:newlength]/continuum[0:newlength] # numpy arrays