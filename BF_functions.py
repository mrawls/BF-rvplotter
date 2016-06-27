from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl
from scipy import ndimage
import pandas as pd
import gaussfitter as gf
'''
Functions used in BF_python.py
Read the damn comments
(I'm sorry there aren't more objects)
'''

def logify_spec(isAPOGEE=False, w00=5400, n=38750, stepV=1.7, m=171):
    # The new log-wavelength array will be w1. it will have equal spacing in velocity.
    # Specify reasonable values when you call this function or else bad things will happen.
### GUIDELINES FOR CHOOSING GOOD INPUT VALUES ###
#    good APOGEE values
#        w00 = 15145         # starting wavelength of the log-wave array in Angstroms
#        n = 20000           # desired length of the log-wave vector in pixels (must be EVEN)
#    good ARCES values
#        w00 = 5400          # starting wavelength of the log-wave array in Angstroms
#        n = 38750           # desired length of the log-wave vector in pixels (must be EVEN) 
#    stepV = 1.7             # step in velocities in the wavelength vector w1
#    m = 171                 # length of BF (must be ODD)
### GUIDELINES FOR CHOOSING GOOD INPUT VALUES ###
    r = stepV/299792.458     # put stepV in km/s/pix
    w1 = w00 * np.power((1+r), np.arange(float(n)))
    print('The new log-wavelength scale will span %d - %d A with stepsize %f km/s.' % (w1[0], w1[-1], stepV))
    print(' ')
    return w1, m, r

def read_one_specfile(infile = 'myspectrum.txt', isAPOGEE = False):
    '''
    Read in a single FITS or txt spectrum file
    (Bare-bones version of read_specfiles, below)
    Requires    infile, isAPOGEE
    Returns     wave, spec
    '''
    if infile[-3:] == 'txt':
        try:
            wave, spec = np.loadtxt(open(infile), comments='#', usecols=(0,1), unpack=True)
            print('Text file {0}, isAPOGEE = {1}'.format(infile[-15:], isAPOGEE))
        except:
            raise FileNotFoundError('The file {0} was not found or cannot be opened'.format(infile))
        if isAPOGEE == True: # we need to sort by wavelength
            spec = spec[np.argsort(wave)]
            wave = wave[np.argsort(wave)]
    elif infile[-4:] == 'fits' or infile[-4:] == 'FITS':
        # assume it's a FITS file
        # Read in the FITS file with all the data in the primary HDU
        try:
            hdu = fits.open(infile)
        except:
            print('{0} not found or cannot be opened'.format(infile))
        else:
            head = hdu[0].header
            try: datetime = head['date-obs']
            except: datetime = head['date']
            print('FITS file {0}, isAPOGEE = {1}, header date {2}'.format(infile[-17:], isAPOGEE, datetime))                            
        if isAPOGEE == True: # APOGEE: the data is in a funny place and backwards
            spec = hdu[1].data
            spec = spec.flatten()
            spec = spec[::-1]
        else: # non-APOGEE (regular) option
            spec = hdu[0].data
        # Define the original wavelength scale
        if isAPOGEE == True: # APOGEE: read wavelength values straight from FITS file
            wave = hdu[4].data
            wave = wave.flatten()
            wave = wave[::-1]
        else: # non-APOGEE (linear): create wavelength values from header data
            headerdwave = head['cdelt1']
            headerwavestart = head['crval1']
            headerwavestop = headerwavestart + headerdwave*len(spec)
            wave = np.arange(headerwavestart, headerwavestop, headerdwave)
        if len(wave) != len(spec): # The wave array is sometimes 1 longer than it should be?
            minlength = min(len(wave), len(spec))
            wave = wave[0:minlength]
            spec = spec[0:minlength]
        try: # check to see if we have a file with log angstroms
            logcheck = head['dispunit'] 
        except:
            logcheck = 'linear' # assume linear if no 'dispunit' is in header
        if logcheck == 'log angstroms':
            wave = np.power(10,wave) # make it linear
    else:
        print('File does not end in \'txt\' or \'fits\', no spectrum loaded.')
        wave = []; spec = []
    return wave, spec
    
def read_specfiles(infiles = 'infiles_BF.txt', bjdinfile = 'bjds_baryvels.txt', isAPOGEE = False):
    '''
    Read in some FITS or TXT files that are spectra and may or may not be APOGEE
    Requires     infiles, bjdinfile, isAPOGEE
    Returns     nspec, filenamelist, datetimelist, wavelist, speclist
    '''
    f1 = open(infiles)
    print('Reading files from 1st column in %s' % infiles)
    #print('The first one had better be your template spectrum.')
    print(' ')
    speclist = []; wavelist = []
    filenamelist = []; datetimelist = []
    if isAPOGEE == False:
        checkAPOGEE = True #notallinfiles are APOGEE, but let's check in case *some* are
    else:
        checkAPOGEE = False #all the infiles are APOGEE so we don't have to search
    i = 0
    for line in f1: # This loop happens once for each spectrum
        infile = line.rstrip()
        if checkAPOGEE == True: # check to see if a subset of infiles are from APOGEE or not
            if 'apogee' in infile or 'APOGEE' in infile: isAPOGEE = True
            else: isAPOGEE = False
        if infile[-4:] == 'fits' or infile[-4:] == 'FITS':
            # assume it's a FITS file
            try:
                hdu = fits.open(infile)
                head = hdu[0].header
                filenamelist.append(infile)
                try: datetime = head['date-obs']
                except: datetime = head['date']
                datetimelist.append(Time(datetime, scale='utc', format='isot'))
                print('FITS file {0}, isAPOGEE = {1}, header date {2}'.format(infile[-17:], isAPOGEE, datetime))
            except:
                raise FileNotFoundError('The file {0} was not found or cannot be opened'.format(filename))
            # it's time to dig out the spectral (flux) data and the wavelength scale!
            if isAPOGEE == True: # APOGEE: the data is in a funny place and backwards
                wave, spec = ProcessAPOGEEFITS(hdu)
            else: # not APOGEE
                spec = hdu[0].data # hope the info we want is in the zeroth HDU
                try:
                    headerdwave = head['cdelt1']
                    headerwavestart = head['crval1']
                    headerwavestop = headerwavestart + headerdwave*len(spec)
                    wave = np.arange(headerwavestart, headerwavestop, headerdwave)
                except:
                    raise RuntimeError('Cannot find wavelength info in FITS header')
            if len(wave) != len(spec): # the wave array is sometimes 1 longer than it should be?
                minlength = min(len(wave), len(spec))
                wave = wave[0:minlength]
                spec = spec[0:minlength]
            try: # check to see if we have a file with log angstroms
                logcheck = head['dispunit'] 
            except:
                logcheck = 'linear' # hopefully, at least
            if logcheck == 'log angstroms':
                wave = np.power(10, wave) # make it linear
                #spec = spec / np.median(spec) # WARNING really basic, possibly bad normalization
        else: # treat it like a text file
            filenamelist.append(infile)
            datetime = np.loadtxt(bjdinfile, comments='#', usecols=(1,), unpack=True)[i]
            datetimelist.append(Time(datetime, scale='utc', format='jd'))
            try:
                wave, spec = np.loadtxt(open(infile), comments='#', usecols=(0,1), unpack=True)
                print('Text file {0}, isAPOGEE = {1}, bjdinfile date {2}'.format(infile[-17:], isAPOGEE, datetime))
            except:
                raise FileNotFoundError('The file {0} was not found or cannot be opened'.format(filename))
            if isAPOGEE == True: # we need sort by wavelength, just in case it hasn't been
                spec = spec[np.argsort(wave)]
                wave = wave[np.argsort(wave)]
#            if infile[0:5] == 'trans': # you have a model telluric spectrum in nm, not A
#                print("Assuming this is a telluric spectrum in nm, not A, proceed with caution")
#                wave = wave*10
        # at the end of this mess, we have one file's WAVE and corresponding SPEC - save it!
        wavelist.append(wave)
        speclist.append(spec)
        i = i + 1    
    # save the total number of spectra
    nspec = i
    f1.close()
    return nspec, filenamelist, datetimelist, wavelist, speclist

def ProcessAPOGEEFITS(hdu):
    '''
    Turns an APOGEE FITS hdu into a pair of wavelength and spectrum ndarrays
    '''
    spec = hdu[1].data
    spec = spec.flatten()
    spec = spec[::-1]
    spec = spec / np.median(spec) # WARNING really basic, possibly bad normalization
    wave = hdu[4].data
    wave = wave.flatten()
    wave = wave[::-1]
    return wave, spec

def gaussparty(gausspars, nspec, filenamelist, bfsmoothlist, bf_ind, threshold=10):
    '''
    Fits 2 or 3 gaussians to some data
    '''
    param = []
    with open(gausspars) as f1:
        for line in f1:
            if line[0] != '#':
                param.append( line.rstrip() )
    #param = np.loadtxt(gausspars, comments='#')
    bffitlist = []
    bffitlist.append(0)
    gauss1 = [[] for i in range(nspec)]
    gauss2 = [[] for i in range(nspec)]
    gauss3 = [[] for i in range(nspec)]
    gauss1[0] = [0,0]
    gauss2[0] = [0,0]
    gauss3[0] = [0,0]
    error_array = np.ones(len(bfsmoothlist[0]))*0.01 # dummy array with 0.01 error values
    print(' ')
    print('Gaussian fit results: peak amplitude, width, rvraw, rvraw_err')
    print ('-------------------------------------------------------------')
    for i in range(1, nspec):
        # check to see if we are fitting a third gaussian, i.e., one near zero
        # don't print out the result of this fit, but do return it for plotting
        # handle comments in gausspars file without exploding
        if '#' in param[i]:
            commentbegin = param[i].find('#')
            partest = param[i][0:commentbegin].split()
        else:
            partest = param[i].split()
        if len(partest) == 6: ngauss = 2
        elif len(partest) == 9: ngauss = 3
        else: print('something is wrong with your gausspars file!')       
        minpars=[0.002, float(partest[1])-threshold, 0,  0.002, float(partest[4])-threshold, 0]
        maxpars=[1.2,   float(partest[1])+threshold, 15, 1.2,   float(partest[4])+threshold, 15]
        if ngauss == 2:
            bffit = gf.multigaussfit(bf_ind, bfsmoothlist[i], ngauss=ngauss, 
                    params=partest, err=error_array,
                    limitedmin=[True,True,True], limitedmax=[True,True,True], 
                    minpars=minpars, maxpars=maxpars, quiet=True, shh=True)
        elif ngauss == 3:
            minpars.extend([0.002, float(partest[7])-threshold, 0])
            maxpars.extend([1.2,   float(partest[7])+threshold, 15])
            bffit = gf.multigaussfit(bf_ind, bfsmoothlist[i], ngauss=ngauss, 
                    params=partest, err=error_array,
                    limitedmin=[True,True,True], limitedmax=[True,True,True], 
                    minpars=minpars, maxpars=maxpars, quiet=True, shh=True)
        newbffit = [[] for x in range(len(bffit))]
        # Sometimes bffit[2] is None, or contains None. Set it to zeros instead.
        try:
            if not any(bffit[2]): # this will fail if bffit[2] = None
                newbffit[0] = bffit[0]
                newbffit[1] = bffit[1]
                newbffit[2] = [0, 0, 0, 0, 0]
            else:
                newbffit = bffit
        except:
            print('WARNING - gaussfit is acting up, fit may not have happened, adjust gausspars file:')
            if not bffit[2]: # this catches the case where bffit[2] = None
                newbffit[0] = bffit[0]
                newbffit[1] = bffit[1]
                newbffit[2] = [0, 0, 0, 0, 0]
            else:
                newbffit = bffit
        bffitlist.append(newbffit)
        # NOTE: to get the gaussian fit corresponding to bfsmoothlist[i], use bffitlist[i][1].
        # RV1 for observation i is bffitlist[i][0][1] +/- bffitlist[i][2][1].
        # RV2 for observation i is bffitlist[i][0][4] +/- bffitlist[i][2][4].
        # (note: need to check if bffit[2] == None before calling bffit[2][1] or bffit[2][4])
        print('{0:s}    {1:.3f} {2:.2f} {3:.4f} {4:.4f} \t {5:.3f} {6:.2f} {7:.4f} {8:.4f}'.format(
            filenamelist[i][-20:], newbffit[0][0], newbffit[0][2], newbffit[0][1], newbffit[2][1],
            newbffit[0][3], newbffit[0][5], newbffit[0][4], newbffit[2][4]))
    print(' ')
    print('You MUST manually guesstimate the location of each Gaussian\'s peak in %s!' % gausspars)
    print('Until you do, the above values will be WRONG and the plot will look TERRIBLE.')
    print(' ')
    return bffitlist

def rvphasecalc(bjdinfile, outfile, nspec, period, BJD0, rvraw1, rvraw1_err, rvraw2, rvraw2_err, rvstd, bcvstd):
    rv1 = []; rv2 = []
    rv1.append(0); rv2.append(0)
    rv1_err = []; rv2_err = []
    rv1_err.append(0); rv2_err.append(0)
    g1 = open(bjdinfile)
    g2 = open(outfile, 'w')
    print('Calculating RVs...')
    bjdmid, bcv = np.loadtxt(g1, comments='#', usecols=(1,2), unpack=True)
    bjdfunny = bjdmid - 2454833.
    phase = []
    phase.append(0)
    for i in range(1, nspec):
        fracP = (bjdmid[i] - BJD0) / period
        if fracP < 0:
            phase.append(1 + (fracP % 1))
            cycle = int(fracP) - 1
        else:
            phase.append((fracP % 1))
            cycle = int(fracP)
        rv1.append(rvraw1[i] + bcv[i] - rvstd - bcvstd) # DON'T MESS UP THE +/- SIGNS
        rv2.append(rvraw2[i] + bcv[i] - rvstd - bcvstd)
        rv1_err.append(rvraw1_err[i])
        rv2_err.append(rvraw2_err[i])
        print ('%.9f %.9f %.9f %.5f %.5f %.5f %.5f' % (bjdmid[i], phase[i], bjdfunny[i], 
                rv1[i], rv1_err[i], rv2[i], rv2_err[i]), file=g2)
    g1.close()
    g2.close()
    print(' ')
    print('BJD, phase, and RVs written to %s.' % outfile)
    print('Use rvplotmaker.py to plot the RV curve.')
    return phase, bjdfunny, rv1, rv2, rv1_err, rv2_err
