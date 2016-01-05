from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
#from PyAstronomy import pyasl
#from scipy import ndimage
#import pandas as pd
#import gaussfitter as gf
'''
Beginning of a better suite of object-oriented functions to replace BF_functions.py
Use with specload_test.py

Work in progress, Jan 2016
'''

class Arces1dSpec(object):
    '''
    retrieves time of an ARCES spectral observation from header
    '''
    def __init__(self, filename = None):
        if filename == None:
            raise FileNotFoundError('The file {0} was not found'.format(filename))
        self.filename = filename
        self._dateobs = None
    
    @property
    def dateobs(self):
        if self._dateobs == None: # use private var to avoid reopening if dateobs already set
            hdu = fits.open(self.filename)
            head = hdu[0].header
            datetime = head['date-obs']
            self._dateobs = Time(datetime, scale='utc', format='isot')
        return self._dateobs
    


class ListOfSpectra(object):
    '''
    takes a text file listing FITS and/or TXT spectra; retrieves wavelengths and fluxes
    '''
    def __init__(self, filename = None):
        if filename == None:
            raise FileNotFoundError('The file {0} was not found'.format(filename))
        self.filename = filename
        self._specfiles = None
        self._waves = None
        self._fluxes = None
    
    @property
    def specfiles(self):
        if self._specfiles == None: # hasn't been loaded
            self._specfiles = self.loadspec()[0] # this loads it
        return self._specfiles
    
    @property
    def waves(self):
        if self._waves == None: # hasn't been loaded
            self._waves = self.loadspec()[1] # this loads it
        return self._waves

    @property
    def fluxes(self):
        if self._fluxes == None: # hasn't been loaded
            self._fluxes = self.loadspec()[2] # this loads it
        return self._fluxes
    
    def loadspec(self):
        specfiles = []; waves = []; fluxes = []
        f1 = open(self.filename)
        for line in f1:
            specfiles.append( line.rstrip() )
        for file in specfiles:
            if file[-3:] == 'txt':
                wave, spec = np.loadtxt(open(file), comments='#', usecols=(0,1), unpack=True)
            elif file[-4:] == 'fits':
                hdu = fits.open(file)
                head = hdu[0].header
                spec = hdu[0].data
                headerdwave = head['cdelt1']
                headerwavestart = head['crval1']
                headerwavestop = headerwavestart + headerdwave*len(spec)
                wave = np.arange(headerwavestart, headerwavestop, headerdwave)
                #dateobs = Arces1dSpec.dateobs
            else:
                print('File type not recognized')
                wave = []
                spec = []
            waves.append(wave)
            fluxes.append(spec)
        f1.close()      
        return specfiles, waves, fluxes
