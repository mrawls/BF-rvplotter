from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
'''
My first real foray into proper object oriented python
SciCoder2 Workshop, Jan 2016
'''

#from BF_functions_better import Arces1dSpec
from BF_functions_better import ListOfSpectra

dir = '../../RG_spectra/7037405_1/'
fitsfile = 's_lspec150506.0033.ec.fits'

#myspec = Arces1dSpec(dir+fitsfile)
#print(myspec.dateobs)

infile = dir+'infiles_arcesBF_updated.txt'
speclist = ListOfSpectra(infile)

specfiles = speclist.specfiles
waves = speclist.waves
fluxes = speclist.fluxes

for file, wave, flux in zip(specfiles, waves, fluxes):
    print(file[-20:], wave[0], flux[0])