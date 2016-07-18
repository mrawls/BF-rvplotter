from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
#import MakeModel
from PyAstronomy import pyasl
import pandas as pd
import sys
sys.path.insert(0, '../../github/BF-rvplotter')
import BF_functions as bff
'''
Make a basic model of the atmospheric O2 A-band in the 7650 A region.
Write out this model at STP.
Fit this model to a FITS spectrum.
Based on the examples provided in the TelFit-1.2 package, available
here: http://www.as.utexas.edu/~kgulliks/projects.html
This is kind of hacked together and loosely based on BF_python.py...

Input: three text files - telfitin, timefile, and gausspar
(a list of FITS infiles, a list of JD timestamps, and guesses for velocity shifts)
The timefile contents don't matter one bit here

Output: two plots in succession (like BF_python.py), and RV shift info printed to screen
(manually save the shift info if you want, and then use specshift.py to remove the shift)
'''

telfitin =  '../../RG_spectra/9970396/infiles_arcesBF_telfit.txt'
timefile =  '../../RG_spectra/9970396/bjdinfile_arcesBF.txt'
gausspar =  '../../RG_spectra/9970396/gausspars_telfit.txt'
shiftfile = '../../RG_spectra/9970396/shifts_TEST.txt'

# Parameters for the broadening function (don't change w00 or n if you want the A-band!)
w00 = 7595
n = 1000
stepV = 1.5
mbf = 171

##### first, the option to make a basic model...
# once this file exists, you don't have to re-create it each time

temp = 273.15 # Kelvin
pres = 1013.25 # hPa (1 hectoPascal = 100 Pa = 1 mbar)
humid = 50 # percent relative humidity
o2val = 2.0e5
angle = 45.0
co2 = 368.50
o3 = 3.90/100
ch4 = 1.80
co = 1.40/10

# Set the start and end wavelength, **in nm**
# This is for the O2 A-band
#wavestart = 758
#waveend = 772

# Make an instance of the Modeler class
#modeler = MakeModel.Modeler()
#modeler.MakeModel(pressure=pres, temperature=temp, humidity=humid, o2=o2val,
#	save=True, libfile=None, lowfreq=1e7/waveend, highfreq=1e7/wavestart, resolution=31000)

# Manually access the model file
#modelfile = "transmission" + "-%.2f" % pres + "-%.2f" % temp + "-%.1f" % humid + "-%.1f" % angle + "-%.2f" % (co2) + "-%.2f" % (o3 * 100) + "-%.2f" % ch4 + "-%.2f" % (co * 10)
#wavemod, fluxmod = np.loadtxt(modelfile, unpack=True)
#wavemod = wavemod*10
#print(wavemod, fluxmod)

##### next, compare the model with a real spectrum and find the RV shift...

# set variables and read stuff in
isAPOGEE = False # a necessary bff.logify_spec parameter
period = 1.0 # nobody cares
BJD0 = 2450000.0 # nobody cares
rvstd = 0 # there's no standard star
bcvstd = 0 # there's no standard star
amp = 5.0 # arbitrary multiplication factor
w1, m, r = bff.logify_spec(isAPOGEE, w00, n, stepV, mbf)
nspec, filenamelist, datetimelist, wavelist, speclist, source = bff.read_specfiles(infiles=telfitin, bjdinfile=timefile)

# interpolate everything onto the same log-spaced grid
newspeclist = []
yoffset = 1
plt.axis([w1[0], w1[-1], 0, 27])
plt.xlabel('wavelength')
for i in range (0, nspec):
	newspec = np.interp(w1, wavelist[i], speclist[i])
	newspeclist.append(newspec)
	plt.plot(w1, newspeclist[i]+yoffset) #, label=datetimelist[i].iso[0:10])
	yoffset = yoffset + 1
#plt.legend()
plt.show()

# single value decomposition
svd = pyasl.SVD()
svd.decompose(newspeclist[0], m)
bflist = []
bfsmoothlist = []
for i in range (0, nspec):
	# Obtain the broadening function
	bf = svd.getBroadeningFunction(newspeclist[i]) # this one is like a matrix
	bfarray = svd.getBroadeningFunction(newspeclist[i], asarray=True)
	# Smooth the array-like broadening function
	bfsmooth = amp*pd.rolling_window(bfarray, window=5, win_type='gaussian', std=1.5, center=True)
	# The rolling window makes nans at the start because it's a punk.
	for j in range(0,len(bfsmooth)):
		if np.isnan(bfsmooth[j]) == True:
			bfsmooth[j] = 0
	bflist.append(bf)
	bfsmoothlist.append(bfsmooth)
# Obtain the indices in RV space that correspond to the BF
bf_ind = svd.getRVAxis(r, 1)

# plot the smoothed BFs
# this plot is boring, skip it
#plt.axis([-100, 70, -0.2, 12])
#plt.xlabel('Velocity (km s$^{-1}$)')
#plt.ylabel('Broadening Function (arbitrary amplitude)')
#yoffset = 0.0
#for i in range(1, nspec):
#	plt.plot(bf_ind, bfsmoothlist[i]+yoffset, color='b')
#	yoffset = yoffset + 0.4
#plt.show()

# fit the smoothed BFs with two gaussians (only one is really used here)
bffitlist = bff.gaussparty(gausspar, nspec, filenamelist, bfsmoothlist, bf_ind)
rvraw1 = []; rvraw1_err = []
rvraw1.append(0), rvraw1_err.append(0)
for i in range(1, len(bffitlist)):
    rvraw1.append(bffitlist[i][0][1])
    rvraw1_err.append(bffitlist[i][2][1])

# write results of the first fit gaussian only to an outfile
outfile = open(shiftfile, 'w')
for (image, rv, err) in zip(filenamelist[1:], rvraw1[1:], rvraw1_err[1:]):
    print(image, rv, err, file=outfile)
outfile.close()
print('Shifts saved to {0}'.format(shiftfile))

# plot final BFs in individual panels
windowcols = 4		                        # how many window columns there should be
windowrows = np.rint(nspec/windowcols)+1	# how many window rows there should be
xmin = -19
xmax = 19
ymin = -0.05
ymax = 0.7
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Velocity (km s$^{-1}$)', ha='center', va='center', size=26)
fig.text(0.07, 0.5, 'Broadening Function', ha='center', va='center', size=26, rotation='vertical')
for i in range (1,nspec):
	ax = fig.add_subplot(windowrows, windowcols,i)
	ax.yaxis.set_major_locator(MultipleLocator(0.2))
	if i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21 and i!=25:
		ax.set_yticklabels(())
	if i!=20 and i!=21 and i!=22 and i!=23 and i!=24 and i!=25:
		ax.set_xticklabels(())
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.axis([xmin, xmax, ymin, ymax])
	plt.tick_params(axis='both', which='major', labelsize=14)
	plt.plot(bf_ind, bfsmoothlist[i], color='k', lw=1.5, label='Smoothed BF')
	plt.plot(bf_ind, bffitlist[i][1], color='#e34a33', lw=1.5, label='Two-Gaussian fit')	
	if i==23: ax.legend(bbox_to_anchor=(2.1,0.7), loc=1, borderaxespad=0., 
						frameon=False, prop={'size':20})
plt.show()