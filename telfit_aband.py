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
'''

##### first, the option to make a basic model...

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
isAPOGEE = False
period = 1.0
BJD0 = 2450000.0
rvstd = 0
bcvstd = 0
amp = 5.0
w1, m, r = bff.logify_spec(isAPOGEE, w00=7595, n=1000, stepV=1.7, m=171)
nspec, filenamelist, datetimelist, wavelist, speclist, source = bff.read_specfiles(infiles='infiles_telfit.txt', bjdinfile='times.txt')

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
	bfsmooth = amp*pd.rolling_window(bfarray, window=5, win_type='gaussian', std=1.5)
	# The rolling window makes nans at the start because it's a punk.
	for j in range(0,len(bfsmooth)):
		if np.isnan(bfsmooth[j]) == True:
			bfsmooth[j] = 0
	bflist.append(bf)
	bfsmoothlist.append(bfsmooth)
# Obtain the indices in RV space that correspond to the BF
bf_ind = svd.getRVAxis(r, 1)

# plot the smoothed BFs
plt.axis([-100, 70, -0.2, 12])
plt.xlabel('Velocity (km s$^{-1}$)')
plt.ylabel('Broadening Function (arbitrary amplitude)')
yoffset = 0.0
for i in range(1, nspec):
	plt.plot(bf_ind, bfsmoothlist[i]+yoffset, color='b')
	yoffset = yoffset + 0.4
plt.show()

# fit the smoothed BFs with two gaussians (only one is really used here)
bffitlist, rvraw1, rvraw1_err, rvraw2, rvraw2_err = bff.gaussparty('gaussfit_pars.txt', nspec, filenamelist, bfsmoothlist, bf_ind)

# plot final BFs in individual panels
# manually adjust this multi-panel plot based on how many spectra you have
windowcols = 4		# how many window columns there should be
windowrows = 6		# how many window rows there should be
xmin = -79
xmax = 79
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
	#plt.text(0.9*xmin, 0.8*ymax, '%.3f $\phi$' % (phase[i]), size=12)
	#plt.text(0.9*xmin, 0.6*ymax, '%s' % (datetimelist[i].iso[0:10]), size=12)
	if source[i] == 'arces': plt.text(0.4*xmax, 0.8*ymax, 'ARCES', color='#0571b0', size=12)
	elif source[i] == 'tres': plt.text(0.4*xmax, 0.8*ymax, 'TRES', color = '#008837', size=12)
	elif source[i] == 'arces': plt.txt(0.4*xmax, 0.8*ymax, 'APOGEE', color = 'k', size=12)
	else: plt.text(0.9*xmin, 0.4*ymax, 'SOURCE?', color = 'k', size=12)
	plt.plot(bf_ind, bfsmoothlist[i], color='k', lw=1.5, label='Smoothed BF')
	plt.plot(bf_ind, bffitlist[i][1], color='#e34a33', lw=1.5, label='Two-Gaussian fit')	
	if i==23: ax.legend(bbox_to_anchor=(2.1,0.7), loc=1, borderaxespad=0., 
						frameon=False, prop={'size':20})
plt.show()