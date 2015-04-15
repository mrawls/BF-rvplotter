from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
'''
Radial velocity plotter
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase
You need to have the text file 'rvs_out1.txt' in the working directory

Columns should be as follows
BJD, PHASE, KEPLER_BJD, RV1, RV1_ERR, RV2, RV2_ERR, SOURCE

all of these are numbers, except 'source' is text

Update Dec 2014... this is getting awfully specific for my particular binary.
Now there are three different symbols for three different telescope/instrument combos.
'''

doShift = False

# File containing times, phases, and RVs as specified below
#filename = '../../TelFit/9246715_telfit/rvs_out_STPshift.txt'
filename = '../../../Dropbox/KIC9246715/rvs_out_STPshift_smoothfix.txt'
#filename = 'rvs_out1.txt'
#filename = '../../RG_spectra/rvs_out_arces.txt'
#filename = '../../RG_spectra/APOGEE/KIC7037405/rvs_out_test.txt'
#filename = '../../RG_spectra/APOGEE/KIC7037405/rvs_out_model.txt'
#filename = '../../RG_spectra/rvs_out_arces.txt'
sysname = '9246715'
#sysname = '7037405'

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2
timestart = 1117 #1600
timeend = 1970 #2300
phasemin = 0.5
phasemax = 1.5
RVmin = -70 #-80
RVmax = 50 #0

f1 = open(filename)
bjd, phase, rv1, rverr1, rv2, rverr2, source = np.loadtxt(f1, comments='#', 
	dtype={'names': ('bjd', 'phase', 'rv1', 'rverr1', 'rv2', 'rverr2', 'source'),
	'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, '|S15')},
	usecols=(0,1,3,4,5,6,7), unpack=True)
f1.close()

if doShift == True:
	try:
		shift = np.loadtxt('../../TelFit/9246715_telfit/telfit_RVshifts.txt', comments='#', usecols=(1,), unpack=True)
		rv1 = rv1 - shift
		rv2 = rv2 - shift
		print('Shift applied')
	except:
		print('Could not find shift file')

# Double the arrays so we can plot from phase 0 to phase 2 for clarity, if desired
rv1_double = np.concatenate((rv1,rv1), axis=0)
rv2_double = np.concatenate((rv2,rv2), axis=0)
phase_double = np.concatenate((np.array(phase),np.array(phase)+1.0), axis=0)
rverr1_double = np.concatenate((rverr1,rverr1), axis=0)
rverr2_double = np.concatenate((rverr2,rverr2), axis=0)

# Set up the figure
fig = plt.figure(1, figsize=(15,10))

# Unfolded RV vs time (BJD-2454833)
ax2 = plt.subplot(2,1,1)
plt.axis([timestart, timeend, RVmin, RVmax])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major', labelsize=20)
# dotted lines to guide the eye
plt.plot(bjd-2454833, rv1, color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd-2454833, rv2, color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
for idx, label in enumerate(source):
	if label == 'arces':
		plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='ko:', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
		plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='ko:', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
	elif label == 'tres':
		plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='kv:', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
		plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='kv:', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
	elif label == 'apogee':
		plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='ks:', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
		plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='ks:', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)				
##plt.errorbar(bjd-2454833, rv1, yerr=rverr1, fmt='ko:', color='0.75', mfc=red, mec=red, ms=10, lw=1.5)
##plt.errorbar(bjd-2454833, rv2, yerr=rverr2, fmt='ko:', color='0.75', mfc=yel, mec=yel, ms=10, lw=1.5)
plt.xlabel("Time (BJD -- 2454833)", size=24, labelpad=10)

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major', labelsize=20)
for idx, label in enumerate(source):
	idx2 = idx + len(source)
	if label == 'arces':
		plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='ARCES, Star 1' if idx==6 or idx==0 else '')
		plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5, label='Star 2' if idx==6 else '')
		plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
		plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
	elif label == 'tres':
		plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='v', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='TRES, Star 1' if idx==0 else '')
		plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='v', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5, label='Star 2' if idx==0 else '')
		plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='v', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
		plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='v', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
	elif label == 'apogee':
		plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='s', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='APOGEE, Star 1' if idx==23 else '')
		plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='s', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5, label='Star 2' if idx==23 else '')		
		plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='s', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
		plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='s', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
plt.xlabel("Orbital Phase", size=24)

# Draw lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Create the legend and y-axis label
plt.legend(ncol=3, loc=1, fontsize=20, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35), columnspacing=0.7)
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=24, rotation='vertical')
fig.text(0.14, 0.115, 'Folded', size=24)
fig.text(0.14, 0.55, 'Unfolded', size=24)
fig.text(0.2, 0.9, sysname, size=32)

plt.show()