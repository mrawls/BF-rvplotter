from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Radial velocity plotter
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase
You need to have the text file 'rvs_out1.txt' in the working directory
'''

# File containing times, phases, and RVs as specified below
filename = 'rvs_out1.txt'

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2
timestart = 1117
timeend = 1970
phasemin = 0.5
phasemax = 1.5
RVmin = -59
RVmax = 59

f1 = open(filename)
bjd, phase, rv1, rverr1, rv2, rverr2 = np.loadtxt(f1, comments='#', 
	dtype=np.float64, usecols=(0,1,3,4,5,6), unpack=True)
f1.close()

# MANUAL ADJUSTMENT PARTY
# the last two APOGEE points are off by about 15 km/s for an unknown reason
rv1[-1] = rv1[-1] - 15
rv2[-1] = rv2[-1] - 15
rv1[-2] = rv1[-2] - 15
rv2[-2] = rv2[-2] - 15

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
plt.errorbar(bjd-2454833, rv1, yerr=rverr1, fmt='ko:', color='0.75', mfc=red, mec=red, ms=8, lw=1.5)
plt.errorbar(bjd-2454833, rv2, yerr=rverr2, fmt='ko:', color='0.75', mfc=yel, mec=yel, ms=8, lw=1.5)
plt.xlabel("Time (BJD -- 2454833)", size=24, labelpad=10)

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major', labelsize=20)
plt.errorbar(phase_double, rv1_double, yerr=rverr1_double, marker='o', color=red, mec=red, ecolor=red, ms=8, ls='None', lw=1.5, label='Star 1')
plt.errorbar(phase_double, rv2_double, yerr=rverr2_double, marker='o', color=yel, mec=yel, ecolor=yel, ms=8, ls='None', lw=1.5, label='Star 2')
plt.xlabel("Orbital Phase", size=24)

# Draw lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Create the legend and y-axis label
plt.legend(ncol=2, loc=1, fontsize=24, numpoints=1, frameon=False, bbox_to_anchor=(1,2.2))
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=24, rotation='vertical')

plt.show()