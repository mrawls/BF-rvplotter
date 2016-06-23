from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
'''
Radial velocity plotter!
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase

Columns should be as follows
BJD, PHASE, KEPLER_BJD, RV1, RV1_ERR, RV2, RV2_ERR

Update September 2015:
Has two flag options
1. apply some shift to the RVs before plotting them
2. read in another set of calculated RVs from cols 8,9,10,11 for comparison and plot 
    RV1 = RV_col8 - RV_col3, RV2 = RV_col10 - RV_col5 instead, plus a line at RV = 0.
    
Update June 2016:
Simplified some options; no longer manually sets point shape as a function of "source" string.
If you want that functionality, use an older version of this code... it was messy.
NOTE that any RV value with an error bar = 0 is not plotted!
'''

doShift = False
compareRVs = False
dateoffset = 2454833. # this value will be subtracted from bjds in pane vs. time

### Input EB parameters and RV filenames ###

sysname = 'Joni OB'; filename = '../../Joni_EBs/OBrvoutfile.txt'
timestart = 1550; timeend = 2000
phasemin = 0.48; phasemax = 1.48
RVmin = -50; RVmax = 50

#sysname = '9970396'; filename = '../../RG_spectra/9970396/rvs_patrick.txt'
#timestart = 1550; timeend = 2500
#phasemin = 0.48; phasemax = 1.48
#RVmin = -59; RVmax = 19

#sysname = '3955867'; filename = '../../RG_spectra/3955867/rvs_final.txt'
#timestart = 1650; timeend = 2350
#phasemin = 0.5; phasemax = 1.5
#RVmin = -49; RVmax = 69

#sysname = '8702921'; filename = '../../RG_spectra/8702921/rvs_final.txt'
#timestart = 1200; timeend = 2200
#phasemin = 0.5; phasemax = 1.5
#RVmin = -29; RVmax = 9

#sysname = '9291629'; filename = '../../RG_spectra/9291629/rvs_final.txt'
#timestart = 1650; timeend = 2350
#phasemin = 0.5; phasemax = 1.5
#RVmin = -99; RVmax = 34

# for 8848288 (not one of our RG/EBs)
#sysname = '8848288'; filename = '../../KIC_8848288/rvs_BF.txt'
#timestart = 1310; timeend = 1350
#phasemin = 0.5; phasemax = 1.5
#RVmin = -24; RVmax = -13

#sysname = '5786154'; filename = '../../RG_spectra/5786154_1/rvs_patrick.txt'
#timestart = 1200; timeend = 2200
#phasemin = 0.5; phasemax = 1.5
#RVmin = -55; RVmax = 39

#sysname = '8430105'
#timestart = 1200; timeend = 2200
#phasemin = 0.5; phasemax = 1.5
#RVmin = -50; RVmax = 60

#sysname = '10001167'; filename = '../../RG_spectra/10001167/rvs_jean.txt'
#timestart = 1500; timeend = 2300
#phasemin = 0.49; phasemax = 1.49
#RVmin = -149; RVmax = -61

#sysname = '7037405'; filename = '../../RG_spectra/7037405_1/rvs_final.txt'
#timestart = 1550; timeend = 2550
#phasemin = 0.5; phasemax = 1.5
#RVmin = -79; RVmax = -1

#sysname = '8054233'
#timestart = 1600; timeend = 2200
#phasemin = 0.5; phasemax = 1.5
#RVmin = -35; RVmax = 20

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# usecols=(0,1,3,4,5,6) # this is the default, with RVs in 3,4,5,6 not 8,9,10,11
bjd, phase, rv1, rverr1, rv2, rverr2 = np.loadtxt(filename, comments='#', 
    dtype={'names': ('bjd', 'phase', 'rv1', 'rverr1', 'rv2', 'rverr2'),
    'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.float64)},
    usecols=(0,1, 3,4,5,6), unpack=True)

if doShift == True:
    # apply a shift to each RV before plotting it
    shiftfile = '../../TelFit/9246715_telfit/telfit_RVshifts.txt'
    try:
        shift = np.loadtxt(shiftfile, comments='#', usecols=(1,), unpack=True)
        rv1 = rv1 - shift
        rv2 = rv2 - shift
        print('Shift applied')
    except:
        print('Could not find shift file')

# Skip any RV values that have 0 for error bars
for idx, err in enumerate(rverr1):
    if err == 0:
        rv1[idx] = None
        rverr1[idx] = None
for idx, err in enumerate(rverr2):
    if err == 0:
        rv2[idx] = None
        rverr2[idx] = None
rv1mask = np.isfinite(rv1)
rv2mask = np.isfinite(rv2)

## OPTION to compare one set of RVs against another
if compareRVs == True:
    rvjean1, errjean1, rvjean2, errjean2 = np.loadtxt(filename, comments='#', usecols=(8,9,10,11), unpack=True)
    for idx, err in enumerate(errjean1):
        if err == 0:
            rvjean1[idx] = None
            errjean1[idx] = None
    for idx, err in enumerate(errjean2):
        if err == 0:
            rvjean2[idx] = None
            errjean2[idx] = None
    rvjean1mask = np.isfinite(rvjean1)
    rvjean2mask = np.isfinite(rvjean2)
    rv1 = rvjean1 - rv1 # REDEFINE RV ARRAYS AS DIFFERENCES
    rv2 = rvjean2 - rv2 # REDEFINE RV ARRAYS AS DIFFERENCES
    RVmin = -5; RVmax = 5

# Double the arrays so we can plot any phase from 0 to phase 2... assuming phase is in range (0,1)
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
plt.plot(bjd[rv1mask]-dateoffset, rv1[rv1mask], color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2mask]-dateoffset, rv2[rv2mask], color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
for idx, date in enumerate(bjd):
    plt.errorbar(date-dateoffset, rv1[idx], yerr=rverr1[idx], fmt='ko', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
    plt.errorbar(date-dateoffset, rv2[idx], yerr=rverr2[idx], fmt='ko', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset), size=24, labelpad=10)

# Draw horizontal line at RV = 0
if compareRVs == True: plt.axhline(y=0, color='k', ls=':')

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major', labelsize=20)
for idx, ph in enumerate(phase_double):
    plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
    plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
plt.xlabel("Orbital Phase", size=24)

# Draw vertical lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Draw horizontal line at RV = 0
if compareRVs == True: plt.axhline(y=0, color='k', ls=':')

# Option for a legend and labels (note: for a legend you will need to add a label to the plt.errorbar commands)
#plt.legend(ncol=2, loc=1, fontsize=20, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35), columnspacing=0.7)
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=24, rotation='vertical')
fig.text(0.14, 0.115, 'Folded', size=24)
fig.text(0.14, 0.55, 'Unfolded', size=24)
fig.text(0.2, 0.9, sysname, size=32)

plt.show()