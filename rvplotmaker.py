from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
'''
Radial velocity plotter!
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase

Columns should be as follows
BJD, PHASE, KEPLER_BJD, RV1, RV1_ERR, RV2, RV2_ERR, SOURCE
(all of these are numbers, except 'source' is text)

Update September 2015:
More generalized for plotting RVs from two instruments for one target.
Has two flag options
1. apply some shift to the RVs before plotting them
2. read in another set of calculated RVs from cols 8,9,10,11 for comparison and plot 
    RV1 = RV_col8 - RV_col3, RV2 = RV_col10 - RV_col5 instead, plus a line at RV = 0.
'''

doShift = False
compareRVs = False

# File containing times, phases, and RVs as specified below
#filename = '../../RG_spectra/8054233/rvoutfile1_arcesBF.txt'
filename = '../../RG_spectra/8702921/rvoutfile2_arcesBF.txt'
#filename = '../../RG_spectra/rvs_final_10001167.txt'
#filename = '../../KIC_8848288/rvs_BF.txt'

#sysname = '9970396'
#timestart = 1550
#timeend = 2500
#phasemin = 0.4
#phasemax = 1.4
#RVmin = -59
#RVmax = 19

#sysname = '3955867'
#timestart = 1650
#timeend = 2350
#phasemin = 0.4
#phasemax = 1.4
#RVmin = -49
#RVmax = 69

sysname = '8702921'
timestart = 1200
timeend = 2200
phasemin = 0.5
phasemax = 1.5
RVmin = -29
RVmax = 9

#sysname = '9291629'
#timestart = 1650
#timeend = 2350
#phasemin = 0.4
#phasemax = 1.4
#RVmin = -99
#RVmax = 34

# for 8848288 (not one of our RG/EBs)
#sysname = '8848288'
#timestart = 1310
#timeend = 1350
#phasemin = 0.45
#phasemax = 1.45
#RVmin = -24
#RVmax = -13

#sysname = '9246715'
#sysname = '5786154'

# for 8430105
#sysname = '8430105'
#timestart = 1200
#timeend = 2200
#phasemin = 0.5
#phasemax = 1.5
#RVmin = -50
#RVmax = 60

# for 10001167
#sysname = '10001167'
#timestart = 1500
#timeend = 2300
#phasemin = 0.48
#phasemax = 1.48
#RVmin = -149
#RVmax = -49

# for 7037405
#sysname = '7037405'
#timestart = 1550
#timeend = 2550
#phasemin = 0.4
#phasemax = 1.4
#RVmin = -79
#RVmax = -1

# for 8054233
#sysname = '8054233'
#timestart = 1600
#timeend = 2200
#phasemin = 0.5
#phasemax = 1.5
#RVmin = -35
#RVmax = 20

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# usecols=(0,1,3,4,5,6,7) # this is the default, with RVs in 3,4,5,6 not 8,9,10,11 (7=source!)
bjd, phase, rv1, rverr1, rv2, rverr2, source = np.loadtxt(filename, comments='#', 
    dtype={'names': ('bjd', 'phase', 'rv1', 'rverr1', 'rv2', 'rverr2', 'source'),
    'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, '|S15')},
    usecols=(0,1, 3,4,5,6, 7), unpack=True)

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

##################################
## OPTIONAL FOR COMPARISON PURPOSES
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
    RVmin = -5
    RVmax = 5
## OPTIONAL FOR COMPARISON PURPOSES
##################################

# Double the arrays so we can plot any phase from 0 to phase 2, if desired
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
plt.plot(bjd[rv1mask]-2454833, rv1[rv1mask], color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2mask]-2454833, rv2[rv2mask], color='0.75', mfc=None, mec=None, lw=1.5, ls=':')
for idx, label in enumerate(source):
    if label == 'arces':
        plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='ko', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
        plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='ko', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
#    elif label == 'tres':
#        plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='kv', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
#        plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='kv', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
    elif label == 'apogee':
        plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='ks', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
        plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='ks', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)                
    else:
        plt.errorbar(bjd[idx]-2454833, rv1[idx], yerr=rverr1[idx], fmt='ko', color='0.75', mfc=red, mec='k', ms=10, lw=1.5)
        plt.errorbar(bjd[idx]-2454833, rv2[idx], yerr=rverr2[idx], fmt='ko', color='0.75', mfc=yel, mec='k', ms=10, lw=1.5)
##plt.errorbar(bjd-2454833, rv1, yerr=rverr1, fmt='ko:', color='0.75', mfc=red, mec=red, ms=10, lw=1.5)
##plt.errorbar(bjd-2454833, rv2, yerr=rverr2, fmt='ko:', color='0.75', mfc=yel, mec=yel, ms=10, lw=1.5)
plt.xlabel("Time (BJD -- 2454833)", size=24, labelpad=10)

# Draw horizontal line at RV = 0
if compareRVs == True: plt.axhline(y=0, color='k', ls=':')

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
#ax1.set_xlim([0.5, 1.5])
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major', labelsize=20)
for idx, label in enumerate(source):
    idx2 = idx + len(source)
    if label == 'arces':
        plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='ARCES, Star 1' if idx==0 else '')
        plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5, label='Star 2' if idx==0 else '')
        plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
        plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
#    elif label == 'tres':
#        plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='v', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='TRES, Star 1' if idx==0 else '')
#        plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='v', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5, label='Star 2' if idx==0 else '')
#        plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='v', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
#        plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='v', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
    elif label == 'apogee':
        plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='s', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5, label='APOGEE, Star 1' if idx==16 else '')
        plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='s', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)#, label='Star 2' if idx==16 else '')        
        plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='s', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
        plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='s', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
    else:
        plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
        plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
        plt.errorbar(phase_double[idx2], rv1_double[idx2], yerr=rverr1_double[idx2], marker='o', color=red, mec='k', ecolor=red, ms=10, ls='None', lw=1.5)
        plt.errorbar(phase_double[idx2], rv2_double[idx2], yerr=rverr2_double[idx2], marker='o', color=yel, mec='k', ecolor=yel, ms=10, ls='None', lw=1.5)
plt.xlabel("Orbital Phase", size=24)

# Draw vertical lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Draw horizontal line at RV = 0
if compareRVs == True: plt.axhline(y=0, color='k', ls=':')

# Create the legend and y-axis label
plt.legend(ncol=2, loc=1, fontsize=20, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35), columnspacing=0.7)
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=24, rotation='vertical')
fig.text(0.14, 0.115, 'Folded', size=24)
fig.text(0.14, 0.55, 'Unfolded', size=24)
fig.text(0.2, 0.9, sysname, size=32)

plt.show()