from __future__ import print_function
import numpy as np
from astropy.time import Time
'''
Read in Kepler BJD, RV, RVerror for boths tars
Also get orbital period and zeropoint to calculate phases

Return formatted LaTeX table entries as follows:
UTC Date, Midpoint in BJD-2450000, Phase, RV1(err), RV2(err)
'''
dir = '../../RG_ELCmodeling/'

#KIC = '3955867'
#trialresultfile = 'trials5_6_results.txt'

#KIC = '10001167'
#trialresultfile = 'trial1_results.txt'

KIC = '7037405'
trialresultfile = 'trial7_finalresults.txt'

MSinfile = dir + KIC + '/rvs_' + KIC +'_MS.txt'
RGinfile = dir + KIC + '/rvs_' + KIC +'_RG.txt'
ELCresults = dir + KIC + '/' + trialresultfile

def phasecalc(times, period=100, BJD0=2454833):
	'''
	Function to calculate orbitsl phase for a list of times given period and zeropoint
	'''
	phases = []
	#cycles = []
	for i in range(0, len(times)):
		fracP = (times[i] - BJD0) / period
		if fracP < 0:
			phases.append(fracP % 1)
			#cycles.append(int(fracP))
		else:
			phases.append(fracP % 1)
			#cycles.append(int(fracP) + 1)
		#print(fracP, phases[i])
	return np.array(phases)

def numToTxtMonth(numstring='01'):
    '''
    Function to turn a 2-'digit' string into a three-character month
    e.g., 01 = Jan, 02 = Feb, etc.
    '''
    if numstring == '01': return('Jan')
    elif numstring == '02': return('Feb')
    elif numstring == '03': return('Mar')
    elif numstring == '04': return('Apr')
    elif numstring == '05': return('May')
    elif numstring == '06': return('Jun')
    elif numstring == '07': return('Jul')
    elif numstring == '08': return('Aug')
    elif numstring == '09': return('Sep')
    elif numstring == '10': return('Oct')
    elif numstring =='11': return('Nov')
    elif numstring == '12': return('Dec')
    else: return('???')

# Get period and t0 from ELC results file
ELCoutdata = open(ELCresults)
for row in ELCoutdata:
    row = row.split()
    if row[0] == 'period':
        period = float(row[2])
    if row[0] == 't0v2':
        t0 = float(row[2])

print('KIC ' + KIC +':', period, t0) # this zeropoint is in units of Kepler BJD!

MStimes, MSrvs, MSrverrs = np.loadtxt(MSinfile, unpack=True)
RGtimes, RGrvs, RGrverrs = np.loadtxt(RGinfile, unpack=True)

# calculate orbital phases in Kepler BJD units
MSphases = phasecalc(MStimes, period, t0)
RGphases = phasecalc(RGtimes, period, t0)

# put into real BJD after calculating orbital phases
MStimes = Time(MStimes + 2454833, format='jd', scale='utc')
RGtimes = Time(RGtimes + 2454833, format='jd', scale='utc')

MStimestrings = []
RGtimestrings = []
for time in MStimes:
    month = numToTxtMonth(str(time.iso[5:7]))
    MStimestrings.append(str(time.iso[0:4]) + ' ' + month + ' ' + str(time.iso[8:10]))
for time in RGtimes:
    month = numToTxtMonth(str(time.iso[5:7]))
    RGtimestrings.append(str(time.iso[0:4]) + ' ' + month + ' ' + str(time.iso[8:10]))


for time, BJD, phase, rvRG, rverrRG in zip(RGtimestrings, RGtimes, RGphases, RGrvs, RGrverrs):
    errstringRG = format(rverrRG, '.2f')[2:]
    if errstringRG[0] == '0': errstringRG = errstringRG[1:]
    if BJD in MStimes:
        idx = np.where(MStimes == BJD)
        rvMS = MSrvs[idx[0][0]]
        rverrMS = MSrverrs[idx[0][0]]
        errstringMS = format(rverrMS, '.2f')[2:]
    else:
        rvMS = None
        rverrMS = None
        errstringMS = None
    BJD = BJD.jd - 2450000
    if rvMS != None:
        print('{0} & {1:.6f} & {2:.3f} & {3:.2f}({4}) & {5:.2f}({6}) \\\\'.format(time, BJD, phase, rvRG, errstringRG, rvMS, errstringMS))
    else:
        nodata = '\\nodata'
        print('{0} & {1:.6f} & {2:.3f} & {3:.2f}({4}) & {5} \\\\'.format(time, BJD, phase, rvRG, errstringRG, nodata))
    
    
