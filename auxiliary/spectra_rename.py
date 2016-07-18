from __future__ import print_function
import glob
from shutil import copy2
'''
Specialized program to copy files sorted by date into directories sorted by target.

It's really more of a spectrum MOVER than a spectrum renamer...

My situation: a set of observations in directories named yymmdd have images corresponding
to various stars. All the filenames are 'yymmdd.xxxx.ec.fits.' Star directories exist,
and have some files with the corresponding 'yymmdd.xxxx' string in them.
This program: identifies which files should go where and COPY-PASTA!
'''

#stardir = '../../../../../Volumes/Rocinante/3.5m_spectra/combined_bystar'
stardir = '../../RG_spectra/tequila_201604'
timedir = '../../../../../Volumes/Rocinante/3.5m_spectra'

#timedirs = []
#for time in glob.glob(timedir+'/1*'):
#    timedirs.append(time)

stardirs = []
for star in glob.glob(stardir+'/kplr*'):
    stardirs.append(star)

for dir in stardirs:
    star = dir[-11:]
    fullspecs = glob.glob(dir+'/fullspec*ec.fits')
    for spec in fullspecs:
        datestring = spec[-19:-8]
        #print(star, datestring)
        filetocopy = timedir + '/' + datestring[0:6] + '/' + datestring + '.ec.fits'
        #print('Copying {0} to {1}'.format(filetocopy, dir+'/.'))
        try:
            copy2(filetocopy, dir+'/.')
        except:
            print('ERROR could not copy {0} to {1}'.format(filetocopy, dir+'/.'))
print('Files done copying.')