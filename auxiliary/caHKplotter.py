import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import glob
from echelle_funcs import getOrderWavelengths
'''
Plots the CaII H&K line region from a series of 2D (multi-order) ARCES spectra. Hooray!

Note: This program assumes a spectrum of 'KIC 1234567' will have the string '1234567' 
somewhere in its filename. If not, you're going to have a bad time.
'''
# Define where important things live, like your spectra and a text file listing your targets
rootdir = '../../../RG_spectra/tequila_201604/'
targets = np.loadtxt('starnames.txt', usecols=(0,), comments='#', unpack=True)
strtargets = []
for target in targets:
    strtargets.append(str(int(target)))

#########
## OPTION 1: create a list of spectrum files via glob
#dirlist = []
#for target in strtargets:
#    dirlist.append(rootdir + 'kplr' + target + '/')
#filelist = []
#for dir in dirlist:
#    for file in glob.glob(dir + '1*ec.fits'):
#        filelist.append(file)
#        #print(file)
#########

#########
## OPTION 2: Once a text file listing the names of the good 2D spectra exists, use that
infiles = '../../../Rawls_tidesactivity/good_caHK_spectra.txt'
filelist = []
with open(infiles) as f1:
    for line in f1:
        filelist.append(rootdir + line.rstrip())
useRVdata = True # set True to overplot where each star's Ca line contribution should be
#########

carest = [3968.468, 3933.663]
if useRVdata == True:
    # (1) Get timestamps for each observation from jean's [target].info.txt files
    rtargets = []; reffiles = []; bjdsaves = []; bcvsaves = []
    for file in filelist:
        ### the below is SUPER GREAT but actually not relevant here ###
        #fitsspec = fits.open(file)
        #fitsspec.verify('silentfix')
        #head = fitsspec[0].header
        #utmiddle = head['UTMIDDLE']
        #if '=' in utmiddle:
        #    time = Time(utmiddle[1::], format='fits', scale='utc')
        #else:
        #    time = Time(head['DATE-OBS'][0:11] + utmiddle, format='fits', scale='utc')
        #target = head['OBJNAME']
        #print(target, time.jd)
        #fitsspec.close()
        ### the above is SUPER GREAT but actually not relevant here ###
        fitsspec = fits.open(file)
        target = fitsspec[0].header['OBJNAME']
        fitsspec.close()
        infofile = rootdir + 'kplr' + target + '/' + target + '.info.txt'
        fullspecfiles = np.loadtxt(infofile, usecols=(0,), comments='#', dtype={'names':('fullspecfiles',),'formats':('|S29',)}, unpack=True)
        bjds, bcvs = np.loadtxt(infofile, usecols=(3,4), comments='#', unpack=True)
        for fullspecfile, bjd, bcv in zip(fullspecfiles, bjds, bcvs):
            if file[-19::] in str(fullspecfile):
                rtargets.append(target)
                reffiles.append(file)
                bjdsaves.append(bjd)
                bcvsaves.append(bcv)  
    if len(filelist) != len(reffiles):
        raise IOError('Something bad happened when you tried to assign bjds and bcvs with the .info.txt files.')
    # (2) use bjds to look up both stars' RVs in patrick's LaTeX file of doom
    texin = '../../../Rawls_tidesactivity/rvdata_gaulme2016.txt'
    obsdates, rv1s, rv2s = np.loadtxt(texin, comments='%', usecols=(0,1,3), unpack=True)
    targetlist = []; kepjdlist = []; bcvlist = []; rv1list = []; rv2list = []; newfilelist = []
    for file, target, bjd, bcv in zip(filelist, rtargets, bjdsaves, bcvsaves):
        kepjd = bjd - 2454833.
        for obsdate, rv1, rv2 in zip(obsdates, rv1s, rv2s):
            timecompare = kepjd - obsdate
            if np.abs(timecompare) < 0.0001:
                diff = kepjd - obsdate
                #print(target, diff, rv1, rv2, bcv) # this looks good!
                # save parallel lists
                targetlist.append(target)
                kepjdlist.append(kepjd)
                bcvlist.append(bcv)
                rv1list.append(rv1)
                rv2list.append(rv2)
                newfilelist.append(file)
    # (3) properly combine RVs and BCVs to get predicted location of CaII H and CaII K
    castar1list = []; castar2list = []
    for bcv, rv1, rv2 in zip(bcvlist, rv1list, rv2list):
        drv1 = rv1 - bcv
        drv2 = rv2 - bcv
        castar1 = [carest[0] + carest[0]*drv1/299792.458, carest[1] + carest[1]*drv1/299792.458]
        castar2 = [carest[0] + carest[0]*drv2/299792.458, carest[1] + carest[1]*drv2/299792.458]
        castar1list.append(castar1)
        castar2list.append(castar2)
    # (4) check and see if filelist from infiles and newfilelist from texin agree
    print(len(filelist), len(newfilelist))
    print('These files were in the original filelist but now are not: ', 
            [x for x in filelist if x not in newfilelist])
    filelist = newfilelist


def plotThreeAtATime(filelist):
    '''
    The original plotter is now a function!
    Use this if you're not interested to sort by star yet and just want to plot all
    the spectra in figures, three at a time, for manual inspection.
    It builds up a list of spectra with S/N > 5, defined as having a y-axis values > 25.
    '''
    ## Plot three panes of spectra at a time
    ## (The original idea was to plot a range of orders, but then I realized I only wanted two)
    #orderstart = 90 #89
    #orderstop = 93 #94
    #colors = ['#fc8d59','#ef6548','#d7301f','#b30000','#7f0000']
    ## LOL NOPE let's just manually plot order 90 and 91, which contain the Ca H&K lines
    goodfilelist = []
    for i in np.arange(0, len(filelist)-3, 3):
        fig = plt.figure(figsize=(10,18))
        # LET'S PLOT ONE SPECTRUM in the top of three panels
        ax = fig.add_subplot(311)
        ax.set_xlim([3910, 3990])
        fitsspec = fits.open(filelist[i])
        fitsfluxes = fitsspec[0].data
        fitswaves = getOrderWavelengths(fitsspec)
        # sort the damn orders so they're not backwards
        for idx, (waves, fluxes) in enumerate(zip(fitswaves, fitsfluxes)):
            fitsfluxes[idx] = fluxes[np.argsort(waves)]
            fitswaves[idx] = waves[np.argsort(waves)]
        # truncate the spectra somewhat arbitrarily (one order for H line region, one order for K line)
        idxtrunc90 = np.where(fitswaves[90] > 3950)[0][0]
        idxtrunc91 = np.where(fitswaves[91] > 3950)[0][0]
        waves90 = fitswaves[90][idxtrunc90::]; fluxes90 = fitsfluxes[90][idxtrunc90::]
        waves91 = fitswaves[91][0:idxtrunc91]; fluxes91 = fitsfluxes[91][0:idxtrunc91]
        if (np.max(fluxes90) < 25) or (np.max(fluxes91) < 25) or (len(fluxes90) < 50):
            colors = ['0.75', '0.75'] # this is a bad spectrum so we plot it in gray and feel bad for it
        else:
            colors = ['#3182bd','#08519c'] # this is a good spectrum so it gets pretty blue colors
            goodfilelist.append(filelist[i])
        plt.plot(waves90, fluxes90, color=colors[0])
        plt.plot(waves91, fluxes91, color=colors[1])
        ## (You could use this with the original idea of plotting a range of orders!)
    #    for idx, (waves, fluxes) in enumerate(zip(fitswaves[orderstart:orderstop], fitsfluxes[orderstart:orderstop])):
    #        plt.plot(waves, fluxes, color=colors[idx])
        plt.axvline(x=3968.468, ls=':', color='k')
        plt.axvline(x=3933.663, ls=':', color='k')
        plt.title(filelist[i][-28:], size=20)
        fitsspec.close()
        # DO IT ALL AGAIN for the middle panel
        ax = fig.add_subplot(312)
        ax.set_xlim([3910, 3990])
        fitsspec = fits.open(filelist[i+1])
        fitsfluxes = fitsspec[0].data
        fitswaves = getOrderWavelengths(fitsspec)
        for idx, (waves, fluxes) in enumerate(zip(fitswaves, fitsfluxes)):
            fitsfluxes[idx] = fluxes[np.argsort(waves)]
            fitswaves[idx] = waves[np.argsort(waves)]
        idxtrunc90 = np.where(fitswaves[90] > 3950)[0][0]
        idxtrunc91 = np.where(fitswaves[91] > 3950)[0][0]
        waves90 = fitswaves[90][idxtrunc90::]; fluxes90 = fitsfluxes[90][idxtrunc90::]
        waves91 = fitswaves[91][0:idxtrunc91]; fluxes91 = fitsfluxes[91][0:idxtrunc91]
        if (np.max(fluxes90) < 25) or (np.max(fluxes91) < 25) or (len(fluxes90) < 50):
            colors = ['0.75', '0.75']
        else:
            colors = ['#3182bd','#08519c']
            goodfilelist.append(filelist[i+1])
        plt.plot(waves90, fluxes90, color=colors[0])
        plt.plot(waves91, fluxes91, color=colors[1])
        plt.axvline(x=3968.468, ls=':', color='k')
        plt.axvline(x=3933.663, ls=':', color='k')
        plt.title(filelist[i+1][-28:], size=20)
        fitsspec.close()
        # DO IT ALL ONE LAST TIME for the bottom panel
        ax = fig.add_subplot(313)
        ax.set_xlim([3910, 3990])
        fitsspec = fits.open(filelist[i+2])
        fitsfluxes = fitsspec[0].data
        fitswaves = getOrderWavelengths(fitsspec)
        for idx, (waves, fluxes) in enumerate(zip(fitswaves, fitsfluxes)):
            fitsfluxes[idx] = fluxes[np.argsort(waves)]
            fitswaves[idx] = waves[np.argsort(waves)]
        idxtrunc90 = np.where(fitswaves[90] > 3950)[0][0]
        idxtrunc91 = np.where(fitswaves[91] > 3950)[0][0]
        waves90 = fitswaves[90][idxtrunc90::]; fluxes90 = fitsfluxes[90][idxtrunc90::]
        waves91 = fitswaves[91][0:idxtrunc91]; fluxes91 = fitsfluxes[91][0:idxtrunc91]
        if (np.max(fluxes90) < 25) or (np.max(fluxes91) < 25) or (len(fluxes90) < 50):
            colors = ['0.75', '0.75']
        else:
            colors = ['#3182bd','#08519c']
            goodfilelist.append(filelist[i+2])
        plt.plot(waves90, fluxes90, color=colors[0])
        plt.plot(waves91, fluxes91, color=colors[1])
        plt.axvline(x=3968.468, ls=':', color='k')
        plt.axvline(x=3933.663, ls=':', color='k')
        plt.title(filelist[i+2][-28:], size=20)
        plt.xlabel('Wavelength (\AA)')
        fitsspec.close()
        plt.show()
    for file in goodfilelist:
        print(file)

## Use this if you want to plot spectra three at a time via the function above
#plotThreeAtATime(filelist)

## Use this if you want to make one figure per star with all its spectra stacked and normalized
#fig = plt.figure() # UNCOMMENT FOR MULTIPANE PLOT
for axdx, star in enumerate(strtargets):
    yoffset = 0
    fig = plt.figure(figsize=(9,12)); ax = fig.add_subplot(111) # UNCOMMENT FOR ONE PLOT PER STAR
    #ax = fig.add_subplot(3, 5, axdx+1) # UNCOMMENT FOR MULTIPANE PLOT
    plt.title(star)
    for fidx, file in enumerate(filelist):
        if star in file:
            #print(star, file) # I can write nested loops yes I can
            fitsspec = fits.open(file)
            fitsfluxes = fitsspec[0].data
            fitswaves = getOrderWavelengths(fitsspec)
            for idx, (waves, fluxes) in enumerate(zip(fitswaves, fitsfluxes)):
                fitsfluxes[idx] = fluxes[np.argsort(waves)]
                fitswaves[idx] = waves[np.argsort(waves)]
            idxtrunc90 = np.where(fitswaves[90] > 3950)[0][0]
            idxtrunc91 = np.where(fitswaves[91] > 3950)[0][0]
            waves90 = fitswaves[90][idxtrunc90::]; fluxes90 = fitsfluxes[90][idxtrunc90::]
            waves91 = fitswaves[91][0:idxtrunc91]; fluxes91 = fitsfluxes[91][0:idxtrunc91]
            waves = np.concatenate((waves91, waves90), axis=0)
            fluxes = np.concatenate((fluxes91, fluxes90), axis=0)
            # NORMALIZATION! look out kittens, these fluxes all go from 0 to 1 now
            fluxes = (fluxes-np.min(fluxes))/(np.max(fluxes)-np.min(fluxes))
            plt.plot(waves, fluxes+yoffset, color='#3182bd')
            #print(file[-19:], star, fidx) # loop index consistency check!!
            if useRVdata == True:
                plt.plot(castar1list[fidx][0], yoffset+0.5, ls='None', marker='|', ms=25, mew=1.5, color='r')
                plt.plot(castar1list[fidx][1], yoffset+0.5, ls='None', marker='|', ms=25, mew=1.5, color='r')
                plt.plot(castar2list[fidx][0], yoffset+0.5, ls='None', marker='|', ms=25, mew=1.5, color='k')
                plt.plot(castar2list[fidx][1], yoffset+0.5, ls='None', marker='|', ms=25, mew=1.5, color='k')
            yoffset += 1
            fitsspec.close()
    if useRVdata == False:
        plt.axvline(x=carest[0], ls=':', color='k')
        plt.axvline(x=carest[1], ls=':', color='k')
    plt.axis([3925,3975,-0.1,yoffset+0.1])
    ax.set_yticklabels([])
    ax.set_yticks([])
    plt.xlabel('Wavelength (\AA)')
    plt.show() # UNCOMMENT FOR ONE PLOT PER STAR
#plt.show() # UNCOMMENT FOR MULTIPANE PLOT

