import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
import sys


xMin=8.
xMax=10.	

nbBeam=4
nbPeak=nbBeam*(nbBeam-1)/2
nbTrip=nbBeam*(nbBeam-1)*(nbBeam-2)/6
maxCol=3
nbRow=(nbBeam-1) / maxCol + 1
if nbBeam < maxCol+1:
	nbCol=(nbBeam-1) % maxCol + 1
else:
	nbCol=maxCol


def getData(fileName):
	hdu = fits.open(fileName)
	data=hdu['OI_T3'].data['T3PHI']
	hdu.close()
	return(data)

def getWvl(fileName):
	hdu = fits.open(fileName)
	wvl=hdu['OI_WAVELENGTH'].data['EFF_WAVE']
	hdu.close()
	return(wvl*1e6)

dataName='/data-matisse/ComWorkshop/2017-04-27/TransFunc-01/pcr/MERGE_BCDOUT_RAW_DPHASE_CORRFLUX_CAL_MATISSE_GEN_LAMP_N117_0001.fits'
wave=getWvl(dataName)
nbWave=np.size(wave)
print(str(nbWave)+' values in input wavelength list from '+str(np.amin(wave))+' to '+str(np.amax(wave))+' mu')

index=np.where( (wave > xMin) & (wave < xMax) )
idxWvl=index[0]
nbWvl=np.size(idxWvl)
print('...among which '+str(nbWvl)+' values between '+str(xMin)+' and '+str(xMax)+' mu')
wvl=wave[idxWvl]

nbWin=nbTrip
idxWin=np.arange(nbWin)
data=np.ndarray((nbWin,nbWvl))
inDat = getData(dataName)
dat_=inDat[idxWin,:]
data=dat_[:,idxWvl]

xLab=r'$\lambda$ ($\mu$m)'
yMin=np.amin(data)
yMax=np.amax(data)

yLab=r'phase closure ($^\circ$)'

subTitle=r'baseline\#'

supTitle=dataName.replace('_','\_')

med=np.median(data,axis=1)
print(r'medians of clos data:')
print(med)

iqr=np.subtract(*np.percentile(data, [75, 25],axis=1))
print(r'interquartile ranges:')
print(iqr)

plt.figure(num=0,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	histo=data[numWin,:]
	plt.hist(histo, 10, normed=0, alpha=0.3,label=str(1))	
	plt.xlabel(yLab)
	plt.ylabel('Density')
	plt.title("%s%d: med=%7.1f$^\circ$, iqr=%6.1f$^\circ$" % (subTitle, numWin+1, med[numWin], iqr[numWin]))
	plt.grid()
plt.suptitle(supTitle)


plt.show()



