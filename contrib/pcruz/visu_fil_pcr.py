import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
import sys


fname=sys.argv[-1]
f=open(fname,"r")
lines=f.read().splitlines()
f.close()

dataType=lines[0]
wBand=lines[1]

if wBand == 'L':
	xMin=2.8
	xMax=4.2
	xMin=3.5
	xMax=4.
elif wBand == 'M':
	xMin=4.5
	xMax=5
elif wBand == 'LM':
	xMin=2.8
	xMax1=4.2
	xMin1=4.5
	xMax=5
elif wBand == 'L-':
	xMin=2.8
	xMax=3.5
elif wBand == 'N':
	xMin=8
	xMax=10
else:
	print('second argument must be L or M or LM or N or L-')
	sys.exit(0)
	

nbBeam=4
nbPeak=nbBeam*(nbBeam-1)/2
nbTrip=nbBeam*(nbBeam-1)*(nbBeam-2)/6
maxCol=3
nbRow=(nbBeam-1) / maxCol + 1
if nbBeam < maxCol+1:
	nbCol=(nbBeam-1) % maxCol + 1
else:
	nbCol=maxCol
if dataType=='flux':
	nbWin=nbBeam
elif dataType=='clos':
	nbWin=nbTrip
else:
	nbWin=nbPeak

def getData(fileName):
	hdu = fits.open(fileName)
	if dataType=='vis2'or dataType=='cflux' or dataType=='vis':
		vis2=hdu['OI_VIS2'].data['VIS2DATA']
		data=vis2
		if dataType=='vis' or dataType=='cflux':
			data=np.sqrt(np.absolute(vis2))*np.sign(vis2)
	elif dataType=='flux':
		data=hdu['OI_FLUX'].data['FLUXDATA']
	elif dataType=='dphas':
		data=hdu['OI_VIS'].data['VISPHI']
	elif dataType=='clos':
		data=hdu['OI_T3'].data['T3PHI']
	else:
		print('first argument must be vis2 or vis or flux or clos or dphas or cflux')
		sys.exit(0)
	hdu.close()
	return(data)

def getWvl(fileName):
	hdu = fits.open(fileName)
	wvl=hdu['OI_WAVELENGTH'].data['EFF_WAVE']
	hdu.close()
	return(wvl*1e6)

def getCorflag(fileName):
	hdu = fits.open(fileName)
	flag=hdu[0].header['HIERARCH ESO PRO REC1 PARAM9 VALUE']
	hdu.close()
	return(flag)	

def getFlux(fileName):
	hdu = fits.open(fileName)
	data=hdu['IMAGING_DATA'].data['DATA3']
	tartyp=hdu['IMAGING_DATA'].data['TARTYP']
	hdu.close()
	indexTarget=np.where(tartyp=='T')
	indexSky=np.where(tartyp=='S')
	fluxTarget=np.mean(data[indexTarget])
	fluxSky=np.mean(data[indexSky])
	return(fluxTarget-fluxSky)	

nbFile=int(lines[2])
print(str(nbFile)+' input data files')
if nbFile < 2:
	print('number of files must be greater than 1')
	sys.exit(0)

dataName1=lines[3]
corflag=getCorflag(dataName1)
if corflag=='false' and dataType=='cflux' :
	print('no correlated flux calculated => use another argument') 
	sys.exit(0)	
if corflag=='true' and (dataType=='vis' or dataType=='vis2' or dataType=='flux') :
	print('vis2 table contains quadratic correlated flux and no flux is calculated => use another argument') 
	sys.exit(0)	

wave=getWvl(dataName1)
nbWave=np.size(wave)
print(str(nbWave)+' values in input wavelength list from '+str(np.amin(wave))+' to '+str(np.amax(wave))+' mu')
#print(wave)

if wBand == 'LM':
	index=np.where( ((wave > xMin) & (wave < xMax1)) | ((wave > xMin1) & (wave < xMax)) )
else:
	index=np.where( (wave > xMin) & (wave < xMax) )
idxWvl=index[0]
nbWvl=np.size(idxWvl)
print(str(nbWvl)+' values in '+wBand+' band between '+str(xMin)+' and '+str(xMax)+' mu')
wvl=wave[idxWvl]
if nbWvl < 1:
	print('change band or input files')
	sys.exit(0)

idxWin=np.arange(nbWin)
data=np.ndarray((nbFile,nbWin,nbWvl))
dat=np.ndarray((nbFile,nbWin,nbWave))
for numFile in range(nbFile):
	inDat = getData(lines[3+numFile])
	dat_=inDat[idxWin,:]
	data_=dat_[:,idxWvl]
	for numWin in range(nbWin):
		dat[numFile,numWin,:]=dat_[numWin,:]
		data[numFile,numWin,:]=data_[numWin,:]

if dataType=='cflux':
	flux=np.ndarray((nbFile))
	ratio=np.ndarray((nbFile))
	for numFile in range(nbFile):
#		print(lines[3+numFile])
		flux[numFile] = getFlux(lines[3+nbFile+numFile])
		ratio[numFile]=	flux[numFile]/flux[0]		
		for numWin in range(nbWin):
			data[numFile,numWin,:]=data[numFile,numWin,:]/ratio[numFile]
	print(r'mean transmissions vs fil/pol=open (first file)')
	print(ratio)

meanData=np.mean(data,axis=0)
maxData=np.amax(data,axis=0)
minData=np.amin(data,axis=0)

relDiff=np.absolute((maxData-minData)/meanData)

if np.shape(data[0,:,:]) != np.shape(relDiff) :
	print('data[0,:,:]',np.shape(data[0,:,:]))
	print('relDiff',np.shape(relDiff))

nbDif=nbFile-1
rdifData=np.ndarray((nbDif,nbWin,nbWvl))
for numDif in range(nbDif):
	for numWin in range(nbWin):
		rdifData[numDif,numWin,:]=(data[numDif+1,numWin,:]-data[0,numWin,:])/data[0,numWin,:]


ddata=np.ndarray((nbFile,nbWin,nbWvl))
rddata=np.ndarray((nbFile,nbWin,nbWvl))
for numFile in range(nbFile):
	for numWin in range(nbWin):
		ddata[numFile,numWin,:]=data[numFile,numWin,:]-meanData[numWin,:]
		rddata[numFile,numWin,:]=ddata[numFile,numWin,:]/meanData[numWin,:]
maxDdata=np.amax(ddata,axis=0)
minDdata=np.amin(ddata,axis=0)
maxRddata=np.amax(rddata,axis=0)
minRddata=np.amin(rddata,axis=0)

mean=np.mean(data,axis=2)
var=np.sum((data-dat[:,:,idxWvl+1])**2,axis=2)/(nbWvl-1)/2.
sigma=np.sqrt(var)
print(r'sigma:')
print(sigma)
snr=mean/sigma
print(r'SNR:')
print(snr)
covar=np.mean(data*dat[:,:,idxWvl+1],axis=2)
print(r'covar/var:')
print(covar/var)

xLab=r'$\lambda$ ($\mu$m)'
outName='visu_'+wBand+'_'+dataType
yMin=np.amin(data)
yMax=np.amax(data)

yLab1=r''+dataType+''
yLab2=r''+yLab1+'-$<$'+yLab1+'$>_{fil}$'
yLab3=r'('+yLab1+'-$<$'+yLab1+'$>_{fil}$)/$<$'+yLab1+'$>_{fil}$'
yLab4=r'('+yLab1+'-$<$'+yLab1+'$>_{fil}$)/$\sigma_{'+yLab1+'}$'
yLab5=r'('+yLab1+'-'+yLab1+'$_{1}$)/'+yLab1+'$_{1}$'

subTitle=r'Base'
if dataType=='vis' :
	numLoc=4
elif dataType=='vis2':
	numLoc=4
elif dataType=='cflux' :
	yLab1 += r' (ADU)'
	yLab2 += r' (ADU)'
	numLoc=1
elif dataType=='flux':
	yLab1 += r' (ADU)'
	yLab2 += r' (ADU)'
	numLoc=1
	subTitle=r'Beam'
elif dataType=='dphas':
	yLab1 += r' ($^\circ$)'
	yLab2 += r' ($^\circ$)'
	numLoc=1
else:
	yLab1 += r' ($^\circ$)'
	yLab2 += r' ($^\circ$)'
	numLoc=4
	subTitle=r'Closure'

supTitle_=''
for numFile in range(nbFile):
 	supTitle_ += 'FILE '+str(numFile+1)+': '+ lines[3+numFile]
	if numFile != nbFile-1:
		supTitle_ += r'\newline '
supTitle=supTitle_.replace('_','\_')

plt.figure(num=0,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	for numFile in range(nbFile):
		plt.plot(wvl,data[numFile,numWin,:],'o',label=numFile+1, linestyle='solid', marker="None")
	plt.xlim([xMin,xMax])
#	plt.ylim([yMin,yMax])
	plt.ylim([np.amin(minData[numWin,:]),np.amax(maxData[numWin,:])])
	plt.xlabel(xLab)
	plt.ylabel(yLab1)
	plt.title(subTitle+str(numWin+1))		
	plt.legend(loc=numLoc)
	plt.grid()	
plt.suptitle(supTitle)
plt.savefig(outName+'.jpg')
print('saving plots in file '+outName+'.jpg')

plt.figure(num=1,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	for numFile in range(nbFile):
		plt.plot(wvl,rddata[numFile,numWin,:],'o',label=numFile+1, linestyle='solid', marker="None")
	plt.xlim([xMin,xMax])
	plt.ylim([np.amin(minRddata[numWin,:]),np.amax(maxRddata[numWin,:])])
	plt.xlabel(xLab)
	plt.ylabel(yLab3)
	plt.title(subTitle+str(numWin+1))		
	plt.legend(loc=numLoc)
	plt.grid()	
plt.suptitle(supTitle)
plt.savefig(outName+'_dif.jpg')
print('saving plots in file '+outName+'_dif.jpg')

plt.figure(num=2,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	for numFile in range(nbFile):
		histo=rddata[numFile,numWin,:]
		plt.hist(histo, 10, normed=0, alpha=0.3,label=str(numFile+1))	
	plt.xlabel(yLab3)
	plt.ylabel('Density')
	plt.title(subTitle+str(numWin+1))		
	plt.legend(loc=1)
	plt.grid()
plt.suptitle(supTitle)
plt.savefig(outName+'_histodif.jpg')
print('saving histograms in file '+outName+'_histodif.jpg')

med=np.median(rddata,axis=2)
print(r'median:')
print(med)

med=np.absolute(med)

print(r'mean and stdev of absolute medians:')
globalMean=np.mean(med)
globalStd=np.std(med,ddof=1)
print(globalMean,globalStd)

iqr=np.subtract(*np.percentile(rddata, [75, 25],axis=2))
print(r'interquartile range:')
print(iqr)

plt.figure(num=3,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	for numFile in range(nbFile):
		histo=ddata[numFile,numWin,:]/sigma[numFile,numWin]
		plt.hist(histo, 10, normed=0, alpha=0.3,label=str(numFile+1))	
	plt.xlabel(yLab4)
	plt.ylabel('Density')
	plt.title(subTitle+str(numWin+1))		
	plt.legend(loc=1)
	plt.grid()
plt.suptitle(supTitle)
plt.savefig(outName+'_histosig.jpg')
print('saving histograms in file '+outName+'_histosig.jpg')

plt.figure(num=4,figsize=(8*nbCol,6*nbRow))
for numWin in range(nbWin):
	numString=str(nbRow)+str(nbCol)+str(numWin+1)
	numVal=int(numString)
	plt.subplot(numVal)
	for numDif in range(nbDif):
		lab=str(numDif+2)+'-1'
		plt.plot(wvl,rdifData[numDif,numWin,:],'o', label=lab, linestyle='solid', marker="None")
		plt.ylim([np.amin(rdifData[numDif,numWin,:]),np.amax(rdifData[numDif,numWin,:])])
	plt.xlim([xMin,xMax])
	plt.xlabel(xLab)
	plt.ylabel(yLab5)
	plt.title(subTitle+str(numWin+1))		
	plt.legend(loc=numLoc)
	plt.grid()	
plt.suptitle(supTitle)
plt.savefig(outName+'_dif2.jpg')
print('saving plots in file '+outName+'_dif2.jpg')

#plt.show()



