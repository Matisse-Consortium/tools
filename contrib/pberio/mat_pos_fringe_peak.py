# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 09:09:57 2017

@author: pbe
"""
def gaus(x,a,b,x0,sigma):
    return b+a*np.exp(-(x-x0)**2/(2*sigma**2))

def wavelength(px,corner,c):
    return c[0]+c[1]*(px+corner)+c[2]*(px+corner)**2

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

c=np.zeros(3,dtype=np.float)
hdu = fits.open('/data/users/pbe/RunNovember2019/Iter1/TEST/Mfil/OBJ_CORR_FLUX_0001.fits')
ftreal=hdu['OBJ_CORR_FLUX'].data['CORRFLUXREAL1']
ftimag=hdu['OBJ_CORR_FLUX'].data['CORRFLUXIMAG1']
corner=hdu['IMAGING_DETECTOR'].data['CORNER'][0][1]
naxis=hdu['IMAGING_DETECTOR'].data['NAXIS'][0]
c[0]=hdu[0].header['HIERARCH PRO DISP COEF0']
c[1]=hdu[0].header['HIERARCH PRO DISP COEF1']
c[2]=hdu[0].header['HIERARCH PRO DISP COEF2']
hdu.close()

dsp=(ftreal**2+ftimag**2).mean(axis=0)


posFringePeak=np.zeros((6,naxis[1]),dtype=np.float)

wave=np.zeros(naxis[1],dtype=np.float)

for iWave in range(naxis[1]):
    wave[iWave]=wavelength(iWave,corner,c)
    if (wave[iWave] > 6.):
        dOverLambda=(naxis[0]/72.)*(8.15/wave[iWave])
    else:
        dOverLambda=(naxis[0]/72.)*(3.25/wave[iWave])
    for iBase in range(6):
        uPixelMin=np.around((3*(iBase+1)-1)*dOverLambda,0)+(naxis[0]/2)
        uPixelMax=np.around((3*(iBase+1)+1)*dOverLambda,0)+(naxis[0]/2)
        if uPixelMax >=naxis[0]:
            uPixelMax=naxis[0]-1
            uPixelMin=naxis[0]-2
        if uPixelMax==(naxis[0]/2):
            uPixelMax+=1
        uPixelMin=np.int(uPixelMin)
        uPixelMax=np.int(uPixelMax)
            
        try:
            n=uPixelMax-uPixelMin
            x=np.arange(n)+uPixelMin
            y=dsp[iWave,uPixelMin:uPixelMax]
            
            popt,pcov = curve_fit(gaus,x,y,p0=[np.max(y),0.,(uPixelMax+uPixelMin)/2. ,3.])
            posFringePeak[iBase,iWave]=popt[2]
        except RuntimeError:
            posFringePeak[iBase,iWave]=np.argmax(dsp[iWave,uPixelMin:uPixelMax])+uPixelMin
        #print iBase,popt[2],np.argmax(dsp[iWave,uPixelMin:uPixelMax])+uPixelMin

deg=3
coefPol=np.zeros((6,deg+1),dtype=np.float)
vecw=np.zeros(naxis[1],dtype=np.float)
for iWave in range(naxis[1]):
#    if ((wave[iWave] > 3.9 and wave[iWave] <4.1)  ):
    if ((wave[iWave] > 4.6 and wave[iWave] <5.1)  ):
#    if ((wave[iWave] > 2.8 and wave[iWave] <4) or (wave[iWave] > 4.4 and wave[iWave] < 5.2) ):
#    if ( (wave[iWave] > 8. and wave[iWave] < 13.) ):
        vecw[iWave]=1.
for iBase in range(6):
    #plt.plot(wave[:],posFringePeak[iBase,:])
    coefPol[iBase,:]=np.polyfit(wave,posFringePeak[iBase,:],deg,w=vecw)
    print coefPol[iBase,:]
#plt.show()
plt.imshow(np.log(dsp))
x=np.arange(naxis[1])
for iBase in range(6):
    p = np.poly1d(coefPol[iBase,:])
    plt.scatter(p(wave),x,marker='.',s=0.1,color='red')
plt.show()
    
for iBase in range(6):
    plt.figure(iBase)
    p = np.poly1d(coefPol[iBase,:])
    plt.plot(wave,posFringePeak[iBase,:],'o')
    plt.plot(wave,p(wave))
    plt.xlim([4.6,5.1])
    #plt.xlim([3.9,4.1])
plt.show()

#lam=np.zeros(naxis[1],dtype=np.float)
#x=np.arange(9)+7
##x=np.arange(5)+2
#for iWave in range(naxis[1]):
#    p = np.poly1d(coefPol[5,:])
#    lam[iWave]=wave[123]*(p(wave[123])-(naxis[0]/2))/(p(wave[iWave])-(naxis[0]/2))
#plt.figure(10)
#plt.plot(lam,wave)
#plt.plot(x,x)
##plt.xlim([7,14])
#plt.ylim([6.5,14.5])


#flux=np.zeros(naxis[1],dtype=np.float)
#f = open('spectrum_plasticfoil_N_LOW.dat', 'r')
#cpt=0
#for line in f:
#    line = line.strip()
#    columns = line.split()
#    flux[cpt] = float(columns[1])
#    cpt=cpt+1
#f.close()
#
#plt.figure(11)
#plt.plot(flux)
#plt.figure(12)
#plt.plot(lam)



# Calcul loi N LOW a partir de N HIGH
#ix=np.where(np.abs(wave-10.5)<=2.5)
#dim=np.size(ix)
#longu=np.zeros(dim,dtype=np.float)
#coefHigh=np.zeros(4,dtype=np.float)
#coefHigh[0]=-8.73819754e-02
#coefHigh[1]=3.63087513e+00
#coefHigh[2]=-5.56164682e+01
#p = np.poly1d(coefPol[5,:])
#for i in range(dim):
#    coefHigh[3]=6.06546595e+02
#    coefHigh[3]=coefHigh[3]-p(wave[ix[0][i]])
#    longu[i]=np.float(np.roots(coefHigh)[2])
#coef=np.polyfit(385+ix[0][:],longu[:],2)
#print coef

#plt.show()
