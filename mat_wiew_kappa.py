# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 09:59:00 2016

@author: pbe
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages

hdulist = fits.open('C:\Users\pbe\Desktop\Python-TestPlan\KAPPA_N_LOW.fits')

nBeam=4
regName=[]
szSpectral=[]
coefKappa=[]
cornerSpectral=[]
for i in range(nBeam):
    iReg=hdulist['KAPPA_MATRIX'].data['REGION'][i]
    ix=np.where(hdulist['IMAGING_DETECTOR'].data['REGION'][:] == iReg)
    regName.append(hdulist['IMAGING_DETECTOR'].data['REGNAME'][ix[0]][0])
    szSpectral.append(hdulist['IMAGING_DETECTOR'].data['NAXIS'][ix[0]][0][1])
    coefKappa.append(hdulist['KAPPA_MATRIX'].data['MATRIX'][i])
    cornerSpectral.append(hdulist['IMAGING_DETECTOR'].data['CORNER'][ix[0]][0][1])
    
detName=hdulist[0].header['HIERARCH ESO DET CHIP NAME']
curName=hdulist[0].header['HIERARCH ESO DET READ CURNAME']
if (detName == 'AQUARIUS'):
    specResolution=hdulist[0].header['HIERARCH ESO INS DIN ID']  
    filterType=hdulist[0].header['HIERARCH ESO INS FIN ID']  
    polarType=hdulist[0].header['HIERARCH ESO INS PON ID']  
else:
    specResolution=hdulist[0].header['HIERARCH ESO INS DIL ID']  
    filterType=hdulist[0].header['HIERARCH ESO INS FIL ID']  
    polarType=hdulist[0].header['HIERARCH ESO INS POL ID']  
    
dispCoef0=hdulist[0].header['HIERARCH PRO DISP COEF0']
dispCoef1=hdulist[0].header['HIERARCH PRO DISP COEF1']
dispCoef2=hdulist[0].header['HIERARCH PRO DISP COEF2']

x=np.arange(szSpectral[0])+cornerSpectral[0]
wavelength=dispCoef0+dispCoef1*x+dispCoef2*x*x
 
configType=detName+"/"+curName+" --- "+specResolution+"/"+filterType+"/"+polarType
hdulist.close()  
    
f, axarr = plt.subplots(3, sharex=True, figsize=(8,8))
f.suptitle('KAPPA MATRIX',fontsize=20)
for i in range(nBeam):
    axarr[0].plot(wavelength,coefKappa[i][0:szSpectral[i]],label=regName[i],marker=(4,i,0),linestyle='')
axarr[0].plot([wavelength[0],wavelength[szSpectral[0]-1]],[1.5,1.5],linestyle='--',color='black')
axarr[0].plot([wavelength[0],wavelength[szSpectral[0]-1]],[3,3],linestyle='--',color='black')
axarr[0].set_ylim([1,3.5])
axarr[0].set_xlim([wavelength[0],wavelength[szSpectral[0]-1]])
axarr[0].set_ylabel('Ratio')
axarr[0].set_title(configType,fontsize=12)

for i in range(nBeam):
    axarr[1].plot(wavelength,coefKappa[i][szSpectral[i]:2*szSpectral[i]],label=regName[i],marker=(4,i,0),linestyle='')
axarr[1].set_ylim([5,7])
axarr[0].set_xlim([wavelength[0],wavelength[szSpectral[0]-1]])
#axarr[1].legend(prop={'size':8},loc=2)
axarr[1].set_ylabel('Zoom')

for i in range(nBeam):
    axarr[2].plot(wavelength,coefKappa[i][2*szSpectral[i]:3*szSpectral[i]],label=regName[i],marker=(4,i,0),linestyle='')
axarr[2].set_ylim([-0.1,0.1])
axarr[0].set_xlim([wavelength[0],wavelength[szSpectral[0]-1]])
axarr[2].legend(prop={'size':8},loc=1)
axarr[2].set_ylabel('Shift')
axarr[2].set_xlabel('Wavelength ($\mu$m)')

pp=PdfPages('foo.pdf')
pp.savefig(f)
pp.close()