# -*- coding: utf-8 -*-
"""
Created on Tue Sep 06 09:04:52 2016

@author: pbe
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages

hdulist = fits.open('C:\Users\pbe\Desktop\Python-TestPlan\MATISSE_GEN_cal_shift_L_SPCAL_OPEN_HIGH_210_0001.fits')
nBeam=5
szSpectral=[]
cornerSpectral=[]
data=[]
for i in range(nBeam):
    szSpectral.append([hdulist['IMAGING_DETECTOR'].data['NAXIS'][i+8][0],hdulist['IMAGING_DETECTOR'].data['NAXIS'][i+8][1]])
    cornerSpectral.append([hdulist['IMAGING_DETECTOR'].data['CORNER'][i+8][0],hdulist['IMAGING_DETECTOR'].data['CORNER'][i+8][1]])

data.append(hdulist['IMAGING_DATA'].data['DATA9'][0])
data.append(hdulist['IMAGING_DATA'].data['DATA10'][0])
data.append(hdulist['IMAGING_DATA'].data['DATA11'][0])
data.append(hdulist['IMAGING_DATA'].data['DATA12'][0])
data.append(hdulist['IMAGING_DATA'].data['DATA13'][0])
hdulist.close()  


imgFullFrame=np.zeros((2048,2048),dtype=np.float)

for i in range(nBeam):
    imgFullFrame[cornerSpectral[i][1]-1:cornerSpectral[i][1]-1+szSpectral[i][1],cornerSpectral[i][0]-1:cornerSpectral[i][0]-1+szSpectral[i][0]]=data[i]
plt.figure(1)
plt.imshow(imgFullFrame[1000:1100,1000:1100])

hdulist = fits.open('C:\Users\pbe\Desktop\Python-TestPlan\matisse_bpm_det1_rid2.fits')
bpm=hdulist[0].data
hdulist.close()  
plt.figure(2)
plt.imshow(bpm[1000:1100,1000:1100])


hdulist = fits.open('C:\Users\pbe\Desktop\Python-TestPlan\imageraw.fits')
raw=hdulist[0].data
imgraw=np.zeros((2048,2048),dtype=np.float)
imgraw[41:41+1960,750:750+625]=raw[:,300:300+625]
hdulist.close()  
