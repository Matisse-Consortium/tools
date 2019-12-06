# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:37:15 2017

@author: pbe
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


corner=35
dimSpectral=1960
dimSpatialInterf=625
dimSpatialPhoto=150

xcoef=np.zeros(3*dimSpectral,dtype=np.double)
xcoef[dimSpectral:2*dimSpectral]=1.

ycoefInterf=np.zeros(3*dimSpatialInterf,dtype=np.double)
ycoefInterf[dimSpatialInterf:2*dimSpatialInterf]=1.
ycoefPhoto1=np.zeros(3*dimSpatialPhoto,dtype=np.double)
ycoefPhoto1[dimSpatialPhoto:2*dimSpatialPhoto]=1.
ycoefPhoto1[0:dimSpatialPhoto]=-4.287257202132794
ycoefPhoto2=np.zeros(3*dimSpatialPhoto,dtype=np.double)
ycoefPhoto2[dimSpatialPhoto:2*dimSpatialPhoto]=1.
ycoefPhoto2[0:dimSpatialPhoto]=-1.722299663118945
ycoefPhoto3=np.zeros(3*dimSpatialPhoto,dtype=np.double)
ycoefPhoto3[dimSpatialPhoto:2*dimSpatialPhoto]=1.
ycoefPhoto3[0:dimSpatialPhoto]=-1.842854589888842
ycoefPhoto4=np.zeros(3*dimSpatialPhoto,dtype=np.double)
ycoefPhoto4[dimSpatialPhoto:2*dimSpatialPhoto]=1.
ycoefPhoto4[0:dimSpatialPhoto]=-3.7748156822590966

xerror=np.zeros(3*dimSpectral,dtype=np.double)
yerrorInterf=np.zeros(3*dimSpatialInterf,dtype=np.double)
yerrorPhoto=np.zeros(3*dimSpatialPhoto,dtype=np.double)

hdu = fits.open('/data/CalibMap_New_Nov19/SHIFT_MAP_HAWAII-2RG_HIGH.fits')
#hdu = fits.open('SHIFT_MAP_HAWAII-2RG_HIGH.fits')
hdu[0].header.set('HIERARCH ESO INS DIL ID','HIGH+')
hdu[0].header.set('HIERARCH ESO INS DIL NAME','HIGH+')
hdu[0].header.set('HIERARCH ESO INS FIL ID','M')
hdu[0].header.set('HIERARCH ESO INS FIL NAME','M')
#hdu[0].header.set('HIERARCH ESO INS FIL NAME','M')

hdu['IMAGING_DETECTOR'].data['CORNER'][:,1]=35
hdu['IMAGING_DETECTOR'].data['NAXIS'][:,1]=1950

del hdu['SHIFT_MAP']
names=np.array(['DETECTOR','MAPX_1','MAPX_2','MAPX_3','MAPX_4','MAPX_5','MAPY_1','MAPY_2','MAPY_3','MAPY_4','MAPY_5','DISP'])
c1 = fits.Column(name='DETECTOR', format='1J', array=[1,1])
c2 = fits.Column(name='MAPX_1', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c3 = fits.Column(name='MAPX_2', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c4 = fits.Column(name='MAPX_3', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c5 = fits.Column(name='MAPX_4', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c6 = fits.Column(name='MAPX_5', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c7 = fits.Column(name='MAPY_1', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto1,yerrorPhoto])
c8 = fits.Column(name='MAPY_2', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto2,yerrorPhoto])
c9 = fits.Column(name='MAPY_3', format=str(3*dimSpatialInterf)+'D', array=[ycoefInterf,yerrorInterf])
c10 = fits.Column(name='MAPY_4', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto3,yerrorPhoto])
c11 = fits.Column(name='MAPY_5', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto4,yerrorPhoto])

#c12 = fits.Column(name='DISP', format='5D', array=[[4.11497338,-2.51977632e-04,1.13811555e-08,0,0.],[0.,0.,0.,0.,0.]])
a=5.14533488
b=-2.94751526e-04
c=3.17555020e-09
alpha=0.198616
beta=-0.0851867
gamma=0.00923958
aN=a+alpha+a*beta+gamma*a**2
bN=b+b*beta+2*a*b*gamma
cN=c+c*beta+gamma*b**2+2*a*c*gamma
dN=2*b*c*gamma
eN=gamma*c**2
print [aN,bN,cN,dN,eN]
x=np.arange(2048)*1.
plt.plot(x,a+b*x+c*x**2)
plt.plot(x,aN+bN*x+cN*x**2+dN*x**3+eN*x**4)
plt.show()
c12 = fits.Column(name='DISP', format='5D', array=[[aN,bN,cN,dN,eN],[0.,0.,0.,0.,0.]])

tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12],name='SHIFT_MAP')
tbhdu.header['INSTRUME'] = ('MATISSE', '')
tbhdu.header['MJD-OBS'] = (58815.709, '')
tbhdu.header['DATE-OBS'] = ('2019-11-29T17:02:04.4615', '')
tbhdu.header['NDETECT'] = (1, '')
tbhdu.header['NREGION'] = (5, '')
tbhdu.header['MAXTEL'] = (4, '')
tbhdu.header['HIERARCH ESO DET DID'] = ('ESO-VLT-DIC.NGCDCS-287412', '')
tbhdu.header['HIERARCH ESO DET ID'] = ('MATISSE-H2RG', '')


hdu.append(tbhdu)

#hdu.writeto('SHIFT_MAP_HAWAII-2RG_HIGH+_L.fits')
hdu.writeto('SHIFT_MAP_HAWAII-2RG_HIGH+_M.fits')
