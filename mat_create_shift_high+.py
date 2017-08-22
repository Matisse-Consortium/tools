# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:37:15 2017

@author: pbe
"""

import numpy as np
from astropy.io import fits

dimSpectral=1960
dimSpatialInterf=625
dimSpatialPhoto=150

xcoef=np.zeros(3*dimSpectral,dtype=np.double)
xcoef[dimSpectral:2*dimSpectral]=1.

ycoefInterf=np.zeros(3*dimSpatialInterf,dtype=np.double)
ycoefInterf[dimSpatialInterf:2*dimSpatialInterf]=1.
ycoefPhoto=np.zeros(3*dimSpatialPhoto,dtype=np.double)
ycoefPhoto[dimSpatialPhoto:2*dimSpatialPhoto]=1.

xerror=np.zeros(3*dimSpectral,dtype=np.double)
yerrorInterf=np.zeros(3*dimSpatialInterf,dtype=np.double)
yerrorPhoto=np.zeros(3*dimSpatialPhoto,dtype=np.double)

hdu = fits.open('C:\Users\pbe\Desktop\Hawaii\shift.fits')
hdu[0].header.set('HIERARCH ESO INS DIL ID','HIGH+')
hdu[0].header.set('HIERARCH ESO INS DIL NAME','HIGH+')
#hdu['IMAGING_DETECTOR'].data['NAXIS'][:,1]=dimSpectral
#hdu['IMAGING_DETECTOR'].data['CORNER'][:,1]=17
del hdu['SHIFT_MAP']
names=np.array(['DETECTOR','MAPX_1','MAPX_2','MAPX_3','MAPX_4','MAPX_5','MAPY_1','MAPY_2','MAPY_3','MAPY_4','MAPY_5','DISP'])
c1 = fits.Column(name='DETECTOR', format='1J', array=[1,1])
c2 = fits.Column(name='MAPX_1', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c3 = fits.Column(name='MAPX_2', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c4 = fits.Column(name='MAPX_3', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c5 = fits.Column(name='MAPX_4', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c6 = fits.Column(name='MAPX_5', format=str(3*dimSpectral)+'D', array=[xcoef,xerror])
c7 = fits.Column(name='MAPY_1', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto,yerrorPhoto])
c8 = fits.Column(name='MAPY_2', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto,yerrorPhoto])
c9 = fits.Column(name='MAPY_3', format=str(3*dimSpatialInterf)+'D', array=[ycoefInterf,yerrorInterf])
c10 = fits.Column(name='MAPY_4', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto,yerrorPhoto])
c11 = fits.Column(name='MAPY_5', format=str(3*dimSpatialPhoto)+'D', array=[ycoefPhoto,yerrorPhoto])
c12 = fits.Column(name='DISP', format='3D', array=[[0.,0.,0.],[0.,0.,0.]])

tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12],name='SHIFT_MAP')
tbhdu.header['INSTRUME'] = ('MATISSE', '')
tbhdu.header['MJD-OBS'] = (57813.709, '')
tbhdu.header['DATE-OBS'] = ('2017-03-01T17:02:04.4615', '')
tbhdu.header['NDETECT'] = (1, '')
tbhdu.header['NREGION'] = (5, '')
tbhdu.header['MAXTEL'] = (4, '')
tbhdu.header['HIERARCH ESO DET DID'] = ('ESO-VLT-DIC.NGCDCS-287412', '')
tbhdu.header['HIERARCH ESO DET ID'] = ('MATISSE-H2RG', '')


hdu.append(tbhdu)

hdu.writeto('C:\Users\pbe\Desktop\Hawaii\SHIFT_HIGH+.fits')
