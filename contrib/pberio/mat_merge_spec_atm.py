# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:02:14 2017

@author: pbe
"""


import numpy as np
from astropy.io import fits
import glob

fic=glob.glob('.\\SPEC_ATM\\Spec_Atmo_*')

nbptmax=10000
arr1=[]
arr2=[]
arr3=[]
arr4=[]
arr5=[]
arr6=[]

#for i in range(0,1):
for i in range(0,len(fic)):
    band=fic[i].split("_")[3]
    if (band=="L"):
        detector="HAWAII-2RG"
    else:
        detector="AQUARIUS"
    resolution=fic[i].split("_")[4].upper()
    pwv=float(((fic[i].split("_")[5]).split("V")[1]).split(".fits")[0])
    hdu = fits.open(fic[i])
    wavelength=hdu[1].data["lam"]
    flux=hdu[1].data["flux"]
    nbpt=np.size(wavelength)
    del hdu
    print detector,resolution,pwv

    wave=np.zeros(nbptmax)
    spec=np.zeros(nbptmax)
    wave[0:nbpt]=wavelength
    spec[0:nbpt]=flux
    arr1.append(detector)
    arr2.append(resolution)
    arr3.append(pwv)
    arr4.append(nbpt)
    arr5.append(wave)
    arr6.append(spec)
        
    
c1 = fits.Column(name='DETECTOR', format='16A', array=arr1)
c2 = fits.Column(name='RESOLUTION', format='16A', array=arr2)
c3 = fits.Column(name='PWV', format='D', array=arr3)
c4 = fits.Column(name='NBPT', format='I', array=arr4)
c5 = fits.Column(name='lam', format=str(nbptmax)+'D', array=arr5)
c6 = fits.Column(name='flux', format=str(nbptmax)+'D', array=arr6)
tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6],name='SPEC_CALC')
tbhdu.writeto('test.fits')
del tbhdu
hdu = fits.open('test.fits')
hdu[0].header.set('HIERARCH ESO DPR CATG','SPEC_ATM')
hdu[0].header.set('HIERARCH ESO DPR TYPE','SKYCALC')
hdu[0].header.set('HIERARCH ESO DPR TECH','ESO')
hdu[0].header.set('DATE-OBS','2017-08-22T12:00:00')
hdu[0].header.set('MJD-OBS',56555.5)
hdu.writeto('test2.fits')


        