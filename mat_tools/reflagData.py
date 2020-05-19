#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 2019 June 24th

Redo the data flagging properly

@author: F. Millour
"""

# Import stuff
import astropy
from astropy.io import fits
import numpy as np
from shutil import copyfile
from os import walk
import glob

# Work on all FITS files in the current directory
dir="*.fits"

#print(glob.glob(dir))

for filenames in glob.glob(dir):
    print(filenames)
    
    data=fits.open(filenames, mode='update')

    WLEN = data['OI_WAVELENGTH'].data['EFF_WAVE']

    # Read VIS2 data
    VIS2     = data['OI_VIS2'].data['VIS2DATA']
    VIS2ERR  = data['OI_VIS2'].data['VIS2ERR']
    VIS2FLAG = data['OI_VIS2'].data['FLAG']

   #wlmin = 3.5e-6;
   #wlmax = 4.1e-6;
   #wlmin = 2.94e-6;
   #wlmax = 4.19e-6;
    wlmin = 8e-6;
    wlmax = 12.5e-6;
    
   #flag = ~ VIS2FLAG
    flag = (WLEN < wlmax)      &\
           (WLEN > wlmin)      &\
           (VIS2 > 0. - VIS2ERR) &\
           (VIS2 < 1. + VIS2ERR) &\
           (VIS2 > -0.1)         &\
           (VIS2 < 1.1)          &\
           (VIS2ERR > 0)         &\
           (VIS2ERR < 0.1)    #   &\
         # (VIS2 / VIS2ERR > 3)

    print(flag)
    flag = ~flag

    data['OI_VIS2'].data['FLAG'] = flag
    
    # Read CP data
    CP     = data['OI_T3'].data['T3PHI']
    CPERR  = data['OI_T3'].data['T3PHIERR']
    CPFLAG = data['OI_T3'].data['FLAG']

    flag = (WLEN < wlmax) &\
           (WLEN > wlmin) &\
           (CPERR > 0.)     &\
           (CPERR < 60)

    print(flag)
    flag = ~ flag

    data['OI_T3'].data['FLAG'] = flag
    
    data.flush()
    data.close()
