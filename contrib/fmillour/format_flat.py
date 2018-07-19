#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Jul 18th 2018

@author: fmillour
"""

import os
import wx
import sys
import numpy as np
from astropy.io import fits
from mat_fileDialog import mat_FileDialog
from shutil import copyfile
from   matplotlib import pyplot as plt


flatFile = "/data/CalibMap/FLATFIELD_AQUARIUS_HIGH_GAIN_DIT_0.075.fits";

header = fits.getheader(flatFile)

#print(header)

hdu   = fits.open(flatFile)
image = hdu[0].data;

plt.imshow(image)

newflat = image;

for i in range(32):
    gain = header["HIERARCH ESO QC DET2 CHANNEL"+str(i+1)+" GAIN2"]
   # print(gain)    
    newflat[512:1024,i*32:(i+1)*32] = gain * image[512:1024,i*32:(i+1)*32]

for i in range(32):
    gain = header["HIERARCH ESO QC DET2 CHANNEL"+str(i+33)+" GAIN2"]
   # print(gain)    
    newflat[0:512,i*32:(i+1)*32] = gain * image[0:512,i*32:(i+1)*32]


newflat = newflat / np.average(newflat[:,360:650])

plt.imshow(newflat)
plt.show()

hdu[0].data = newflat;

hdu.writeto("/data/newflat.fits")
