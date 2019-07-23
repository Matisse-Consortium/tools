#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from astropy.io import fits
import fnmatch


listArg  = sys.argv
name_dir = sys.argv[1]

if len(listArg) > 2:
    keys     = sys.argv[2:]
else:
    keys = ['DATE-OBS',
            'HIERARCH ESO DET SEQ1 DIT',
            'HIERARCH ESO DET NDIT'
            ]

files=os.listdir(name_dir)

print(('File name'), end=' ')
print(keys)

for file in files:
    if fnmatch.fnmatch(file,"*.fits*"):
       print((file), end=' ')
       hdu  = fits.open(os.path.join(name_dir,file))
       for ky in keys:
          vals = hdu[0].header[ky]
          print((vals), end=' ')
       print('')
