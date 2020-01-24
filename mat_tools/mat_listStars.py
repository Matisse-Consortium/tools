#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

@author: fmillour

What star is present in a directory?

This software is a computer program whose purpose is to show oifits
files from the MATISSE instrument.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software. 

You can use, modify and/ or redistribute the software under the
terms of the CeCILL license as circulated by CEA, CNRS and INRIA at
the following URL "http://www.cecill.info". You have a copy of the
licence in the LICENCE.md file.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""

import sys
import os
import numpy as np
from astropy.io import fits
import fnmatch
from libAutoPipeline import matisseType

RES = []
obsTypes = ["TARGET_RAW", "CALIB_RAW","TARGET_RAW_INT", "CALIB_RAW_INT"]

for root, dirs, files in os.walk("."):
	print("Looking for observing files in "+root)
        for files in fnmatch.filter(files, "*.fits*"):
        
            header  = fits.open(os.path.join(root,files))[0].header
            types = matisseType(header)
            
            for typ in obsTypes:
                if types == typ:
          	    print(files)
                    try:
                        starname  = header['HIERARCH ESO OBS TARG NAME']
                        dit  = header['HIERARCH ESO DET SEQ1 DIT']
                        chip = header["HIERARCH ESO DET CHIP NAME"]
                        mod  = header["HIERARCH ESO DET READ CURNAME"]
                    except:
                        print("problem reading file")
			continue
                    if chip == "HAWAII-2RG":
                        try:
                            res     = header["HIERARCH ESO INS DIL NAME"]
                            strin = starname+" L-"+res+" "+types
                            RES   = np.append(RES, strin)
                            #print(strin)
                        except:
                            print("problem reading L file")
                    elif chip == "AQUARIUS":
                        try:
                            res     = header["HIERARCH ESO INS DIN NAME"]
                            strin = starname+" N-"+res+" "+types
                            RES   = np.append(RES, strin)
                            #print(strin)
                        except:
                            print("problem reading N file")
                        

res, idx = np.unique(RES, return_index=True)
print("Stars:")
for i in res:
	print(i)

