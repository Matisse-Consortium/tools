#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits
import fnmatch
from libAutoPipeline import matisseType

RES_LM = []
RES_N = []
obsTypes = ["TARGET_RAW", "CALIB_RAW", "SKY_RAW"]

for root, dir, files in os.walk("."):
        print root
        print ""
        for file in fnmatch.filter(files, "*.fits*"):
            print(file),
        
            header  = fits.open(os.path.join(root,file))[0].header
            type = matisseType(header)
            
            for typ in obsTypes:
                if type == typ:
                    print(typ),
                    try:
                        dit  = header['HIERARCH ESO DET SEQ1 DIT']
                        chip = header["HIERARCH ESO DET CHIP NAME"]
                        mod  = header["HIERARCH ESO DET READ CURNAME"]
                    except:
                        pass
                    if chip == "HAWAII-2RG":
                        try:
                            res     = header["HIERARCH ESO INS DIL NAME"]
                            strin = chip+" "+res+" "+mod+" "+str(dit)
                            RES_LM   = np.append(RES_LM, strin)
                            print(strin)
                        except:
                            pass
                    elif chip == "AQUARIUS":
                        try:
                            res     = header["HIERARCH ESO INS DIN NAME"]
                            strin = chip+" "+res+" "+mod+" "+str(dit)
                            RES_N   = np.append(RES_N, strin)
                            print(strin)
                        except:
                            pass
                        
                        
            print(' ')

res, idx = np.unique(RES_LM, return_index=True)
print("DITs LM",res)
res, idx = np.unique(RES_N, return_index=True)
print("DITs N", res)
