#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits
import fnmatch
from libAutoPipeline import matisseType


listArg  = sys.argv
name_dir = sys.argv[1]

files = os.listdir(name_dir)

DITS_LM = []
DITS_N  = []
RES_LM = []
RES_N  = []
obsTypes = ["TARGET_RAW", "CALIB_RAW", "SKY_RAW"]

for file in files:
    if fnmatch.fnmatch(file,"*.fits*"):
        print(file),
        
        header  = fits.open(os.path.join(name_dir,file))[0].header
        type = matisseType(header)

        for typ in obsTypes:
            if type == typ:
                print(typ),
                try:
                    dit  = header['HIERARCH ESO DET SEQ1 DIT']
                    chip = header["HIERARCH ESO DET CHIP NAME"]
                except:
                    pass
                if chip == "HAWAII-2RG":
                    try:
                        res     = header["HIERARCH ESO INS DIL NAME"]
                        RES_LM  = np.append(RES_LM, res)
                        DITS_LM = np.append(DITS_LM, dit)
                    except:
                        pass
                elif chip == "AQUARIUS":
                    try:
                        res     = header["HIERARCH ESO INS DIN NAME"]
                        RES_N   = np.append(RES_N, res)
                        DITS_N  = np.append(DITS_N, dit)
                    except:
                        pass
                
            
        print('')

res, idx = np.unique(DITS_LM, return_index=True)
print("DITs LM",res)
print("RES  LM",RES_LM[idx])
res, idx = np.unique(DITS_N, return_index=True)
print("DITs N", res)
print("RES  N", RES_N[idx])
