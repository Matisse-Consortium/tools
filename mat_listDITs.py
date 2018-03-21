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
                    dit = header['HIERARCH ESO DET SEQ1 DIT']
                except:
                    pass
                chip = header["HIERARCH ESO DET CHIP NAME"]
                if chip == "HAWAII-2RG":
                    DITS_LM = np.append(DITS_LM, dit)
                elif chip == "AQUARIUS":
                    DITS_N  = np.append(DITS_N, dit)
                
            
        print('')
        
print("DITs LM",np.unique(DITS_LM))
print("DITs N",np.unique(DITS_N))
