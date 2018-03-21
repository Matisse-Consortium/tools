#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits
import fnmatch


listArg  = sys.argv
name_dir = sys.argv[1]

files = os.listdir(name_dir)

DITS = []


for file in files:
    if fnmatch.fnmatch(file,"*.fits*"):
        print(file),
        
        header  = fits.open(os.path.join(name_dir,file))[0].header
        
        res  = ""
        catg = None
        typ  = None
        tech = None
    
        try:
            catg=header['HIERARCH ESO PRO CATG']
        except:
            try:
                catg = header['HIERARCH ESO DPR CATG']
                typ  = header['HIERARCH ESO DPR TYPE']
                tech = header['HIERARCH ESO DPR TECH']
            except:
                pass

        print(catg),
        print(typ),
        print(tech),
            
        if((catg=="SCIENCE" or catg=="TEST") and typ=="OBJECT" and tech=="INTERFEROMETRY"):
            res="TARGET_RAW"
            print(res),
            DITS = np.append(DITS, header['HIERARCH ESO DET SEQ1 DIT'])
        elif(catg == "TEST" and typ == "STD" and tech == "INTERFEROMETRY") or (catg == "CALIB" and typ == "OBJECT" and tech == "INTERFEROMETRY") or (catg == "CALIB" and typ == "OBJECT,FLUX" and tech == "INTERFEROMETRY") or (catg == "CALIB" and typ == "STD" and tech == "INTERFEROMETRY"): 
            res="CALIB_RAW"
            print(res),
            DITS = np.append(DITS, header['HIERARCH ESO DET SEQ1 DIT'])
        elif(catg == "TEST" or catg=="CALIB" or catg=="SCIENCE") and typ == "SKY" and tech == "INTERFEROMETRY": 
            res="SKY_RAW"
            print(res),
            DITS = np.append(DITS, header['HIERARCH ESO DET SEQ1 DIT'])
        else:
            res=catg
            
        print('')
        
print(np.unique(DITS))
