#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on 2019 June 24th
@author: F. Millour

Redo the data flagging properly (apply it to oifits files)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.

You can use, modify and/ or redistribute the software under the terms
of the CeCILL license as circulated by CEA, CNRS and INRIA at the
following URL "http://www.cecill.info". You have a copy of the licence
in the LICENCE.md file.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""

# Import stuff
import os
import argparse
from fnmatch import fnmatch
from tqdm import tqdm
import astropy
from astropy.io import fits
import numpy as np
from shutil import copyfile
from os import walk
import glob

###############################################################################

def reflagData(oifits, keepOldFlags=0):
    #print(oifits)
    
    data=fits.open(oifits, mode='update')

    #print("reading data")
    WLEN = data['OI_WAVELENGTH'].data['EFF_WAVE']

    # Read VIS2 data
    VIS2     = data['OI_VIS2'].data['VIS2DATA']
    VIS2ERR  = data['OI_VIS2'].data['VIS2ERR']
    VIS2FLAG = data['OI_VIS2'].data['FLAG']

    #print("Computing visibility")
    vis    = np.sqrt(np.abs(VIS2)) * np.sign(VIS2)
    viserr = 0.5* VIS2ERR / (vis+(vis==0));

    wlmin = (2.94e-6, 4.5e-6, 8e-6);
    wlmax = (4.1e-6,  5.0e-6, 12.5e-6);

    # print("Getting old vis flags")
    # if(keepOldFlags==1):
    #flag = ~bool(VIS2FLAG);
    # else:
    #     flag = np.array(1,np.shape(VIS2FLAG));
    
    #print("Computing new vis flags")
    flag = (VIS2 > 0. - VIS2ERR) &\
        (VIS2 < 1. + VIS2ERR) &\
        (VIS2 > -0.1)         &\
        (VIS2 < 1.1)          &\
        (VIS2ERR > 0)         &\
        (VIS2ERR < 0.1)       &\
        (vis > 0. - viserr)   &\
        (vis < 1. + viserr)   &\
        (vis > -0.1)          &\
        (vis < 1.1)           &\
        (viserr > 0)          &\
        (viserr < 0.1)
           
    #print("Filtering wavelength")
    # Filter interbands
    wlmask = np.zeros(np.shape(WLEN),dtype=bool)
    for i,wl in enumerate(wlmin):
        wlmaski = ((WLEN < wlmax[i]) &\
                   (WLEN > wlmin[i]))
        wlmask = wlmask | wlmaski;
    flag = flag & wlmask;
        
    #print("Inverting flags")
    # print(flag)
    flag = ~flag

    #print("Setting flags")
    data['OI_VIS2'].data['FLAG'] = flag
    
    #print("Getting CP flags")
    # Read CP data
    CP     = data['OI_T3'].data['T3PHI']
    CPERR  = data['OI_T3'].data['T3PHIERR']
    CPFLAG = data['OI_T3'].data['FLAG']

    #if(keepOldFlags==1):
    #    flag3 = ~CPFLAG;
    #else:
    #    flag3 = 1;
    
    #print("Computing new flags")
    flag3 = (CPERR > 0.) & (CPERR < 60)
    flag3 = flag3 & wlmask;

    #print("Redundant closure phase filtering")
    # Check the redundant closure phase
    fot i,cpi in enumerate()
           
    #print("Inverting flag")
    # print(flag)
    flag3 = ~flag3;

    data['OI_T3'].data['FLAG'] = flag3;
    
    #print("Saving fiel")
    data.flush();
    data.close();
    

###############################################################################

if __name__ == '__main__':
    print("Starting mat_reflagData...")

    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='A small utility tool to reflag MATISSE data.')

    #--------------------------------------------------------------------------
    parser.add_argument('oiFits', default="",  \
    help='The path to the file or directory containing your oifits data.')


    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_reflagData.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        print("\n     Example : python mat_tidyupOiFits.py /data/2018-05-19_OIFITS")
        sys.exit(0)
        
    if os.path.isfile(args.oiFits):
        print("Reading file "+args.oiFits+"...")
        reflagData(args.oiFits)

    elif os.path.isdir(args.oiFits):
        cwd = os.getcwd()

        newdir = os.path.join(cwd, os.path.basename(os.path.abspath(args.oiFits))+"_REFLAGGED");
        
        try:
            print(newdir+" already exists...")
            os.stat(newdir)
        except:
            print("Creating directory "+newdir)
            os.mkdir(newdir)
            
        print("Oifits files will be copied and reflagged in that directory.")
        
        filecounter = 0
        allfiles    = [];
        for root,subfolders,files in os.walk(args.oiFits):
            for fil in files:
                allfiles.append(os.path.join(root,fil))
                filecounter += 1;

        print("Number of files to treat:",len(allfiles))
        
        for fil in tqdm(allfiles, total=filecounter, unit=" files", unit_scale=False, desc="Working on files"):
            matchfilestoavoid = ["TARGET_CAL_0*","OBJ_CORR_FLUX_0*",
                                 "OI_OPDWVPO_*","PHOT_BEAMS_*",
                                 "CALIB_CAL_0*","RAW_DPHASE_*",
                                 "matis_eop*","nrjReal*",
                                 "DSPtarget*","nrjImag*",
                                 "fringePeak*","BSreal*","BSimag*",
                                 "*Image0.fits","*LFF.fits"]
            
            fifil = os.path.basename(fil)
            
            if fifil.endswith('fits'):
                #print(fifil)
                for i in matchfilestoavoid:
                    if fnmatch(fifil, i):
                        print("Breeaaaaaak!")
                        break;
                
                try:
                    #print("copying file to the right path...")
                    copyfile(fil,
                             os.path.join(newdir,fifil))
                    #print("reflagging...",os.path.join(newdir,fifil))
                    reflagData(os.path.join(newdir,fifil))
                    #print("Done!");
                except:
                    do_nothing = 1;
                    #print("Not a fits file!")

    print("I made my job!")
