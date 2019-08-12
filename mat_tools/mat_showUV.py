#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Thu Dec  6 22:22:19 2018
@author: fmillour

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.

You can use, modify and/ or redistribute the software under the terms
of the CeCILL license as circulated by CEA, CNRS and INRIA at the
following URL "http://www.cecill.info". You have a copy of the licence
in the LICENCE.md file.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""

import wx
import argparse
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from   astropy.io import fits as fits
from mat_fileDialog import mat_FileDialog, identifyFile

###############################################################################

def open_hdr(oi_file):
    try:
        hdu = fits.open(oi_file)
    except IOError:
        print(("Unable to read fits file: " + oi_file))
        return {}

    hdr = hdu[0].header

    #print(hdr)

    return hdr

###############################################################################

def get_UV(file):

    BX = []
    BY = []

    res = open_hdr(file)

    try:
        target_name = res['HIERARCH ESO OBS TARG NAME']
    except:
        print("No TARGET name")
        target_name = "-"

    typ    = res['HIERARCH ESO PRO CATG']

    try:
        instru = res['INSTRUME']
    except:
        try:
            instru = res['HIERARCH ESO INS MODE']
        except:
            print('error, unknown instrument!')

    #print(instru)


    if instru == 'MATISSE':
        ntels = 4;
    elif instru == 'MIDI':
        ntels = 2;
    else:
        ntels = res['HIERARCH ESO DET NTEL']

    #read in priority the keywords
    try:
        base=0;
        for i in np.arange(1,ntels+1):
            for j in np.arange(i+1,ntels+1):
                tel1 = i;
                tel2 = j;
                base+=1;

                blen = (res['HIERARCH ESO ISS PBL'+str(tel1)+str(tel2)+' START']+res['HIERARCH ESO ISS PBL'+str(tel1)+str(tel2)+' END'])/2
                bang = (res['HIERARCH ESO ISS PBLA'+str(tel1)+str(tel2)+' START']+res['HIERARCH ESO ISS PBLA'+str(tel1)+str(tel2)+' END'])/2

                bx = blen * np.sin(bang * np.pi / 180.)
                by = blen * np.cos(bang * np.pi / 180.)

                BX = np.append(BX, bx)
                BY = np.append(BY, by)
    # otherwise look into the OI_VIS2 table directly
    except:
        print('No PBL keywords. taking a look into the OI_VIS2 table...')
        hdu = fits.open(file)
        #print(hdu['OI_VIS2'])

    return BX,BY,target_name, typ

###############################################################################

def get_UVs(files):
    BX   = []
    BY   = []
    TARG = []

    for file in files:
        print(file)
        bx, by, targ, typ = get_UV(file)

        if typ == 'CALIB_RAW_INT':
            continue

        if targ == '-':
            continue

        BX = np.append(BX, bx, axis=0)
        BY = np.append(BY, by, axis=0)

        TARG = np.append(TARG, targ)

    return BX, BY, TARG;

###############################################################################

def plot_UV(BX, BY, TARG, marker='o', markersize=4, color="red"):

    uniques = set(TARG)
    print(uniques)
    if len(uniques) < 4:
        nwiny = len(uniques)
        nwinx = 1
    else:
        nwin = int(np.sqrt(len(uniques)))
        nwiny = nwin;
        nwinx = nwin+1;

    for i,base in enumerate(BX):
        for j,uni in enumerate(uniques):
            if uni == TARG[i/6]:
                idx = j

        #print(idx)
        try:
            ax = plt.subplot(nwinx,nwiny,idx+1)
            ax.axis('equal')
            aruni = np.array(list(uniques))
            #print(aruni)
            ax.set_title(aruni[idx])
            #print(aruni[idx])
            
            ax.plot( base, BY[i], marker=marker, markersize=markersize,   color=color)
            ax.plot(-base,-BY[i], marker='o',    markersize=markersize-2, color=color)
        except:
            print("ouille !")

    ax.plot(0,0, marker='+', markersize=markersize, color=color)

###############################################################################

if __name__ == '__main__':
    print("Starting...")
    
    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='(u,v) plane plot.')
    #--------------------------------------------------------------------------
    parser.add_argument('in_dir', metavar='in_dir', type=str, \
    help='The path to the directory containing your oifits data.', default=None)
    #--------------------------------------------------------------------------

    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_autoCalib.py --help to be kind with you:\033[0m\n")
        parser.print_help()
	sys.exit(0)

    ########################################################################

    files = glob.glob(args.in_dir+"/*.fits*")

    plt.figure(1)
    plt.clf();

    BX, BY, TARG = get_UVs(files)
    plot_UV(BX, BY, TARG)
    plt.show();
    
