#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 22:22:19 2018

@author: fmillour
"""

import wx
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

def plot_UV(BX, BY, TARG, marker='o', markersize=6, color="red"):
    plt.axis('equal')

    uniques = set(TARG)
    print(uniques)
    nwin = int(np.sqrt(len(uniques)))

    for i,base in enumerate(BX):
        for j,uni in enumerate(uniques):
            if uni == TARG[i/6]:
                idx = j

        ax = plt.subplot(nwin,nwin+1,idx+1)
        plt.axis('equal')
        aruni = np.array(list(uniques))
        print(aruni)
        ax.set_title(aruni[idx])
        print(aruni[idx])

        plt.plot( base, BY[i], marker=marker, markersize=markersize,   color=color)
        plt.plot(-base,-BY[i], marker='o',    markersize=markersize-3, color=color)

    plt.plot(0,0, marker='+', markersize=markersize, color=color)

###############################################################################

if __name__ == '__main__':
    listArg = sys.argv
    for elt in listArg:
        if ('--help' in elt):
            print ("Usage: mat_fileDialog.py [--dir=start directory]")
            sys.exit(0)

    repBase = []
    for elt in listArg:
        if ('--dir' in elt):
            item=elt.split('=')
            repBase=item[1]

    print(repBase)

    app = wx.App()

    openFileDialog = mat_FileDialog(None, 'Select a directory','')
    if openFileDialog.ShowModal() == wx.ID_OK:
        directory = openFileDialog.GetPaths()
        print(directory)


###############################################################################

    files = glob.glob(directory[0]+"/*.fits*")

    plt.figure(1)
    plt.clf();

    BX, BY, TARG = get_UVs(files)
    plot_UV(BX, BY, TARG)
    plt.show();
    
    openFileDialog.Destroy()
    app.MainLoop()
    app.Destroy()
