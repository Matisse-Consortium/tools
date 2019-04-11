#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 23:06 2018

@author: fmillour
"""

import os
import wx
import sys
from astropy.io import fits
from mat_fileDialog import mat_FileDialog
from shutil import copyfile


# Function to remove all spaces from a given string
def removeSpaces(string):
 
    # To keep track of non-space character count
    count = 0
 
    list = []
 
    # Traverse the given string. If current character
    # is not space, then place it at index 'count++'
    for i in xrange(len(string)):
        if string[i] != ' ':
            list.append(string[i])
 
    return toString(list)

def change_calibmap_name(calibmap):
    direc = os.path.dirname(calibmap)
    
    # First check this is indeed a calib map file
    try:
        hdu = fits.getheader(calibmap)
        if hdu['HIERARCH ESO PRO CATG'] == 'OBS_FLATFIELD':
            print('yay let\'s do it!')
            print(hdu['HIERARCH ESO DET SEQ1 DIT'])
            if hdu['HIERARCH ESO DET CHIP NAME'] == 'HAWAII-2RG':
                # PIL
                obstype = hdu['HIERARCH ESO INS PIL NAME']
                resol   = hdu['HIERARCH ESO INS DIL NAME']
            elif hdu['HIERARCH ESO DET CHIP NAME'] == 'AQUARIUS':
                obstype = hdu['HIERARCH ESO INS PIN NAME']
                resol   = hdu['HIERARCH ESO INS DIN NAME']
                
                
            newName = os.path.join(direc,
                                hdu['DATE-OBS'].replace(':','_') +
                                '_' + hdu['HIERARCH ESO PRO CATG'] +
                                '_' + hdu['HIERARCH ESO DET CHIP NAME'] +
                                '_' + hdu['HIERARCH ESO DET READ CURNAME'] +
                                '_DC9_' + str(hdu['HIERARCH ESO DET CLDC1 DC9']) +
                                '_DIT_' + str(hdu['HIERARCH ESO DET SEQ1 DIT']) +
                                '_' + obstype + '_' + resol +
                                   '.fits');
            
                                   
            print("renaming "+calibmap+" into " +newName)
            #copyfile(src, dst)
            os.rename(calibmap, newName)
        else:
            print("Not a calibmap file!")
    except:
        print("Not a fits file!")
    
if __name__ == '__main__':
    # Get command line arguments
    listArg   = sys.argv
    name_file = []
    for elt in listArg:
        if ('--help' in elt):
            print( "Usage: Tidyup_calibmap.py <directory>")
            sys.exit(0)
        elif len(listArg) == 2:
            name_file = sys.argv[1]
            print(name_file)
        
    app = wx.App()
    if not name_file:
        print("No input name given, running file selector...")
        openFileDialog = mat_FileDialog(None, 'Open a file',"lmk,")
        if openFileDialog.ShowModal() == wx.ID_OK:
            name_file = openFileDialog.GetPaths()[0]
            print( name_file)
        openFileDialog.Destroy()
    app.MainLoop()
    app.Destroy()
    
    if os.path.isfile(name_file):
        print("Reading file "+name_file+"...")
        dic = change_calibmapFile_name(name_file)
    
    elif os.path.isdir(name_file):
        newdir = os.path.join(name_file,"CALIBMAP");
        try:
            print(newdir+" already exists...")
            os.stat(newdir)
        except:
            print("Creating directory "+newdir)
            os.mkdir(newdir)
        print("Listing files in directory")
        
        print(os.path.join(name_file,"**/*.fits*"))
        for root,subfolders,files in os.walk(name_file):
            for fileToTreat in files:
                #print(fileToTreat)
                if (fileToTreat.startswith('OBS_FLATFIELD') and (fileToTreat.endswith('fits')) ):
                    try:
                        hdu    = fits.getheader(os.path.join(root,fileToTreat))
                        if hdu['HIERARCH ESO PRO CATG'] == 'OBS_FLATFIELD':
                            print('Found an calibmap file. Copying it...')
                            #print(fileToTreat)
                            fil = os.path.basename(fileToTreat)
                            
                            copyfile(os.path.join(root,fileToTreat),
                                     os.path.join(newdir,fil))
                            change_calibmap_name(os.path.join(newdir,fileToTreat))
                    except:
                        print("Not a proper fits file!")
        
print("I made my job, baby!")
