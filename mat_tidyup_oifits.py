#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:38:46 2018

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

def change_oifitsFile_name(oifits):
    direc = os.path.dirname(oifits)
    
    # First check this is indeed an oifits file
    try:
        hdu     = fits.getheader(oifits)
        if hdu['HIERARCH ESO PRO CATG'] == 'CALIB_RAW_INT' or hdu['HIERARCH ESO PRO CATG'] == 'TARGET_RAW_INT':
            print('yay let\'s do it!')
            
            try:
                targ = hdu['HIERARCH ESO OBS TARG NAME']
            except:
                targ = hdu['OBJECT']
                    
            newName = os.path.join(direc,
                                hdu['DATE-OBS'].replace(':','_') +
                                '_' + targ.replace(" ","") +
                                '_' + hdu['HIERARCH ESO DET CHIP TYPE'] + '.fits');
                                   
            print("renaming "+oifits+" into " +newName)
            #copyfile(src, dst)
            os.rename(oifits, newName)
    except:
        print("Not a fits file!")
    
if __name__ == '__main__':
    # Get command line arguments
    listArg   = sys.argv
    name_file = []
    for elt in listArg:
        if ('--help' in elt):
            print( "Usage: mat_show_rawdata.py [--dir=start directory]")
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
        dic = change_oifitsFile_name(name_file)
    
    elif os.path.isdir(name_file):
        cwd = os.getcwd()
        newdir = os.path.join(cwd, os.path.basename(os.path.abspath(name_file))+"_OIFITS");
        try:
            print(newdir+" already exists...")
            os.stat(newdir)
        except:
            print("Creating directory "+newdir)
            os.mkdir(newdir)
        print("Listing files in directory")
        
        print(os.path.join(name_file,"**/*.fits*"))
        for root,subfolders,files in os.walk(name_file):
            for file in files:
                #print(file)
                if ( (file.startswith('CALIB_RAW_INT') or file.startswith('TARGET_RAW_INT')) and (file.endswith('fits')) ):
                    try:
                        hdu    = fits.getheader(os.path.join(root,file))
                        if hdu['HIERARCH ESO PRO CATG'] == 'CALIB_RAW_INT' or hdu['HIERARCH ESO PRO CATG'] == 'TARGET_RAW_INT':
                            print('Found an oifits file. Copying it...')
                            print(file)
                            fil = os.path.basename(file)
                            
                            copyfile(os.path.join(root,file),
                                     os.path.join(newdir,fil))
                            change_oifitsFile_name(os.path.join(newdir,file))
                    except:
                        print("Not a fits file!")
        
print("I made my job, baby!")
