#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
  $Id$

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
  
  Created on Tue Dec 20 13:22:21 2016

  @author: asoulain
"""

# Import necessary libraries
import matplotlib as mpl
mpl.use('TkAgg')
import sys
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
import os
import pickle as cPickle
import numpy as np
import matplotlib.gridspec as gridspec
from mat_fileDialog import mat_FileDialog
import wx
import fnmatch
import tempfile


# Save the file content in a temporary python binary file
_dBfile_matisse = tempfile.gettempdir()+'/data_matisse.dpy'
    
def readDb_matisse(name_file):
    global _dBfile_matisse
    if not os.path.exists(_dBfile_matisse):
        return None
    f = open(_dBfile_matisse, 'rb')
    db = cPickle.load(f, encoding="utf8")
    f.close()
    if name_file in db.keys():
        return db[name_file]
    else:
        #print db.keys()
        return None

def writeDb_matisse(name_file, data):
    global _dBfile_matisse
    if not os.path.exists(_dBfile_matisse):
        print (' > getData: creating', _dBfile_matisse)
        db = {name_file:data}
    else:
        print (' > getData: updating', _dBfile_matisse)
        f  = open(_dBfile_matisse, 'rb')
        db = cPickle.load(f, encoding="utf8")
        f.close()
    db[name_file] = data
    f = open(_dBfile_matisse, 'wb')
    cPickle.dump(db, f, 2)
    f.close()
    return
    
def open_mat(name_file):
    tmp = readDb_matisse(name_file)
    if not tmp is None:
        print (' > getData: data already in database MATISSE '+_dBfile_matisse)
        return tmp
    hdu      = pyfits.open(name_file)
    img_det  = hdu['IMAGING_DETECTOR']
    img_data = hdu['IMAGING_DATA'].data
        
    # read region names in imaging_detector and select photometries and interferometries only
    region_name = img_det.data['REGNAME']
    pht  = []
    intf = []
    for i,item in enumerate(region_name):
        if fnmatch.fnmatch(item, '*PHOT*'):
            pht.append(img_data.field(i))
        if fnmatch.fnmatch(item, '*INTERF*'):
            intf.append(img_data.field(i))
    print("tutu")
    print(len(region_name))
    print(region_name)
    
    #print(pht)
    #print(intf)
    
    # Select data with interferometries and photometries
    if intf:
        interf = intf
        # convert to plain regular numpy array
        interf2 = interf#.view(interf.dtype.fields or interf.dtype, np.ndarray)
        print(type(interf))
        #print(np.shape(interf))
        print(len(interf))
        print(type(interf2))
        #print(np.shape(interf2))
        print(len(interf2))
    else:
        interf2 = [];
    if pht:
        phot  = pht
        # convert to plain regular numpy array
        phot2 = phot#.view(phot.dtype.fields or phot.dtype, np.ndarray)
    else:
        phot2 = [];
        
    dic = {'INTERF':interf2, 'PHOT':phot2}
    writeDb_matisse(name_file, dic)
    return dic

def show_mat(dic):
    interf = dic['INTERF'];
    print("titi")
    print(type(interf))
    #print(np.shape(interf))
    print(len(interf))
    
    phot = dic['PHOT']
    print(type(phot))
    #print(np.shape(interf))
    print(len(phot))
    
    plt.close('all')
    plt.figure(figsize = (9, 6))
    G = gridspec.GridSpec(1, len(interf))
    
    
    for i,it in enumerate(interf):
        axes_i = plt.subplot(G[:,i])
        plt.xticks([]), plt.yticks([])
        axes_i.imshow(it[0], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_i.set_title('INTERF')
        
        
    for p,ph in enumerate(phot):
        axes_p = plt.subplot(G[:,i])
        axes_p.imshow(ph[0], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_p.set_title('PHOT')
    
    
    
    
    plt.xticks([]), plt.yticks([])
    
    plt.tight_layout()
    G.update(wspace=0.1)
    plt.show()

###############################################################################

listArg = sys.argv
name_file = []
for elt in listArg:
    if ('--help' in elt):
        print( "Usage: mat_show_rawdata.py [--dir=start directory]")
        sys.exit(0)
    elif len(listArg) == 2:
        name_file = sys.argv[1]
        print(name_file)
        

if __name__ == '__main__':

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
    print("Reading file "+name_file+"...")
    dic = open_mat(name_file)
    print("Plotting data "+name_file+"...")
    show_mat(dic)


