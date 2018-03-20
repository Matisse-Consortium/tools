#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
  $Id$

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur
  
  Created on Tue Dec 20 13:22:21 2016
  @author: asoulain

  This software is a computer program whose purpose is to show raw data
  files from the MATISSE instrument.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
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
from mat_fileDialog import identifyFile
import wx
import fnmatch
import tempfile
from PIL import Image
from matplotlib.figure import Figure
from numpy import arange, sin, pi
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

#from skimage import data, img_as_float
#from skimage import exposure

# Save the file content in a temporary python binary file
_dBfile_matisse = tempfile.gettempdir()+'/data_matisse.dpy'
    
########################################################################
class displayRawData(wx.Frame):
    #----------------------------------------------------------------------
    def __init__(self,guiTitle="RAW DATA GUI",mainpath=".",fileTypes=[],checkPresent=[]):
        wx.Frame.__init__(self, None, wx.ID_ANY, guiTitle)
        self.panel = wx.Panel(self, wx.ID_ANY)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.AddSpacer(10)
        
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        t = arange(0.0, 3.0, 0.01)
        s = sin(2 * pi * t)

        self.axes.plot(t, s)
        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()
        
        # Recipe execution and SOF file management
        grid=wx.GridBagSizer(2,4)        
        btnChop       = wx.Button    (self.panel, label='Chop Frames',size=(100, -1))
        grid.Add(btnChop, (0, 0))    
        
        self.vbox.Add(grid,border=10,flag=wx.LEFT|wx.RIGHT|wx.EXPAND) 
        self.vbox.AddSpacer(20) 
    
    
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
    #tmp = readDb_matisse(name_file)
    #if not tmp is None:
    #    print (' > getData: data already in database MATISSE '+_dBfile_matisse)
    #    return tmp
    hdu      = pyfits.open(name_file)
    img_det  = hdu['IMAGING_DETECTOR']
    img_data = hdu['IMAGING_DATA'].data
        
    # read region names in imaging_detector and select photometries and interferometries only
    region_name = img_det.data['REGNAME']
    pht       = []
    pht_name  = []
    intf      = []
    intf_name = []
    img       = []
    img_name  = []
    for i,item in enumerate(region_name):
        if fnmatch.fnmatch(item, '*PHOT*'):
            pht_name.append(item)
            pht.append(img_data.field(i+1))
        if fnmatch.fnmatch(item, '*INTERF*'):
            intf_name.append(item)
            intf.append(img_data.field(i+1))
        if fnmatch.fnmatch(item, '*IMG*'):
            img_name.append(item)
            img.append(img_data.field(i+1))
            
    tartyp = img_data.field('TARTYP')
        
    dic = {'INTERF':intf, 'PHOT':pht, 'IMG':img, 'TARTYP':tartyp, 
           'INTERF_name':intf_name, 'PHOT_name':pht_name, 'IMG_name':img_name}
    #writeDb_matisse(name_file, dic)
    return dic

def plotSomestuff(dic,Grid,i0,key='INTERF'):
    data      = dic[key];
    data_name = dic[key+'_name'];
    for i in data:
        nPix = np.shape(i[0])[1]
    
    for i,it in enumerate(data):
        # get image histogram
        axes_his = plt.subplot(Grid[0,i+i0])
        #print(int(np.sqrt(len(it[0].flatten()))))
        image_histogram, bins = np.histogram(it[0].flatten(),
                                             int(np.sqrt(len(it[0].flatten()))), normed=False)
        axes_his.plot(bins[1:],image_histogram - np.median(image_histogram))
        plt.ylim([0,10])
        maximg = image_histogram - np.median(image_histogram) > 15
        maxim = np.max((bins[1:])[maximg])
        minim = np.median(it[0].flatten())
        #print(minim)
        #print(maxim)
        
        # Plot image with correct luminance
        axes_i = plt.subplot(Grid[1,i+i0])
        plt.xticks([]), plt.yticks([])
        axes_i.imshow(it[0], interpolation = 'nearest', cmap = 'afmhot', origin = 'down',
                      vmin=minim, vmax =maxim )
        axes_i.set_title(data_name[i])
    

def show_mat(dic):
    # Read the content of the dictionary
    interf = dic['INTERF'];
    phot   = dic['PHOT']
    img    = dic['IMG']
    tartyp = dic['TARTYP']
    
    nPix = [];
    for i in interf:
        nPix = np.append(nPix,np.shape(i[0])[1])
    for i in phot:
        nPix = np.append(nPix,np.shape(i[0])[1])
    for i in img:
        nPix = np.append(nPix,np.shape(i[0])[1])
       
    nFrames = np.shape(tartyp)[0]
    #print(nFrames)
    #print('titi')
        
    for i in interf:
        nFrames = np.shape(i)[0]
    for i in phot:
        nFrames = np.shape(i)[0]
    for i in img:
        nFrames = np.shape(i)[0]
    #print('toto')
    #print(nFrames)
        
    plt.close('all')
    plt.figure(figsize = (9, 6))
    G = gridspec.GridSpec(3, len(interf)+len(phot)+len(img),width_ratios=nPix,height_ratios=[2,10,1])
        
    plotSomestuff(dic,G,0,                    key='INTERF')
    plotSomestuff(dic,G,len(interf),          key='PHOT')
    plotSomestuff(dic,G,len(interf)+len(phot),key='IMG')
    
    print(nFrames)
    
    img = Image.new('RGB', (nFrames, 1), (0, 255, 0))
    for i,t in enumerate(tartyp):
        if t=='T':
            img.putpixel((i,0), (0,255,0))
        elif t=='U':
            img.putpixel((i,0), (255,0,0))
        elif t=='S':
            img.putpixel((i,0), (0,255,255))
            
    axes_c = plt.subplot(G[2,:])
    plt.imshow(img,aspect='auto')
    
    
    plt.xticks([]), plt.yticks([])
    
    plt.tight_layout()
    G.update(wspace=0.1)
    plt.show()

###############################################################################

if __name__ == '__main__':
    listArg = sys.argv
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
    
    print("Reading file "+name_file+"...")
    dic = open_mat(name_file)
    print("Plotting data "+name_file+"...")
    show_mat(dic)


