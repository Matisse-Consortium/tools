#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
  $Id: $

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur
  
  Created on Sat Mar 17 06:39:49 2018
  @author: fmillour
  fmillour@oca.eu

  This software is a computer program whose purpose is to show oifits
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


import sys
import wx
from mat_fileDialog import mat_FileDialog
from mat_fileDialog import identifyFile
from astropy.io import fits as fits
    
def open_oi(oi_file):
    #tmp = readDb_matisse(name_file)
    #if not tmp is None:
    #    print (' > getData: data already in database MATISSE '+_dBfile_matisse)
    #    return tmp
    hdu      = fits.open(oi_file)
    wl   = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    
    vis2 = hdu['OI_VIS2'].data['VIS2DATA']
    vis2e = hdu['OI_VIS2'].data['VIS2ERR']
    u = hdu['OI_VIS2'].data['UCOORD']
    v = hdu['OI_VIS2'].data['VCOORD']
    
    cp   = hdu['OI_T3'].data['T3PHI']
    u1 = hdu['OI_T3'].data['U1COORD']
    v1 = hdu['OI_T3'].data['V1COORD']
    u2 = hdu['OI_T3'].data['U2COORD']
    v2 = hdu['OI_T3'].data['V2COORD']
                
    dic = {'WLEN':wl, 'VIS2':vis2, 'VIS2ERR':vis2e, 'CP':cp, 'U':u, 'V':v, 'U1':u1, 'V1':v1, 'U2':u2, 'V2':v2}
    return dic

def show_oi(dic):
    wl = dic['WLEN'];
    vis2 = dic['VIS2'];
    vis2e = dic['VIS2ERR'];
    u = dic['U'];
    v = dic['V'];
    cp = dic['CP'];
    
    print(np.shape(vis2))
    
    for i,j in enumerate(u):
        r = np.sqrt(u[i]**2 + v[i]**2);
        freq = r/wl;
        
        test = np.logical_and(vis2[i,:] >= 0, vis2e[i,:] / vis2[i,:] < 0.5)
        print(vis2[i,:])
        #plt.plot(freq, vis2[i,:])
        plt.figure(1)
        #plt.plot(freq[test], vis2[i,test])
        plt.semilogy(freq[test], vis2[i,test])
        plt.ylim([-0.2,1.1])
        
        plt.figure(2)
        plt.plot(wl, vis2e[i,:])
        

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
    dic = open_oi(name_file)
    print("Plotting data "+name_file+"...")
    show_oi(dic)


