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
    hdu      = fits.open(oi_file)
    wl   = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    dic = {'WLEN':wl}
    
    try:
        dic['CFLUX']    = hdu['OI_VIS'].data['VISAMP'];
        dic['CFLUXERR'] = hdu['OI_VIS'].data['VISAMPERR'];
        dic['DPHI']     = hdu['OI_VIS'].data['VISPHI']
        dic['DPHIERR']  = hdu['OI_VIS'].data['VISPHIERR']
        dic['DPHI_U']   = hdu['OI_VIS'].data['UCOORD']
        dic['DPHI_V']   = hdu['OI_VIS'].data['VCOORD']
    except:  
        print("WARNING: No OI_VIS table!")
        
    try:
        dic['VIS2']    = hdu['OI_VIS2'].data['VIS2DATA']
        dic['VIS2ERR'] = hdu['OI_VIS2'].data['VIS2ERR']
        dic['VIS2_U']  = hdu['OI_VIS2'].data['UCOORD']
        dic['VIS2_V']  = hdu['OI_VIS2'].data['VCOORD']
    except:    
        print("WARNING: No OI_VIS2 table!")
        
    try:    
        dic['T3AMP']    = hdu['OI_T3'].data['T3AMP']
        dic['T3AMPERR'] = hdu['OI_T3'].data['T3AMPERR']
        dic['CLOS']     = hdu['OI_T3'].data['T3PHI']
        dic['CLOSERR']  = hdu['OI_T3'].data['T3PHIERR']
        dic['U1']       = hdu['OI_T3'].data['U1COORD']
        dic['V1']       = hdu['OI_T3'].data['V1COORD']
        dic['U2']       = hdu['OI_T3'].data['U2COORD']
        dic['V2']       = hdu['OI_T3'].data['V2COORD']
    except:
        print("WARNING: No OI_T3 table!")
        
    return dic

def show_oi_vs_freq(dic):
    wl    = dic['WLEN'];
    vis2  = dic['VIS2'];
    vis2e = dic['VIS2ERR'];
    u     = dic['VIS2_U'];
    v     = dic['VIS2_V'];
    cp    = dic['CLOS'];
    cpe   = dic['CLOSERR'];
    u1    = dic['U1'];
    v1    = dic['V1'];
    u2    = dic['U2'];
    v2    = dic['V2'];
    
    plt.figure(figsize = (9, 6))
    G = gridspec.GridSpec(2, 1)
    
    axes_v2 = plt.subplot(G[0,:])
    for i,j in enumerate(u):
        r = np.sqrt(u[i]**2 + v[i]**2);
        freq = r/wl;
        
        test = np.logical_and(vis2[i,:] >= 0, vis2e[i,:] / vis2[i,:] < 1)
        
        #plt.plot(freq, vis2[i,:])
        #plt.plot(freq[test], vis2[i,test])
        plt.semilogy(freq, vis2[i,:],color='lightgray')
        plt.semilogy(freq[test], vis2[i,test])
        plt.ylim([1e-4,1.1])
        axes_v2.set_title('Squared visibilities vs frequencies')
        
    axes_cp = plt.subplot(G[1,:])
    for i,j in enumerate(u1):
        r1 = np.sqrt(u1[i]**2 + v1[i]**2);
        r2 = np.sqrt(u2[i]**2 + v2[i]**2);
        r3 = np.sqrt((u1[i]+u2[i])**2 + (v1[i]+v2[i])**2);
        freq = np.maximum(np.maximum(r1,r2),r3)/wl;
        
        test = np.absolute(cpe[i,:]) < 180/pi/3
        
        #plt.plot(freq, vis2[i,:])
        #plt.plot(freq[test], vis2[i,test])
        plt.plot(freq, cp[i,:],color='lightgray')
        plt.plot(freq[test], cp[i,test])
        plt.ylim([-200,200])
        axes_cp.set_title('Closure phase vs frequencies')
        

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
    show_oi_vs_freq(dic)


