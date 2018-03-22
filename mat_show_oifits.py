#!/usr/bin/env python
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
import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from mat_fileDialog import mat_FileDialog
from mat_fileDialog import identifyFile
from astropy.io import fits as fits
import os
import glob

###############################################################################

def open_oi(oi_file):
    try:
        hdu = fits.open(oi_file)
    except IOError:
        print "Unable to read fits file: "+oi_file
        return {}

    hdr = hdu[0].header

    wl  = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    dic = {'WLEN':wl}

    target_name = hdu['OI_TARGET'].data['TARGET'][0]
    if not target_name:
        try:
            target_name = hdr['HIERARCH ESO OBS TARG NAME']
        except KeyError:
            print "Target name not found."
            target_name = ""
    dic['TARGET'] = target_name
    try:
        dateobs = hdr['DATE-OBS']
    except KeyError:
        dateobs = ""
    dic['DATEOBS'] = dateobs
    try:
        det_name = hdr['HIERARCH ESO DET CHIP NAME']
    except KeyError:
        print "Detector name not found."
        det_name = ""
    if (det_name == 'AQUARIUS'):
        band = 'N'
    elif (det_name == 'HAWAII-2RG'):
        band = 'LM'
    else:
        band = ''
    dic['BAND'] = band
    try:
        dispersion_name = hdr['HIERARCH ESO DET DISP NAME']
    except KeyError:
        print "Dispersion name not found."
        dispersion_name = ""
    dic['DISP'] = dispersion_name
    try:
        DIT = hdr["HIERARCH ESO DET SEQ1 DIT"]  # (s)
    except KeyError:
        DIT = np.nan
        print "DIT not found"
    dic['DIT'] = DIT

    try:
        dic['VIS'] = {}
        dic['VIS']['CFLUX']    = hdu['OI_VIS'].data['VISAMP'];
        dic['VIS']['CFLUXERR'] = hdu['OI_VIS'].data['VISAMPERR'];
        dic['VIS']['DPHI']     = hdu['OI_VIS'].data['VISPHI']
        dic['VIS']['DPHIERR']  = hdu['OI_VIS'].data['VISPHIERR']
        dic['VIS']['U']        = hdu['OI_VIS'].data['UCOORD']
        dic['VIS']['V']        = hdu['OI_VIS'].data['VCOORD']
        dic['VIS']['TIME']     = hdu['OI_VIS'].data['MJD'];
    except:
        print("WARNING: No OI_VIS table!")

    try:
        dic['VIS2'] = {}
        dic['VIS2']['VIS2']    = hdu['OI_VIS2'].data['VIS2DATA']
        dic['VIS2']['VIS2ERR'] = hdu['OI_VIS2'].data['VIS2ERR']
        dic['VIS2']['U']       = hdu['OI_VIS2'].data['UCOORD']
        dic['VIS2']['V']       = hdu['OI_VIS2'].data['VCOORD']
        dic['VIS2']['TIME']    = hdu['OI_VIS2'].data['MJD'];
    except:
        print("WARNING: No OI_VIS2 table!")

    try:
        dic['TF2'] = {}
        dic['TF2']['VIS2']    = hdu['OI_TF2'].data['TF2']
        dic['TF2']['VIS2ERR'] = hdu['OI_TF2'].data['TF2ERR']
        #dic['TF2']['U']       = hdu['OI_TF2'].data['UCOORD']
        #dic['TF2']['V']       = hdu['OI_TF2'].data['VCOORD']
        dic['TF2']['TIME']    = hdu['OI_TF2'].data['MJD'];
    except:
        print("WARNING: No OI_TF2 table!")

    try:
        dic['T3'] = {}
        dic['T3']['T3AMP']    = hdu['OI_T3'].data['T3AMP']
        dic['T3']['T3AMPERR'] = hdu['OI_T3'].data['T3AMPERR']
        dic['T3']['CLOS']     = hdu['OI_T3'].data['T3PHI']
        dic['T3']['CLOSERR']  = hdu['OI_T3'].data['T3PHIERR']
        dic['T3']['U1']       = hdu['OI_T3'].data['U1COORD']
        dic['T3']['V1']       = hdu['OI_T3'].data['V1COORD']
        dic['T3']['U2']       = hdu['OI_T3'].data['U2COORD']
        dic['T3']['V2']       = hdu['OI_T3'].data['V2COORD']
        dic['T3']['TIME']     = hdu['OI_T3'].data['MJD'];
    except:
        print("WARNING: No OI_T3 table!")

    try:
        dic['FLUX'] = {}
        dic['FLUX']['FLUX']    = hdu['OI_FLUX'].data['FLUXDATA']
        dic['FLUX']['FLUXERR'] = hdu['OI_FLUX'].data['FLUXERR']
        dic['FLUX']['TIME']    = hdu['OI_FLUX'].data['MJD'];
    except:
        print("WARNING: No OI_FLUX table!")

    return dic

###############################################################################

def show_oi_vs_freq(dic, log=False):
    wl    = dic['WLEN'];
    vis2  = dic['VIS2']['VIS2'];
    vis2e = dic['VIS2']['VIS2ERR'];
    u     = dic['VIS2']['U'];
    v     = dic['VIS2']['V'];
    cp    = dic['T3']['CLOS'];
    cpe   = dic['T3']['CLOSERR'];
    u1    = dic['T3']['U1'];
    v1    = dic['T3']['V1'];
    u2    = dic['T3']['U2'];
    v2    = dic['T3']['V2'];

    plt.figure(figsize = (9, 6))
    G = gridspec.GridSpec(2, 1)

    axes_v2 = plt.subplot(G[0,:])
    
    # Plot all data first
    for i,j in enumerate(u):
        r = np.sqrt(u[i]**2 + v[i]**2);
        freq = r/wl;
        if log:
            axes_v2.semilogy(freq, vis2[i,:],color='lightgray')
            plt.ylim([1e-4,1.1])
        else:
            axes_v2.plot(freq, vis2[i,:],color='lightgray')
            
    # Plot valid data
    for i,j in enumerate(u):
        r = np.sqrt(u[i]**2 + v[i]**2);
        freq = r/wl;
        test = np.logical_and(vis2[i,:] >= 0, vis2e[i,:] / vis2[i,:] < 1)
        if log:
            axes_v2.semilogy(freq[test], vis2[i,test])
            plt.ylim([1e-4,1.1])
        else:
            axes_v2.plot(freq[test], vis2[i,test])
            
    plt.ylim([-0.1,1.1])
    plt.ylabel('V2')
    #plt.xlabel('Spatial Frequency (B/$\lambda$)')
    axes_v2.set_title('Squared visibilities vs frequencies')

    # Plot all data first
    axes_cp = plt.subplot(G[1,:])
    for i,j in enumerate(u1):
        r1 = np.sqrt(u1[i]**2 + v1[i]**2);
        r2 = np.sqrt(u2[i]**2 + v2[i]**2);
        r3 = np.sqrt((u1[i]+u2[i])**2 + (v1[i]+v2[i])**2);
        freq = np.maximum(np.maximum(r1,r2),r3)/wl;
        axes_cp.plot(freq, cp[i,:],color='lightgray')
        
    # Plot valid data only
    for i,j in enumerate(u1):
        r1 = np.sqrt(u1[i]**2 + v1[i]**2);
        r2 = np.sqrt(u2[i]**2 + v2[i]**2);
        r3 = np.sqrt((u1[i]+u2[i])**2 + (v1[i]+v2[i])**2);
        freq = np.maximum(np.maximum(r1,r2),r3)/wl;
        test = np.absolute(cpe[i,:]) < 180/math.pi/3
        axes_cp.plot(freq[test], cp[i,test])
        
    plt.ylim([-200,200])
    axes_cp.set_title('Closure phase vs frequencies')
    plt.ylabel('Closure phase')
    plt.xlabel('Spatial Frequency (B/$\lambda$)')

    plt.show()

###############################################################################

def show_oi_vs_wlen(dic,wlen,datatype="VIS2"):
    data  = dic[datatype];
    datae = dic[datatype+"ERR"];

    for i,j in enumerate(data):
        plt.plot(wlen*1e6, data[i,:])

###############################################################################
# This function shows the selected oifits data (flux, visibility, closure phase etc.)
# as a function of time. It reads data from multiple oifits files in a given
# directory.
# The data to be plotted can be filtered with the filter_oi_list function.
# Example usage:
# filtered_list_of_dicts = filter_oi_list(oifits_dir,dates=["2018-03-14"],bands=['LM'],spectral_resolutions=['MED'],DIT_range=[0,0.2],targets=['l pup'])
# show_oi_vs_time(filtered_list_of_dicts, [3.5, 3.95], key="VIS2", datatype='VIS2') #[3.5, 3.95] [10.2,10.9]
#
def show_oi_vs_time(list_of_dicts, wlenRange,key="VIS2",datatype="VIS2"):
    plot_colors = ['red','blue','green','gold','magenta','cyan','orange','pink','purple','darkgreen']
    # wl = []
    # data = [] #{datatype: []}
    # datae = [] #{datatype+"ERR": []}
    # datat = [] #{"TIME": []}
    datax = []
    datay = []
    datayerr = []
    j=0
    n_rows = -1
    irregular_data_flag = 0
    for dic in list_of_dicts:
        if dic:
            #wl.append(dic['WLEN'])
            #data.append(dic[key][datatype])
            #datae.append(dic[key][datatype+"ERR"])
            #datat.append(dic[key]["TIME"])
            wl = np.array(dic['WLEN'])
            data = np.array(dic[key][datatype])
            datae = np.array(dic[key][datatype+"ERR"])
            datat = np.array(dic[key]["TIME"])
            if data.shape[0] > 6:
                irregular_data_flag = 1
            #print data.shape[0],irregular_data_flag
            wlenRange_idx = np.logical_and(wl > wlenRange[0]/1.0e6,wl < wlenRange[1]/1.0e6)
            if irregular_data_flag == 0:
                datax.append(datat)
                datay.append(np.nanmean(data[:,wlenRange_idx],axis=1))
                datayerr.append(np.nanmean(datae[:,wlenRange_idx],axis=1))
            else:
                datax.append(datat[0:6])
                datay.append(np.nanmean(data[0:6,wlenRange_idx],axis=1))
                datayerr.append(np.nanmean(datae[0:6,wlenRange_idx],axis=1))
    datax = np.array(datax)
    datay = np.array(datay)
    datayerr = np.array(datayerr)
    #print datay.shape[0]
    if datay.shape[0] > 0:
        n_rows = datay.shape[1]
        plt.figure(figsize=(9, 6))
        for i in range(n_rows):
            plt.errorbar(datax[:,i], datay[:,i],yerr=datayerr[:,i],fmt='o',color=plot_colors[i],elinewidth=1.5)
        if "VIS" in datatype:
            plt.ylim([-0.1, 1.1])
        plt.ylabel(datatype)
        plt.xlabel('MJD')
        plt.title(datatype+' vs time')
        plt.show()
    else:
        print "No data to plot."

###############################################################################
# dates = example format: ["2018-03-16"]
# bands = 'LM', 'N'
# spectral_resolutions: 'LOW','MED','HIGH'
# DIT_range: [min,max] (s)
# targets =
def filter_oi_list(input_dir, dates=[],bands=[],spectral_resolutions=[],DIT_range=[],targets=[]):
    oifits_file_list = glob.glob(input_dir + '/*')
    N_files = len(oifits_file_list)
    filtered_list_of_dicts = []
    for file in oifits_file_list:
        dic = open_oi(file)
        if dic:
            date = dic['DATEOBS'][0:10]
            if dates:
                if date not in dates:
                    continue
            if bands:
                if dic['BAND'] not in bands:
                    continue
            if spectral_resolutions:
                if dic['DISP'] not in spectral_resolutions:
                    continue
            if DIT_range:
                if not(dic['DIT'] >= DIT_range[0] and dic['DIT'] <= DIT_range[1]):
                    continue
            target = dic['TARGET']
            if targets:
                targets = [x.lower().replace("_"," ") for x in targets]
                target = target.lower().replace("_"," ")
                if target not in targets:
                    continue
            print target, date, dic['BAND'], dic['DISP'], dic['DIT']
            filtered_list_of_dicts.append(dic)

    return filtered_list_of_dicts

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

    dic = {};
    if os.path.isfile(name_file):
        print("Reading file "+name_file+"...")
        dic = open_oi(name_file)
        print("Plotting data "+name_file+"...")
        show_oi_vs_freq(dic)
        print("Plotting data "+name_file+"...")
        plt.figure()
        show_oi_vs_wlen(dic['FLUX'] ,dic['WLEN'],datatype='FLUX')
        print("Plotting data "+name_file+"...")
        plt.figure()
        wlen = dic['WLEN']
        print(wlen)
    elif os.path.isdir(name_file):
        name_dir = name_file
        filtered_list_of_dicts = filter_oi_list(name_dir)
        show_oi_vs_time(filtered_list_of_dicts ,[3.5,3.95],key="VIS2", datatype='VIS2')


