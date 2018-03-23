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

  Changelog:
  2018-03-23: new functions: oi_data_select_frame, filter_oi_list, open_oi_dir, show_vis2_tf2_vs_time, show_oi_vs_time (jvarga)
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
        print ("Unable to read fits file: "+oi_file)
        return {}

    hdr = hdu[0].header

    wl  = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    dic = {'WLEN':wl}

    target_name = hdu['OI_TARGET'].data['TARGET'][0]
    if not target_name:
        try:
            target_name = hdr['HIERARCH ESO OBS TARG NAME']
        except KeyError:
            print ("Target name not found.")
            target_name = ""
    dic['TARGET'] = target_name
    target_category= hdu['OI_TARGET'].data['CATEGORY'][0] # "CAL" or "SCI"
    if not target_category:
        print ("Target category not found.")
        target_category = "CAL"
    dic['CATEGORY'] = target_category
    try:
        dateobs = hdr['DATE-OBS']
    except KeyError:
        dateobs = ""
    dic['DATEOBS'] = dateobs
    try:
        det_name = hdr['HIERARCH ESO DET CHIP NAME']
    except KeyError:
        print ("Detector name not found.")
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
        print ("Dispersion name not found.")
        dispersion_name = ""
    dic['DISP'] = dispersion_name
    try:
        DIT = hdr["HIERARCH ESO DET SEQ1 DIT"]  # (s)
    except KeyError:
        DIT = np.nan
        print ("DIT not found")
    dic['DIT'] = DIT
    try:
        BCD1 = hdr["HIERARCH ESO INS BCD1 NAME"]
        BCD2 = hdr["HIERARCH ESO INS BCD2 NAME"]
    except KeyError:
        BCD1 = ""
        BCD2 = ""
        print ("BCD NAME not found")
    dic['BCD1NAME'] = BCD1
    dic['BCD2NAME'] = BCD2
    try:
        dic['TEL_NAME'] = hdu['OI_ARRAY'].data["TEL_NAME"]
        dic['STA_NAME'] = hdu['OI_ARRAY'].data["STA_NAME"]
        dic['STA_INDEX'] = hdu['OI_ARRAY'].data["STA_INDEX"]
    except KeyError:
        dic['TEL_NAME'] = {}
        dic['STA_NAME'] = {}
        dic['STA_INDEX'] = {}
        print ("Key in table OI_ARRAY not found")
    try:
        dic['VIS'] = {}
        dic['VIS']['CFLUX']    = hdu['OI_VIS'].data['VISAMP']
        dic['VIS']['CFLUXERR'] = hdu['OI_VIS'].data['VISAMPERR']
        dic['VIS']['DPHI']     = hdu['OI_VIS'].data['VISPHI']
        dic['VIS']['DPHIERR']  = hdu['OI_VIS'].data['VISPHIERR']
        dic['VIS']['U']        = hdu['OI_VIS'].data['UCOORD']
        dic['VIS']['V']        = hdu['OI_VIS'].data['VCOORD']
        dic['VIS']['TIME']     = hdu['OI_VIS'].data['MJD']
        dic['VIS']['STA_INDEX'] = hdu['OI_VIS'].data['STA_INDEX']
    except:
        print("WARNING: No OI_VIS table!")

    try:
        dic['VIS2'] = {}
        dic['VIS2']['VIS2']    = hdu['OI_VIS2'].data['VIS2DATA']
        dic['VIS2']['VIS2ERR'] = hdu['OI_VIS2'].data['VIS2ERR']
        dic['VIS2']['U']       = hdu['OI_VIS2'].data['UCOORD']
        dic['VIS2']['V']       = hdu['OI_VIS2'].data['VCOORD']
        dic['VIS2']['TIME']    = hdu['OI_VIS2'].data['MJD']
        dic['VIS2']['STA_INDEX'] = hdu['OI_VIS2'].data['STA_INDEX']
    except:
        print("WARNING: No OI_VIS2 table!")

    try:
        dic['TF2'] = {}
        dic['TF2']['TF2']    = hdu['OI_TF2'].data['TF2']
        dic['TF2']['TF2ERR'] = hdu['OI_TF2'].data['TF2ERR']
        #dic['TF2']['U']       = hdu['OI_TF2'].data['UCOORD']
        #dic['TF2']['V']       = hdu['OI_TF2'].data['VCOORD']
        dic['TF2']['TIME']     = hdu['OI_TF2'].data['MJD']
        dic['TF2']['STA_INDEX'] =hdu['OI_TF2'].data['STA_INDEX']
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
        dic['T3']['TIME']     = hdu['OI_T3'].data['MJD']
        dic['T3']['STA_INDEX'] = hdu['OI_T3'].data['STA_INDEX']
    except:
        print("WARNING: No OI_T3 table!")

    try:
        dic['FLUX'] = {}
        dic['FLUX']['FLUX']    = hdu['OI_FLUX'].data['FLUXDATA']
        dic['FLUX']['FLUXERR'] = hdu['OI_FLUX'].data['FLUXERR']
        dic['FLUX']['TIME']    = hdu['OI_FLUX'].data['MJD']
        dic['FLUX']['STA_INDEX'] = hdu['OI_FLUX'].data['STA_INDEX']
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
        try:
            #wl.append(dic['WLEN'])
            #data.append(dic[key][datatype])
            #datae.append(dic[key][datatype+"ERR"])
            #datat.append(dic[key]["TIME"])
            wl    = np.array(dic['WLEN'])
            
            data  = np.array(dic[key][datatype])
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
        except:
            pass
    datax = np.array(datax)
    datay = np.array(datay)
    datayerr = np.array(datayerr)
    #print datay.shape[0]
    if datay.shape[0] > 0:
        n_rows = datay.shape[1]
        fig, axs = plt.subplots(n_rows/2, 2, figsize=(15, 10))
        axs = axs.ravel()
        for i in range(n_rows):
            axs[i].errorbar(datax[:,i], datay[:,i],yerr=datayerr[:,i],fmt='o',color=plot_colors[i],elinewidth=1.5)
            if "VIS" in datatype:
                axs[i].set_ylim([-0.1, 1.1])
            axs[i].set_ylabel(datatype)
            axs[i].set_xlabel('MJD')
        plt.suptitle(datatype+' vs time')
        plt.tight_layout()
        plt.show()
    else:
        print "No data to plot."

###############################################################################
def show_vis2_tf2_vs_time(list_of_dicts, wlenRange):
    plot_colors = ['red', 'blue', 'green', 'gold', 'magenta', 'cyan', 'orange', 'pink', 'purple', 'darkgreen']
    target_names_cal = []
    V2_MJD_arr_cal = []
    V2_arr_cal = []
    V2err_arr_cal = []
    V2_sta_index_cal = []

    target_names_CP_cal = []
    CP_MJD_arr_cal = []
    CP_arr_cal = []
    CPerr_arr_cal = []
    CP_sta_index_cal = []

    target_names_TF2 = []
    TF2_MJD_arr = []
    TF2_arr = []
    TF2err_arr = []
    TF2_sta_index = []

    target_names = []
    V2_MJD_arr = []
    V2_arr = []
    V2err_arr = []
    V2_sta_index = []

    target_names_CP = []
    CP_MJD_arr = []
    CP_arr = []
    CPerr_arr = []
    CP_sta_index = []

    for dic in list_of_dicts:
        wl = np.array(dic['WLEN'])
        wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
        category = dic['CATEGORY'].lower()
        if 'cal' in category:
            try:
                datay = np.array(dic['VIS2']['VIS2'])
                datayerr = np.array(dic['VIS2']['VIS2ERR'])
                datax = np.array(dic['VIS2']["TIME"])
                n_rows = datay.shape[0]
                # print datay.shape
                for i in range(n_rows):
                    target_names_cal.append(dic['TARGET'])
                    sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                    V2_sta_index_cal.append(sta_index)
                    # print np.nanmean(datay[i, wlenRange_idx])
                    V2_MJD_arr_cal.append(datax[i])
                    V2_arr_cal.append(np.nanmean(datay[i, wlenRange_idx]))
                    V2err_arr_cal.append(np.nanmean(datayerr[i, wlenRange_idx]))
            except:
                print (dic['TARGET'],dic['DATEOBS'], "No CAL VIS2 data found.")
            try:
                datay = np.array(dic['T3']['CLOS'])
                datayerr = np.array(dic['T3']['CLOSERR'])
                datax = np.array(dic['T3']["TIME"])
                n_rows = datay.shape[0]
                for i in range(n_rows):
                    target_names_CP_cal.append(dic['TARGET'])
                    sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                    CP_sta_index_cal.append(sta_index)
                    CP_MJD_arr_cal.append(datax[i])
                    CP_arr_cal.append(np.nanmean(datay[i, wlenRange_idx]))
                    CPerr_arr_cal.append(np.nanmean(datayerr[i, wlenRange_idx]))
            except:
                print (dic['TARGET'],dic['DATEOBS'], "No CAL CP data found.")
            try:
                datay = np.array(dic['TF2']['TF2'])
                datayerr = np.array(dic['TF2']['TF2ERR'])
                datax = np.array(dic['TF2']["TIME"])
                n_rows = datay.shape[0]
                # print datay.shape
                for i in range(n_rows):
                    target_names_TF2.append(dic['TARGET'])
                    sta_index = np.sort(dic['TF2']['STA_INDEX'][i])
                    TF2_sta_index.append(sta_index)
                    # print np.nanmean(datay[i, wlenRange_idx])
                    TF2_MJD_arr.append(datax[i])
                    TF2_arr.append(np.nanmean(datay[i, wlenRange_idx]))
                    TF2err_arr.append(np.nanmean(datayerr[i, wlenRange_idx]))
            except:
                print (dic['TARGET'],dic['DATEOBS'],"No CAL TF2 data found.")
        if 'sci' in category:
            try:
                datay = np.array(dic['VIS2']['VIS2'])
                datayerr = np.array(dic['VIS2']['VIS2ERR'])
                datax = np.array(dic['VIS2']["TIME"])
                n_rows = datay.shape[0]
                # print datay.shape
                for i in range(n_rows):
                    target_names.append(dic['TARGET'])
                    sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                    V2_sta_index.append(sta_index)
                    # print np.nanmean(datay[i, wlenRange_idx])
                    V2_MJD_arr.append(datax[i])
                    V2_arr.append(np.nanmean(datay[i, wlenRange_idx]))
                    V2err_arr.append(np.nanmean(datayerr[i, wlenRange_idx]))
            except:
                print (dic['TARGET'],dic['DATEOBS'],"No SCI VIS2 data found.")
            try:
                datay = np.array(dic['T3']['CLOS'])
                datayerr = np.array(dic['T3']['CLOSERR'])
                datax = np.array(dic['T3']["TIME"])
                n_rows = datay.shape[0]
                for i in range(n_rows):
                    target_names_CP.append(dic['TARGET'])
                    sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                    CP_sta_index.append(sta_index)
                    CP_MJD_arr.append(datax[i])
                    CP_arr.append(np.nanmean(datay[i, wlenRange_idx]))
                    CPerr_arr.append(np.nanmean(datayerr[i, wlenRange_idx]))
            except:
                print (dic['TARGET'],dic['DATEOBS'],"No SCI CP data found.")

    sta_names = dic['STA_NAME']

    target_names_cal = np.array(target_names_cal)
    V2_MJD_arr_cal = np.array(V2_MJD_arr_cal)
    V2_arr_cal = np.array(V2_arr_cal)
    V2err_arr_cal = np.array(V2err_arr_cal)
    V2_sta_index_cal = np.array(V2_sta_index_cal)

    target_names_CP_cal = np.array(target_names_CP_cal)
    CP_MJD_arr_cal = np.array(CP_MJD_arr_cal)
    CP_arr_cal = np.array(CP_arr_cal)
    CPerr_arr_cal = np.array(CPerr_arr_cal)
    CP_sta_index_cal = np.array(CP_sta_index_cal)

    target_names_TF2 = np.array(target_names_TF2)
    TF2_MJD_arr = np.array(TF2_MJD_arr)
    TF2_arr = np.array(TF2_arr)
    TF2err_arr = np.array(TF2err_arr)
    TF2_sta_index = np.array(TF2_sta_index)

    target_names = np.array(target_names)
    V2_MJD_arr = np.array(V2_MJD_arr)
    V2_arr = np.array(V2_arr)
    V2err_arr = np.array(V2err_arr)
    V2_sta_index = np.array(V2_sta_index)

    target_names_CP = np.array(target_names_CP)
    CP_MJD_arr = np.array(CP_MJD_arr)
    CP_arr = np.array(CP_arr)
    CPerr_arr = np.array(CPerr_arr)
    CP_sta_index = np.array(CP_sta_index)

    sta_indices = np.unique(V2_sta_index_cal, axis=0)
    n_max_config = np.nanmax([6, sta_indices.shape[0]])
    # print sta_indices.shape

    if len(V2_MJD_arr_cal) > 0 and len(V2_MJD_arr) > 0:
        MJD_range = [np.nanmin([np.nanmin(V2_MJD_arr_cal),np.nanmin(V2_MJD_arr)]),np.nanmax([np.nanmax(V2_MJD_arr_cal),np.nanmax(V2_MJD_arr)])]
    elif len(V2_MJD_arr) > 0:
        MJD_range = [np.nanmin(V2_MJD_arr),np.nanmax(V2_MJD_arr)]
    elif len(V2_MJD_arr_cal) > 0:
        MJD_range = [np.nanmin(V2_MJD_arr_cal), np.nanmax(V2_MJD_arr_cal)]
    else:
        MJD_range = [0.0,1.0]
    text_width_MJD = 0.008
    fig1, axs1 = plt.subplots(3, 2, figsize=(15, 13),sharex=True, sharey=True)
    axs1 = axs1.ravel()
    text_y = 1.15
    for i in range(n_max_config):
        #print i
        if len(V2_sta_index_cal) > 0:
            idxst = np.all(V2_sta_index_cal == sta_indices[i],axis=1)
            if len(V2_arr_cal[idxst]) > 0:
                axs1[i].errorbar(V2_MJD_arr_cal[idxst], V2_arr_cal[idxst], yerr=V2err_arr_cal[idxst], fmt='o', color='gray',
                             elinewidth=1.5,label='V2 cal')
                if i in range(2):
                    text_tag_flag = 1
                    prev_text_MJD = 0.0
                    prev_target_name = ""
                    for j in range(np.sum(idxst)):
                        if V2_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                            text_tag_flag = 1
                        if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                            axs1[i].text(V2_MJD_arr_cal[idxst][j], text_y, target_names_cal[idxst][j].replace('_',' '), rotation=90,
                                     va='bottom')
                            text_tag_flag = 0
                            prev_text_MJD = V2_MJD_arr_cal[idxst][j]
                            prev_target_name = target_names_cal[idxst][j]
        if len(TF2_sta_index) > 0:
            idxst = np.all(TF2_sta_index == sta_indices[i],axis=1)
            if len(TF2_arr[idxst]) > 0:
                axs1[i].errorbar(TF2_MJD_arr[idxst], TF2_arr[idxst], yerr=TF2err_arr[idxst], fmt='o', color='blue',
                                 elinewidth=1.5, label='TF2')

        if len(V2_sta_index) > 0:
            idxst = np.all(V2_sta_index == sta_indices[i], axis=1)
            if len(V2_arr[idxst]) > 0:
                axs1[i].errorbar(V2_MJD_arr[idxst], V2_arr[idxst], yerr=V2err_arr[idxst], fmt='o', color='red',
                                elinewidth=1.5,label='V2 sci')
                if i in range(2):
                    text_tag_flag = 1
                    prev_text_MJD = 0.0
                    prev_target_name = ""
                    for j in range(np.sum(idxst)):
                        if V2_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                            text_tag_flag = 1
                        if text_tag_flag == 1 or (prev_target_name != target_names[idxst][j]):
                            axs1[i].text(V2_MJD_arr[idxst][j], text_y, target_names[idxst][j].replace('_', ' '),
                                         rotation=90, va='bottom', color='darkred')
                            text_tag_flag = 0
                            prev_text_MJD = V2_MJD_arr[idxst][j]
                            prev_target_name = target_names[idxst][j]

        axlabel = sta_names[sta_indices[i,0] == dic['STA_INDEX']][0] + ' - ' + \
                  sta_names[sta_indices[i,1] == dic['STA_INDEX']][0]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        axs1[i].text(0.05, 0.95, axlabel,horizontalalignment='left',verticalalignment='top',
                     transform=axs1[i].transAxes, bbox=props)
        leg = axs1[i].legend(loc='upper right')
        leg.get_frame().set_alpha(0.5)
        axs1[i].set_ylim([-0.1, 1.1])
        axs1[i].set_ylabel('$V^2$')
        axs1[i].set_xlabel('$\mathrm{MJD}$')
    plt.suptitle('$V^2\mathrm{\ vs.\ time}$')
    fig1.subplots_adjust(hspace=0,wspace=0)
    for i in range(4):
        plt.setp(axs1[i].get_xticklabels() , visible=False)
        x_axis = axs1[i].axes.get_xaxis()
        x_axis.get_label().set_visible(False)
    for i in range(1,6,2):
        plt.setp(axs1[i].get_yticklabels() , visible=False)
        y_axis = axs1[i].axes.get_yaxis()
        y_axis.get_label().set_visible(False)
    #plt.tight_layout()

    CP_sta_indices = np.unique(CP_sta_index_cal, axis=0)
    #print CP_sta_indices
    n_max_config = np.nanmax([4,CP_sta_indices.shape[0]])

    fig2, axs = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
    axs = axs.ravel()
    text_y = 60
    for i in range(n_max_config):
        if len(CP_sta_index_cal) > 0:
            idxst = np.all(CP_sta_index_cal == CP_sta_indices[i],axis=1)
            if len(CP_arr_cal[idxst]) > 0:
                axs[i + 0].errorbar(CP_MJD_arr_cal[idxst], CP_arr_cal[idxst], yerr=CPerr_arr_cal[idxst], fmt='o',
                                    color='gray',elinewidth=1.5, label='CP cal')
                if i in range(2):
                    text_tag_flag = 1
                    prev_text_MJD = 0.0
                    prev_target_name = ""
                    for j in range(np.sum(idxst)):
                        if CP_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                            text_tag_flag = 1
                        if text_tag_flag == 1 or (prev_target_name != target_names_CP_cal[idxst][j]):
                            axs[i + 0].text(CP_MJD_arr_cal[idxst][j], text_y, target_names_CP_cal[idxst][j].replace('_',' '), rotation=90,
                                     va='bottom')
                            text_tag_flag = 0
                            prev_text_MJD = CP_MJD_arr_cal[idxst][j]
                            prev_target_name = target_names_CP_cal[idxst][j]
        if len(CP_sta_index) > 0:
            idxst = np.all(CP_sta_index == CP_sta_indices[i],axis=1)
            if len(CP_arr[idxst]) > 0:
                axs[i + 0].errorbar(CP_MJD_arr[idxst], CP_arr[idxst], yerr=CPerr_arr[idxst], fmt='o',
                                    color='darkred',elinewidth=1.5, label='CP sci')
                if i in range(2):
                    text_tag_flag = 1
                    prev_text_MJD = 0.0
                    prev_target_name = ""
                    for j in range(np.sum(idxst)):
                        if CP_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                            text_tag_flag = 1
                        if text_tag_flag == 1 or (prev_target_name != target_names_CP[idxst][j]):
                            axs[i + 0].text(CP_MJD_arr[idxst][j], text_y, target_names_CP[idxst][j].replace('_', ' '),
                                             rotation=90,
                                             va='bottom')
                            text_tag_flag = 0
                            prev_text_MJD = CP_MJD_arr[idxst][j]
                            prev_target_name = target_names_CP[idxst][j]
        axlabel = sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                  sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                  sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        axs[i + 0].text(0.05, 0.95, axlabel,horizontalalignment='left',verticalalignment='top',
                     transform=axs[i + 0].transAxes, bbox=props)
        leg = axs[i + 0].legend(loc='upper right')
        leg.get_frame().set_alpha(0.5)
        axs[i + 0].set_ylabel('$CP\,\left(^\circ\\right)$')
        axs[i + 0].set_xlabel('$\mathrm{MJD}$')
    plt.suptitle('$CP\mathrm{\ vs.\ time}$')
    fig2.subplots_adjust(hspace=0, wspace=0)
    for i in range(2):
        plt.setp(axs[i+0].get_xticklabels() , visible=False)
        x_axis = axs[i+0].axes.get_xaxis()
        x_axis.get_label().set_visible(False)
    for i in range(1,4,2):
        plt.setp(axs[i+0].get_yticklabels() , visible=False)
        y_axis = axs[i+0].axes.get_yaxis()
        y_axis.get_label().set_visible(False)
    #plt.tight_layout()
    plt.show()

###############################################################################
def open_oi_dir(input_dir):
    oifits_file_list = glob.glob(input_dir + '/*fits*')

    N_files = len(oifits_file_list)
    list_of_dicts = []
    for file in oifits_file_list:
        dic = open_oi(file)
        if dic:
            print (dic['TARGET'], dic['DATEOBS'], dic['BAND'], dic['DISP'], dic['DIT'],dic['CATEGORY'])
            list_of_dicts.append(dic)

    return list_of_dicts

###############################################################################
# dates = example format: ["2018-03-16"]
# bands = 'LM', 'N'
# spectral_resolutions: 'LOW','MED','HIGH'
# DIT_range: [min,max] (s)
# targets = []
def filter_oi_list(list_of_dicts, dates=[],bands=[],spectral_resolutions=[],DIT_range=[],targets=[]):
    filtered_list_of_dicts = []
    for dic in list_of_dicts:
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
            print ("Selected: ", target, date, dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'])
            filtered_list_of_dicts.append(dic)

    return filtered_list_of_dicts

###############################################################################
# name_dir = r"D:\jvarga\Dokumentumok\MATISSE\data\OIFITS\2018-03-14"
# list_of_dicts = open_oi_dir(name_dir)
# filtered_list_of_dicts = filter_oi_list(list_of_dicts,dates=["2018-03-14"],bands=['LM'],spectral_resolutions=['MED'],DIT_range=[0.18,0.21],targets=[])
# #show_oi_vs_time(filtered_list_of_dicts, [3.5, 3.95], key="VIS2", datatype='VIS2') #[3.5, 3.95] [10.2,10.9]
# print "Selected",len(filtered_list_of_dicts),"objects"
# show_vis2_tf2_vs_time(filtered_list_of_dicts,wlenRange=[3.5,3.95 ]) # wlenRange=[3.5, 3.95]);
# raise SystemExit

class oi_data_select_frame(wx.Frame):
    # Data type selector GUI - work in progress ...
    def __init__(self,  *args, **kw):
        super(oi_data_select_frame, self).__init__(*args, **kw)
        self.InitUI()

    def InitUI(self):
        self.SetTitle('Select OIDATA')
        self.Centre()
        panel = wx.Panel(self, wx.ID_ANY)
        vbox_ultimate = wx.BoxSizer(wx.VERTICAL)
        vbox = wx.BoxSizer(wx.HORIZONTAL)
        vbox2 = wx.BoxSizer(wx.HORIZONTAL)

        sbox = wx.StaticBox(panel, -1, '')
        text = wx.StaticText(sbox, wx.ID_ANY, "")
        sboxSizer = wx.StaticBoxSizer(sbox, wx.VERTICAL)

        hbox = wx.BoxSizer(wx.VERTICAL)
        b1=wx.Button(panel, 0, label='VISAMP, VISPHI')
        b2=wx.Button(panel, 1, label='VIS2, T3PHI')
        b3=wx.Button(panel, 2, label='T3AMP')
        b4=wx.Button(panel, 3, label='TF2')
        b5=wx.Button(panel, 4, label='FLUX')

        hbox.Add(b1, 0, wx.ALL | wx.LEFT, 3)
        hbox.Add(b2, 0, wx.ALL | wx.LEFT, 3)
        hbox.Add(b3, 0, wx.ALL | wx.LEFT, 3)
        hbox.Add(b4, 0, wx.ALL | wx.LEFT, 3)
        hbox.Add(b5, 0, wx.ALL | wx.LEFT, 3)
        self.Bind(wx.EVT_BUTTON, self.OnButtonClicked)
        sboxSizer.Add(hbox, 0, wx.ALL | wx.LEFT, 3)

        self.statusbar = self.CreateStatusBar(1)
        self.statusbar.SetStatusText('')

        nm = wx.StaticBox(panel, -1, 'Band')
        nmSizer = wx.StaticBoxSizer(nm, wx.VERTICAL)
        nmbox = wx.BoxSizer(wx.VERTICAL)
        nmbox.Add(wx.CheckBox(panel,0, label='LM'), 0, wx.ALL | wx.LEFT, 3)
        nmbox.Add(wx.CheckBox(panel, 1, label='N'), 0, wx.ALL | wx.LEFT, 3)
        nmSizer.Add(nmbox, 0, wx.ALL | wx.UP, 3)

        nm2 = wx.StaticBox(panel, -1, 'Spectral resolution')
        nmSizer2 = wx.StaticBoxSizer(nm2, wx.VERTICAL)
        nmbox2= wx.BoxSizer(wx.VERTICAL)
        nmbox2.Add(wx.CheckBox(panel, 2, label='LOW'), 0, wx.ALL | wx.LEFT, 3)
        nmbox2.Add(wx.CheckBox(panel, 3, label='MED'), 0, wx.ALL | wx.LEFT, 3)
        nmbox2.Add(wx.CheckBox(panel, 4, label='HIGH'), 0, wx.ALL | wx.LEFT, 3)
        nmbox2.Add(wx.CheckBox(panel, 5, label='CRAZY'), 0, wx.ALL | wx.LEFT, 3)
        nmSizer2.Add(nmbox2, 0, wx.ALL | wx.UP, 3)
        self.Bind(wx.EVT_CHECKBOX, self.onChecked)

        #dateobs
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        l1 = wx.StaticText(panel, -1, "Text Field")
        hbox1.Add(l1, 1, wx.EXPAND | wx.ALIGN_LEFT | wx.ALL, 5)
        self.t1 = wx.TextCtrl(panel)
        hbox1.Add(self.t1, 1, wx.EXPAND | wx.ALIGN_LEFT | wx.ALL, 5)
        self.t1.Bind(wx.EVT_TEXT, self.OnKeyTyped)

        #DIT

        #Target

        vbox.Add(sboxSizer, 0, wx.ALL | wx.LEFT, 3)
        vbox.Add(nmSizer, 0, wx.ALL | wx.LEFT, 3)
        vbox.Add(nmSizer2, 0, wx.ALL | wx.LEFT, 3)
        vbox2.Add(hbox1, wx.ALL | wx.LEFT, 3)
        vbox_ultimate.Add(vbox,0, wx.ALL | wx.UP, 3)
        vbox_ultimate.Add(vbox2, 0, wx.ALL | wx.CENTER, 3)
        panel.SetSizer(vbox)
        self.Centre()

        panel.Fit()
        self.Show()

    def OnKeyTyped(self, event):
        print (event.GetString())

    def onChecked(self, e):
        cb = e.GetEventObject()
        print (cb.GetLabel(), ' is clicked', cb.GetValue())

    def OnButtonClicked(self, e):
        eID = e.GetId()
        if eID == 0:
            self.statusbar.SetStatusText('Show VISAMP & VISPHI.')
            print ('Show VISAMP.')
        elif eID == 1:
            self.statusbar.SetStatusText('Show VIS2 & T3PHI.')
            print (r'Show VIS2 & T3PHI.')
        elif eID == 2:
            self.statusbar.SetStatusText('Show T3AMP.')
            print ('Show T3AMP.')
        elif eID == 3:
            self.statusbar.SetStatusText('Show TF2.')
            print ('Show TF2.')
        elif eID == 4:
            self.statusbar.SetStatusText('Show FLUX.')
            print ('Show FLUX.')
        #e.Skip()

#ex = wx.App()
#oi_data_select_frame(None)
#ex.MainLoop()


if __name__ == '__main__':
    listArg = sys.argv
    name_file = []
    typePlot = "VIS2"
    for elt in listArg:
        if ('--help' in elt):
            print( "Usage: mat_show_rawdata.py [--dir=start directory]")
            sys.exit(0)
        elif ('--typePlot' in elt):
            typePlot = elt.split("=")[1]
        #elif len(listArg) == 2:
       #     name_file = sys.argv[1]
       #     print(name_file)

    if typePlot=="VIS2":
        table = "VIS2"
    elif typePlot=="CLOS":
        table="T3"
    elif typePlot=="TF2":
        typePlot="TF2"
        table="TF2"
    elif typePlot=="FLUX":
        table="FLUX"
    elif typePlot=="DPHI":
        table="VIS"
        

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
        list_of_dicts = open_oi_dir(name_dir)
        filtered_list_of_dicts = filter_oi_list(list_of_dicts)
        show_vis2_tf2_vs_time(filtered_list_of_dicts,[3.5,3.95])
        #show_oi_vs_time(filtered_list_of_dicts ,[3.5,3.95],key=table, datatype=typePlot)
        
