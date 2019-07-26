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
  2018-03-23: new functions: oi_data_select_frame, filter_oi_list, open_oi_dir,
              show_vis2_tf2_vs_time, show_oi_vs_time (jvarga)
  2018-03-26: new GUI interface ready: oi_data_select_frame (jvarga)
  2018-04-04: updated GUI and extended functionality: input file/folder
              textbox, filter for target name, more bands (JHK) available
              (for e.g. AMBER data), plot with or without errorbars, plot V or
              V2 (jvarga)
"""

import sys
import wx
import math
import numpy as np
from   matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from   mat_fileDialog import mat_FileDialog
from   mat_fileDialog import identifyFile
from   astropy.io import fits as fits
import os
import glob
import libRobust as robust
#from astroquery.simbad import Simbad
from astropy import coordinates
from os.path import expanduser
from matplotlib.ticker import *

home = expanduser("~")

###############################################################################

def resolve_target(dic):
    try:
        ra  = str(dic['HDR']["RA"])
        dec = str(dic['HDR']["DEC"])
        c = coordinates.SkyCoord(ra, dec, unit=('deg','deg'), frame='icrs')
        result_table = Simbad.query_region(c)
        resolvedTarg = result_table['MAIN_ID'][0]
    except:
        resolvedTarg = dic['TARGET']
    return resolvedTarg;

###############################################################################

def open_oi(oi_file):
    try:
        hdu = fits.open(oi_file)
    except IOError:
        print("Unable to read fits file: " + oi_file)
        return {}

    hdr = hdu[0].header

    wl = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    dic = {'WLEN': wl}

    dic['HDR']  = hdr
    dic['file'] = oi_file;

    try:
        dic['SEEING'] = (hdr['HIERARCH ESO ISS AMBI FWHM START'] + hdr['HIERARCH ESO ISS AMBI FWHM END'])/2.
    except:
        dic['SEEING'] = 0;

    try:
        dic['TAU0']   = (hdr['HIERARCH ESO ISS AMBI TAU0 START'] + hdr['HIERARCH ESO ISS AMBI TAU0 END'])/2.
    except:
        dic['TAU0'] = 0;

    target_name = hdu['OI_TARGET'].data['TARGET'][0]
    if not target_name:
        try:
            target_name = hdr['HIERARCH ESO OBS TARG NAME']
        except KeyError:
            print ("Target name not found.")
            target_name = ""

    dic['TARGET'] = target_name

    # Fix eventual bad target identification
    #dic['TARGET'] = resolve_target(dic)

    try:
        target_category = hdu['OI_TARGET'].data['CATEGORY'][0]  # "CAL" or "SCI"
    except KeyError:
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
        if (det_name == 'AQUARIUS'):
          dispersion_name = hdr['HIERARCH ESO INS DIN NAME']
        else :
          dispersion_name = hdr['HIERARCH ESO INS DIL NAME']
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
        dic['TEL_NAME']  = hdu['OI_ARRAY'].data["TEL_NAME"]
        dic['STA_NAME']  = hdu['OI_ARRAY'].data["STA_NAME"]
        dic['STA_INDEX'] = hdu['OI_ARRAY'].data["STA_INDEX"]
    except KeyError:
        dic['TEL_NAME']  = {}
        dic['STA_NAME']  = {}
        dic['STA_INDEX'] = {}
        print ("Key in table OI_ARRAY not found")
    try:
        dic['VIS'] = {}
        dic['VIS']['VISAMP']    = hdu['OI_VIS'].data['VISAMP']
        dic['VIS']['VISAMPERR'] = hdu['OI_VIS'].data['VISAMPERR']
        dic['VIS']['DPHI']      = hdu['OI_VIS'].data['VISPHI']
        dic['VIS']['DPHIERR']   = hdu['OI_VIS'].data['VISPHIERR']
        dic['VIS']['FLAG']   = hdu['OI_VIS'].data['FLAG']
        try:
            dic['VIS']['CFLUX']    = hdu['OI_VIS'].data['CFXAMP']
            dic['VIS']['CFLUXERR'] = hdu['OI_VIS'].data['CFXAMPERR']
        except:
            print("WARNING: No correlated fluxes in this OI_VIS table!")

        dic['VIS']['U']         = hdu['OI_VIS'].data['UCOORD']
        dic['VIS']['V']         = hdu['OI_VIS'].data['VCOORD']
        dic['VIS']['TIME']      = hdu['OI_VIS'].data['MJD']
        if dic['VIS']['TIME'][0] < 50000:
            dic['VIS']['TIME']      = np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])
        dic['VIS']['STA_INDEX'] = hdu['OI_VIS'].data['STA_INDEX']
    except:
        print("WARNING: No OI_VIS table!")

    try:
        dic['VIS2'] = {}
        dic['VIS2']['VIS2']      = hdu['OI_VIS2'].data['VIS2DATA']
        dic['VIS2']['VIS2ERR']   = hdu['OI_VIS2'].data['VIS2ERR']
        dic['VIS2']['U']         = hdu['OI_VIS2'].data['UCOORD']
        dic['VIS2']['V']         = hdu['OI_VIS2'].data['VCOORD']
        dic['VIS2']['TIME']      = hdu['OI_VIS2'].data['MJD']
        dic['VIS2']['FLAG']   = hdu['OI_VIS2'].data['FLAG']
        if dic['VIS2']['TIME'][0] < 50000:
            print("WARNING: incoherent MJD, picking it up from header")
            print(np.shape(hdu['OI_VIS2'].data['MJD']))
            dic['VIS2']['TIME'] = np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])
            print(np.shape(np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])))
        dic['VIS2']['STA_INDEX'] = hdu['OI_VIS2'].data['STA_INDEX']
    except:
        print("WARNING: No OI_VIS2 table!")

    try:
        dic['TF2'] = {}
        dic['TF2']['TF2']       = hdu['TF2'].data['TF2']
        dic['TF2']['TF2ERR']    = hdu['TF2'].data['TF2ERR']
        # dic['TF2']['U']       = hdu['OI_TF2'].data['UCOORD']
        # dic['TF2']['V']       = hdu['OI_TF2'].data['VCOORD']
        dic['TF2']['TIME']      = hdu['TF2'].data['MJD']
        if dic['TF2']['TIME'][0] < 50000:
            dic['TF2']['TIME'] = np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])
        dic['TF2']['STA_INDEX'] = hdu['TF2'].data['STA_INDEX']
    except:
        print("WARNING: No OI_TF2 table!")

    try:
        dic['T3'] = {}
        dic['T3']['T3AMP']     = hdu['OI_T3'].data['T3AMP']
        dic['T3']['T3AMPERR']  = hdu['OI_T3'].data['T3AMPERR']
        dic['T3']['CLOS']      = hdu['OI_T3'].data['T3PHI']
        dic['T3']['CLOSERR']   = hdu['OI_T3'].data['T3PHIERR']
        dic['T3']['U1']        = hdu['OI_T3'].data['U1COORD']
        dic['T3']['V1']        = hdu['OI_T3'].data['V1COORD']
        dic['T3']['U2']        = hdu['OI_T3'].data['U2COORD']
        dic['T3']['V2']        = hdu['OI_T3'].data['V2COORD']
        dic['T3']['TIME']      = hdu['OI_T3'].data['MJD']
        dic['T3']['FLAG']   = hdu['OI_T3'].data['FLAG']
        if dic['T3']['TIME'][0] < 50000:
            dic['T3']['TIME']      = np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])
        dic['T3']['STA_INDEX'] = hdu['OI_T3'].data['STA_INDEX']
    except:
        print("WARNING: No OI_T3 table!")

    try:
        dic['FLUX'] = {}
        dic['FLUX']['FLUX']      = hdu['OI_FLUX'].data['FLUXDATA']
        dic['FLUX']['FLUXERR']   = hdu['OI_FLUX'].data['FLUXERR']
        dic['FLUX']['TIME']      = hdu['OI_FLUX'].data['MJD']
        dic['FLUX']['FLAG']      = hdu['OI_FLUX'].data['FLAG']
        if dic['FLUX']['TIME'][0] < 50000:
            dic['FLUX']['TIME']      = np.full(len(hdu['OI_VIS'].data['MJD']), hdr['MJD-OBS'])
        dic['FLUX']['STA_INDEX'] = hdu['OI_FLUX'].data['STA_INDEX']
    except:
        print("WARNING: No OI_FLUX table!")

    return dic

###############################################################################

def show_oi_vs_freq(dic, log=False,showvis=False):
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

    plt.figure(figsize=(9, 6))
    G = gridspec.GridSpec(2, 1)

    axes_v2 = plt.subplot(G[0, :])

    # Visibility
    # Plot all data first
    for i, j in enumerate(u):
        r = np.sqrt(u[i] ** 2 + v[i] ** 2);
        freq = r / wl;
        if log:
            if showvis == True:
                axes_v2.semilogy(freq, vis2[i, :], color='lightgray')
            else:
                axes_v2.semilogy(freq, np.sqrt(vis2[i, :]), color='lightgray')
            plt.ylim([1e-4, 1.1])
        else:
            if showvis == True:
                axes_v2.plot(freq, np.sqrt(vis2[i, :]), color='lightgray')
            else:
                axes_v2.plot(freq, vis2[i, :], color='lightgray')

    # Plot valid data
    for i, j in enumerate(u):
        r = np.sqrt(u[i] ** 2 + v[i] ** 2);
        freq = r / wl;
        test = np.logical_and(vis2[i, :] >= 0, vis2e[i, :] / vis2[i, :] < 1)
        if log:
            if showvis == True:
                axes_v2.semilogy(freq[test], np.sqrt(vis2[i, test]))
            else:
                axes_v2.semilogy(freq[test], vis2[i, test])
            plt.ylim([1e-4, 1.1])
        else:
            if showvis == True:
                axes_v2.plot(freq[test], np.sqrt(vis2[i, test]))
            else:
                axes_v2.plot(freq[test], vis2[i, test])

    plt.ylim([-0.1, 1.1])
    if showvis == True:
        plt.ylabel('V')
        axes_v2.set_title('Visibilities vs frequencies')
    else:
        plt.ylabel('V2')
        axes_v2.set_title('Squared visibilities vs frequencies')
    # plt.xlabel('Spatial Frequency (B/$\lambda$)')

    # Closure phase
    # Plot all data first
    axes_cp = plt.subplot(G[1, :])
    for i, j in enumerate(u1):
        r1 = np.sqrt(u1[i] ** 2 + v1[i] ** 2);
        r2 = np.sqrt(u2[i] ** 2 + v2[i] ** 2);
        r3 = np.sqrt((u1[i] + u2[i]) ** 2 + (v1[i] + v2[i]) ** 2);
        freq = np.maximum(np.maximum(r1, r2), r3) / wl;
        axes_cp.plot(freq, cp[i, :], color='lightgray')

    # Plot valid data only
    for i, j in enumerate(u1):
        r1 = np.sqrt(u1[i] ** 2 + v1[i] ** 2);
        r2 = np.sqrt(u2[i] ** 2 + v2[i] ** 2);
        r3 = np.sqrt((u1[i] + u2[i]) ** 2 + (v1[i] + v2[i]) ** 2);
        freq = np.maximum(np.maximum(r1, r2), r3) / wl;
        test = np.absolute(cpe[i, :]) < 180 / math.pi / 3
        axes_cp.plot(freq[test], cp[i, test])

    plt.ylim([-200, 200])
    axes_cp.set_title('Closure phase vs frequencies')
    plt.ylabel('Closure phase')
    plt.xlabel('Spatial Frequency (B/$\lambda$)')

    plt.show()

###############################################################################

def show_oi_vs_wlen(dic, key='VIS2', datatype="VIS2", showvis=False,
                    plot_errorbars=True, correct_polynom=False,
                    timeVertOffset=0, stdevRange=False, normRange=False,
                    stdevTime=False,subplotList=None,colorList=None,useFlag=False):
   # plot_colors = ['red', 'blue', 'green', 'gold', 'magenta', 'cyan', 'orange', 'pink', 'purple', 'darkgreen']

    sta_index_cal = []

    # Get data from the input dictionary
    wl     = dic['WLEN'];
    nbWlen = len(wl)
    data   = dic[key][datatype];
    flag   = dic[key]['FLAG']

    flag   = dic[key]['FLAG']

    idxGood=[]
    for i in range(np.shape(flag)[0]):
        if useFlag:
            idxGood.append(np.where(np.logical_not(flag[i]))[0])
        else:
            idxGood.append(np.arange(np.size(flag[i])))

    datae  = dic[key][datatype + "ERR"];
    time   = dic[key]["TIME"];

    if correct_polynom:
        #print("Correct a polynom")
        if normRange:
            l = normRange[0]
            h = normRange[1]
        else:
            l =  3
            h = -3
        for i,j in enumerate(data):
            fit = np.polyfit(wl[l:h],data[i,l:h],correct_polynom)
            if key == 'VIS2' or key == 'TF2' or key == 'FLUX' or datatype == 'CFLUX':
                data[i,:] = data[i,:] / np.polyval(fit,wl)
            else:
                data[i,:] = data[i,:] - np.polyval(fit,wl)

    if key == 'FLUX':
        sta_index = dic[key]['STA_INDEX']
        sta_index_cal.append([sta_index])
    else:
        sta_index = np.sort(dic[key]['STA_INDEX'])
        sta_index_cal.append(sta_index)

    # Get the unique station indices from the data
    #print("Get the unique station indices from the data")
    sta_indices = np.unique(sta_index_cal, axis=0)
    if key == 'VIS' or key == 'VIS2' or key == 'TF2':
        n_max_config = np.nanmax([6, sta_indices.shape[0]])
        n_plot_rows = 3
    elif key == 'FLUX':
        n_max_config = 4#np.nanmax([4, sta_indices.shape[0]])
        n_plot_rows = 2
    elif key == 'T3':
        n_max_config = np.nanmax([4, sta_indices.shape[0]])
        n_plot_rows = 2

    #print(n_max_config)


    if not(subplotList):
        fig1, axs1 = plt.subplots(n_plot_rows, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
    else:
        axs1 = subplotList

    #print datae
    #print datae.shape
    #print("Plotting...")
    npl=int(len(data)/n_max_config)
    #print(npl)
    for i in range(int(len(data)/n_max_config)):

        label = datatype
        #

        for j in range(n_max_config):
            idx = int(i*n_max_config+j);

            # Take square root of observable
            if showvis == True:
                print("Take square root of observable")
                data[idx, :]  = np.sqrt(data[idx, :])
                datae[idx, :] = 0.5 * datae[idx, :] / np.sqrt(data[idx, :])
                if key == 'VIS2' or key == 'TF2':
                    if key == 'VIS2':
                        label = 'VIS'
                    elif key == 'TF2':
                        label = 'TF'


            if not(colorList):
                axs1[j].plot(wl[idxGood[idx]] * 1e6, data[idx, idxGood[idx]]+timeVertOffset*j)
            else:
                axs1[j].plot(wl[idxGood[idx]] * 1e6, data[idx, idxGood[idx]]+timeVertOffset*j,color=colorList[j],alpha=1-float(i)/float(npl))
            if not(subplotList):
                axs1[j].set_ylabel(label)
                axs1[j].set_xlabel(r"$\lambda\ (\mu\mathrm{m}$)")
                #print(np.shape(data))
                if stdevRange:
                    if datatype == 'CFLUX':
                        off = 1.05;
                    else:
                        off = 0.5;
                    axs1[j].text(wl[stdevRange[1]] * 1e6, timeVertOffset*j+off, "st. dev. in continuum"+str(round(np.std(data[idx, stdevRange[0]:stdevRange[1]]),3)))
            #print("st. dev.",np.std(data[idx, normrange[0]:normrange[1]]))

            if plot_errorbars == True:
                axs1[j].errorbar(wl * 1e6, data[idx, :],
                    yerr=datae[idx, :],
                    ecolor = 'grey', alpha = 0.25, capsize = 0.5,elinewidth = 1)

        if stdevTime == True:
            dat = np.reshape(data,[int(len(data)/n_max_config),n_max_config,nbWlen])
            std = np.std(dat,2)
            if stdevRange:
                axs1[j].text(wl[stdevRange[1]] * 1e6, timeVertOffset*j+off, "st. dev. vs. time"+str(round(np.mean(np.std(dat[:,i, stdevRange[0]:stdevRange[1]],axis=0),3))))

    if (datatype == 'VIS2' or datatype == 'TF2' or datatype == 'CFLUX') and not(subplotList):
        plt.ylim([-0.1,1.1+timeVertOffset*(n_max_config-1)*1.1])
    elif not(subplotList):
        try:
            mn = np.nanmin(data);
            if mn > 0:
                mn = mn*0.9;
            else:
                mn = mn*1.1
            plt.ylim([mn,np.nanmax(data+timeVertOffset*(n_max_config-1))*1.1])
            #plt.ylim([-30,60])
        except:
            pass
    if showvis == True:
        if datatype == 'VIS2':
            datatype = 'VIS'
        elif datatype == 'TF2':
            datatype = 'TF'
    if not(subplotList):
        fig1.suptitle(datatype+' vs. wavelength')
        plt.show()


###############################################################################
# This function shows the selected oifits data (flux, visibility,
# closure phase etc.) as a function of any keywork in the header (by
# default the seeing). It reads data from multiple oifits files in a
# given directory. The data to be plotted can be filtered with the
# filter_oi_list function.
#
def show_oi_vs_anything(list_of_dicts, wlenRange, key="TF2",
datatype="TF2", xaxis="HIERARCH ESO ISS AMBI FWHM START",
showvis=False,plot_errorbars=True,
useStations=True,xlog=False,ylog=False):
    # check if list is not empty:
    if list_of_dicts:
        target_names_cal = []
        XKey_arr_cal     = []
        arr_cal          = []
        err_arr_cal      = []
        sta_index_cal    = []

        for dic in list_of_dicts:
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            category = dic['CATEGORY'].lower()
            if 'cal' in category:
                try:
                    datay    = np.array(dic[key][datatype])
                    datayerr = np.array(dic[key][datatype+'ERR'])
                    #print("ready?")
                    datax    = dic['HDR'][xaxis]
                    #print(datax)
                    n_rows   = datay.shape[0]
                    # print datay.shape
                    for i in range(n_rows):
                        XKey_arr_cal.append(datax)
                        arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                        err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                        target_names_cal.append(dic['TARGET'])

                        if useStations==True:
                            if key == 'FLUX':
                                sta_index = dic[key]['STA_INDEX'][i]
                                sta_index_cal.append([sta_index])
                                print(np.shape(sta_index))
                            else:
                                sta_index = np.sort(dic[key]['STA_INDEX'][i])
                                sta_index_cal.append(sta_index)
                                # print dic[key]['STA_INDEX'][i]
                                print(np.shape(sta_index))
                        else:
                            sta_index = i%6
                            sta_index_cal.append([sta_index])

                except:
                    print(dic['TARGET'], dic['DATEOBS'], "No CAL data found.")

            if useStations==True:
                sta_names = dic['STA_NAME']
            else:
                sta_names = ""

        target_names_cal = np.array(target_names_cal)
        XKey_arr_cal     = np.array(XKey_arr_cal)
        arr_cal          = np.array(arr_cal)
        err_arr_cal      = np.array(err_arr_cal)
        sta_index_cal    = np.array(sta_index_cal)

        sta_indices = np.unique(sta_index_cal, axis=0)
        print(len(sta_indices))
        if key == 'VIS' or key == 'VIS2' or key == 'TF2':
            n_max_config = np.nanmax([6, sta_indices.shape[0]])
            n_plot_rows = 3
        elif key == 'T3' or key == 'FLUX':
            n_max_config = np.nanmax([4, sta_indices.shape[0]])
            n_plot_rows = 2
        # print sta_indices.shape
        if len(XKey_arr_cal) > 0:
            XKey_range = [np.nanmin(XKey_arr_cal), np.nanmax(XKey_arr_cal)]
        else:
            XKey_range = [0.0, 1.0]
        text_width_XKey = (XKey_range[1] - XKey_range[0]) / 20.0

        fig1, axs1 = plt.subplots(n_plot_rows, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        for i in range(n_max_config):
            # print i
            if datatype == 'DPHI' or datatype == 'CLOS':
                axs1[i + 0].plot(XKey_range, [0.0, 0.0], '-', color='gray', lw=1.5)

            if len(sta_index_cal) > 0:
                idxst = np.all(sta_index_cal == sta_indices[i], axis=1)
                if len(arr_cal[idxst]) > 0:
                    label = datatype +' cal'

                    #print(len(XKey_arr_cal))
                    #print(len(arr_cal))

                    if plot_errorbars == True:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' cal'
                                elif key == 'TF2':
                                    label = 'TF' + ' cal'
                                axs1[i].errorbar(XKey_arr_cal[idxst],
                                                 np.sqrt(arr_cal[idxst]),
                                                 yerr=0.5 * err_arr_cal[idxst] / np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color='blue',
                                                 elinewidth=1.0,
                                                 label=label)
                        else:
                            print(i)
                            axs1[i].errorbar(XKey_arr_cal[idxst],
                                             arr_cal[idxst],
                                             yerr=err_arr_cal[idxst],
                                             fmt='o', color='blue',
                                             elinewidth=1.0,
                                             label=label)
                    else:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(XKey_arr_cal[idxst],
                                                 np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color='blue',
                                                 elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(XKey_arr_cal[idxst],
                                             arr_cal[idxst],
                                             fmt='o', color='blue',
                                             elinewidth=1.0,
                                             label=label)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_XKey = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if XKey_arr_cal[idxst][j] > (prev_text_XKey + text_width_XKey):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                ymin, ymax = axs1[i + 0].get_ylim()
                                axs1[i].text(XKey_arr_cal[idxst][j], ymax * 1.05,
                                             target_names_cal[idxst][j].replace('_', ' '), rotation=90,
                                             va='bottom')
                                text_tag_flag = 0
                                prev_text_XKey = XKey_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]

            if key == 'VIS' or key == 'VIS2' or key == 'TF2':
                if useStations==True:
                    axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                              sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
                else:
                    axlabel = ""
            elif key == 'T3':
                if useStations==True:
                    axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                              sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                              sta_names[sta_indices[i, 2] == dic['STA_INDEX']][0]
                else:
                    axlabel = ""
            elif key == 'FLUX':
                if useStations==True:
                    axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0]
                else:
                    axlabel = ""
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            if datatype == 'VIS2' or datatype == 'TF2':
                if ylog == True:
                    axs1[i].set_ylim([0.3, 1.1])
                else:
                    axs1[i].set_ylim([-0.1, 1.1])
            else:
                try:
                    axs1[i].set_ylim([np.nanmin([np.nanmin(arr_cal),np.nanmin(arr_sci)]), np.nanmax([np.nanmax(arr_cal),np.nanmax(arr_sci)])])
                except:
                    pass
            ylabel = datatype
            if showvis == True:
                if datatype == 'VIS2':
                    ylabel = 'VIS'
                    datatype = 'VIS'
                elif datatype == 'TF2':
                    ylabel = 'TF'
                    datatype = 'TF'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel('$\mathrm{'+xaxis.replace("HIERARCH ESO ","").replace(" ","\_")+'}$')

            if xlog == True:
                axs1[i].set_xscale('log')
            if ylog == True:
                axs1[i].set_yscale('log')
            #if ylog == True or xlog == True:
            #    for axis in [axs1[i].xaxis, axs1[i].yaxis]:
            #        axis.set_major_formatter(LogFormatter(minor_thresholds=[1,0.5,0.2],labelOnlyBase=False))
            #        axis.set_minor_formatter(LogFormatter(minor_thresholds=[1,0.5,0.2],labelOnlyBase=False))
                    # axis.set_scientific(False)

        plt.suptitle('$\mathrm{'+datatype+'\ vs.\ '+xaxis.replace("HIERARCH ESO ","").replace(" ","\_")+'}$')

        #plt.tight_layout()
        plt.show()

        plt.savefig(home+"/"+datatype+'_f_'+xaxis.replace(" ","_")+'.png')

###############################################################################
# This function shows the selected oifits data (flux, visibility,
# closure phase etc.)  as a function of time. It reads data from
# multiple oifits files in a given directory.  The data to be plotted
# can be filtered with the filter_oi_list function.
#
# Example usage:
# filtered_list_of_dicts = filter_oi_list(list_of_dicts,
# dates=["2018-03-14"], bands=['LM'], spectral_resolutions=['MED'],
# DIT_range=[0,0.2], targets=['l pup'])
# show_oi_vs_time(filtered_list_of_dicts, [3.5, 3.95], key="VIS2",
# datatype='VIS2') #[3.5, 3.95] [10.2,10.9]
#
def show_oi_vs_time(list_of_dicts, wlenRange, key="VIS2",subplotList=None,
                    datatype="VIS2",showvis=False,plot_errorbars=True,useStations=True, sciColor='red',calColor='blue'):
    print("Starting show_oi_vs_time...")




    # check if list is not empty:
    if list_of_dicts:
        target_names_cal = []
        MJD_arr_cal      = []
        arr_cal          = []
        err_arr_cal      = []
        sta_index_cal    = []

        target_names_sci = []
        MJD_arr_sci      = []
        arr_sci          = []
        err_arr_sci      = []
        sta_index_sci    = []

        for dic in list_of_dicts:
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            category = dic['CATEGORY'].lower()
            if 'cal' in category:
                try:
                    datay    = np.array(dic[key][datatype])
                    datayerr = np.array(dic[key][datatype+'ERR'])
                    datax    = np.array(dic[key]["TIME"])
                    n_rows   = datay.shape[0]
                    # print datay.shape
                    for i in range(n_rows):
                        arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                        err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                        MJD_arr_cal.append(datax[i])
                        target_names_cal.append(dic['TARGET'])
                        if key == 'FLUX':
                            sta_index = dic[key]['STA_INDEX'][i]
                            sta_index_cal.append([sta_index])
                        else:
                            sta_index = np.sort(dic[key]['STA_INDEX'][i])
                            sta_index_cal.append(sta_index)
                        # print dic[key]['STA_INDEX'][i]
                except:
                    print(dic['TARGET'], dic['DATEOBS'], "No CAL data found.")
            elif 'sci' in category:
                try:
                    datay    = np.array(dic[key][datatype])
                    datayerr = np.array(dic[key][datatype+'ERR'])
                    datax    = np.array(dic[key]["TIME"])
                    n_rows   = datay.shape[0]
                    # print datay.shape
                    for i in range(n_rows):
                        arr_sci.append(robust.mean(datay[i, wlenRange_idx]))
                        err_arr_sci.append(robust.mean(datayerr[i, wlenRange_idx]))
                        MJD_arr_sci.append(datax[i])
                        target_names_sci.append(dic['TARGET'])
                        if useStations==True:
                            if key == 'FLUX':
                                sta_index = dic[key]['STA_INDEX'][i]
                                sta_index_sci.append([sta_index])
                            else:
                                sta_index = np.sort(dic[key]['STA_INDEX'][i])
                                sta_index_sci.append(sta_index)
                        else:
                            sta_index = i%6
                            sta_index_cal.append([sta_index])
                except:
                    print(dic['TARGET'], dic['DATEOBS'], "No SCI data found.")

            if useStations==True:
                sta_names = dic['STA_NAME']
            else:
                sta_names = ""

        target_names_cal = np.array(target_names_cal)
        MJD_arr_cal      = np.array(MJD_arr_cal)
        arr_cal          = np.array(arr_cal)
        err_arr_cal      = np.array(err_arr_cal)
        sta_index_cal    = np.array(sta_index_cal)

        target_names_sci = np.array(target_names_sci)
        MJD_arr_sci      = np.array(MJD_arr_sci)
        arr_sci          = np.array(arr_sci)
        err_arr_sci      = np.array(err_arr_sci)
        sta_index_sci    = np.array(sta_index_sci)

        sta_indices = np.unique(sta_index_cal, axis=0)
        if key == 'VIS' or key == 'VIS2' or key == 'TF2':
            n_max_config = np.nanmax([6, sta_indices.shape[0]])
            n_plot_rows = 3
        elif key == 'T3' or key == 'FLUX':
            n_max_config = np.nanmax([4, sta_indices.shape[0]])
            n_plot_rows = 2
        # print sta_indices.shape

        if len(MJD_arr_cal) > 0 and len(MJD_arr_sci) > 0:
            MJD_range = [np.nanmin([np.nanmin(MJD_arr_cal), np.nanmin(MJD_arr_sci)]),
                         np.nanmax([np.nanmax(MJD_arr_cal), np.nanmax(MJD_arr_sci)])]
        elif len(MJD_arr_sci) > 0:
            MJD_range = [np.nanmin(MJD_arr_sci), np.nanmax(MJD_arr_sci)]
        elif len(MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(MJD_arr_cal), np.nanmax(MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0


        if not(subplotList):
            fig1, axs1 = plt.subplots(n_plot_rows, 2, figsize=(15, 16), sharex=True, sharey=True)
            axs1 = axs1.ravel()
        else:
            axs1 = subplotList






        for i in range(n_max_config):
            # print i
            if datatype == 'DPHI' or datatype == 'CLOS':
                axs1[i + 0].plot(MJD_range, [0.0, 0.0], '-', color='gray', lw=1.5)
            if len(sta_index_cal) > 0:
                idxst = np.all(sta_index_cal == sta_indices[i], axis=1)
                if len(arr_cal[idxst]) > 0:
                    label = datatype +' cal'
                    if plot_errorbars == True:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' cal'
                                elif key == 'TF2':
                                    label = 'TF' + ' cal'
                                axs1[i].errorbar(MJD_arr_cal[idxst], np.sqrt(arr_cal[idxst]),
                                                 yerr=0.5 * err_arr_cal[idxst] / np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color=calColor, elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_cal[idxst], arr_cal[idxst],
                                             yerr=err_arr_cal[idxst],
                                             fmt='o', color=calColor, elinewidth=1.0,
                                             label=label)
                    else:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_cal[idxst], np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color=calColor, elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_cal[idxst], arr_cal[idxst],
                                             fmt='o', color=calColor, elinewidth=1.0,
                                             label=label)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                ymin, ymax = axs1[i + 0].get_ylim()
                                axs1[i].text(MJD_arr_cal[idxst][j], ymax * 1.05,
                                             target_names_cal[idxst][j].replace('_', ' '), rotation=90,
                                             va='bottom')
                                text_tag_flag = 0
                                prev_text_MJD = MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]
            if len(sta_index_sci) > 0:
                label = datatype +' sci'
                idxst = np.all(sta_index_sci == sta_indices[i], axis=1)
                if len(arr_sci[idxst]) > 0:
                    if plot_errorbars == True:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_sci[idxst], np.sqrt(arr_sci[idxst]),
                                                 yerr=0.5 * err_arr_sci[idxst] / np.sqrt(arr_sci[idxst]),
                                                 fmt='o', color=sciColor, elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_sci[idxst], arr_sci[idxst], yerr=err_arr_sci[idxst],
                                                 fmt='o', color=sciColor, elinewidth=1.0,
                                                 label=label)
                    else:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_sci[idxst], np.sqrt(arr_sci[idxst]),
                                                 fmt='o', color=sciColor, elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_sci[idxst], arr_sci[idxst],
                                         fmt='o', color=sciColor, elinewidth=1.0,
                                         label=label)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if MJD_arr_sci[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_sci[idxst][j]):
                                ymin, ymax = axs1[i + 0].get_ylim()
                                axs1[i].text(MJD_arr_sci[idxst][j], ymax*1.05, target_names_sci[idxst][j].replace('_', ' '),
                                             rotation=90, va='bottom', color='darkred')
                                text_tag_flag = 0
                                prev_text_MJD = MJD_arr_sci[idxst][j]
                                prev_target_name = target_names_sci[idxst][j]
            if key == 'VIS' or key == 'VIS2' or key == 'TF2':
                if useStations==True:
                    axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
                else:
                    axlabel = ""
            elif key == 'T3':
                if useStations==True:
                    axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                              sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                              sta_names[sta_indices[i, 2] == dic['STA_INDEX']][0]
                else:
                    axlabel = ""
            elif key == 'FLUX':
                axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            if datatype == 'VIS2' or datatype == 'TF2':
                axs1[i].set_ylim([-0.1, 1.1])
            else:
                try:
                    axs1[i].set_ylim([np.nanmin([np.nanmin(arr_cal),np.nanmin(arr_sci)]), np.nanmax([np.nanmax(arr_cal),np.nanmax(arr_sci)])])
                except:
                    pass
            ylabel = datatype
            if showvis == True:
                if datatype == 'VIS2':
                    ylabel = 'VIS'
                    datatype = 'VIS'
                elif datatype == 'TF2':
                    ylabel = 'TF'
                    datatype = 'TF'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel('$\mathrm{MJD}$')
        plt.suptitle('$\mathrm{'+datatype+'\ vs.\ time}$')
        # fig1.subplots_adjust(hspace=0, wspace=0)
        # for i in range(4):
        #     plt.setp(axs1[i].get_xticklabels(), visible=False)
        #     x_axis = axs1[i].axes.get_xaxis()
        #     x_axis.get_label().set_visible(False)
        # for i in range(1, 6, 2):
        #     plt.setp(axs1[i].get_yticklabels(), visible=False)
        #     y_axis = axs1[i].axes.get_yaxis()
        #     y_axis.get_label().set_visible(False)
        if not(subplotList):
            plt.tight_layout()
        if not(subplotList):
            plt.show()


###############################################################################
# showvis: if True, plot visibilities (V) instead of V^2: V is calculated from V^2 (not from the VISAMP table)
def show_vis2_tf2_vs_time(list_of_dicts, wlenRange, showvis=False, saveplots=False, output_path="",plot_errorbars=True):
    starFluxN={'Teta Pyx':27, 'psi Vir':44,'tet Cen':48, \
          'del Vir':160,'eps Sco':52,'del Oph':110, \
          'FZ Lib':17,'del Sco':15,'HD171094':34, \
          'e Aql':27,'29 Cap':31,'alf Ant':19, \
          'HD99333':20,'C01 Cen':22,'alpha Arae':10, \
          'HD138505':18,'HD142527':12,'HD161560':0.4, \
          'HD171094':34,'AV Mic':18,'HD181041':6,  \
          'HD138492':8,'HD147929':8,'HD148255':3, \
          'G Sco':26,'HD156936':3,'HD161849':8, \
          'HD165413':5,'HD186765':9,'Mu Hya':30, \
          'CH Vir':20,'HD 126111':6,'LY TrA':7, \
          'ET Vir':29,'H Sco':30,'HD177507':17, \
          'V345 Tel':9,'RW Lup':18,'HD 150798':144, \
          'BM Sco':85,'RX Tel':82,'HD328913':100, \
          'Eps Oph':16,'HD182669':2.5,'Nu Hya':30}    # check if list is not empty:
    starFluxL={'HD138492':49,'HD147929':49,'HD148255':19, \
               'HD156936':23,'HD161849':33, 'HD165413':35, \
               'HD186765':39,'HD182669':16,'Teta Pyx':187, \
               'tet Cen':390,'Tet Cen':390,'alf Ant':129,'HD99333':112, \
               'C01 Cen':127,'HD138505':105,'HD171094':85, \
               'e Aql':152,'29 Cap':198,'AV Mic':113,'psi Vir':280, \
               'eps Sco':381,'ET Vir':169,'FZ Lib':141,'Eps Oph':116, \
               'HD142527':6,'RW Lup':42,'H Sco':154,'Mu Hya':219, \
               'Nu Hya':236,'CH Vir':63,'HD 126111':43,'LY TrA':35, \
               'HD177507':101,'V345 Tel':154,'HD121730':3,'HD132602':1, \
               'HD142277':3,'HD145278':0.5,'HD165368':50,'HD181925':50, \
               'HD189140':50,'HD 902':10,'TW PsA':10,'nu Phe':10,'HD218619':20, \
               'HD183925':21,'V4443 Sgr':24,'HD217826':5,'HD224936':5, \
               'HD6080':5,'HD216988':3,'HD224644':3,'HD4293':3,'HD457':1, \
               'tet Gru':1,'Pi Hya':157,'V1064 Sco':50,'Ome Lup':64, \
               'HD167818':-1,'B Sgr':54,'HD189831':50,'HD943':-1,'g Aqr':42}    # check if list is not empty:
    if (wlenRange[0] > 6):
        starFlux=starFluxN
    else:
        starFlux=starFluxL
    # check if list is not empty:
    if list_of_dicts:
        # colors: BCD: out-out, in-in, in-out, out-in
        BCD_configs   = np.array([[0, 0], [1, 1], [1, 0], [0, 1]])
        BCD_labels    = np.array(['OUT-OUT', 'IN-IN', 'IN-OUT', 'OUT-IN'])
        BCD_markers   = np.array(['o', 's', 'd', 'p'])
        V2_cal_colors = np.array(['palegoldenrod', 'khaki', 'navajowhite', 'moccasin'])
        V2_colors     = np.array(['red', 'orange', 'salmon', 'pink'])
        TF2_colors    = np.array(['darkseagreen', 'yellowgreen', 'olivedrab', 'darkkhaki'])

        target_names_cal = []
        V2_BCD_arr_cal   = []
        V2_MJD_arr_cal   = []
        V2_arr_cal       = []
        V2err_arr_cal    = []
        V2_sta_index_cal = []

        target_names_CP_cal = []
        CP_BCD_arr_cal      = []
        CP_MJD_arr_cal      = []
        CP_arr_cal          = []
        CPerr_arr_cal       = []
        CP_sta_index_cal    = []

        target_names_TF2 = []
        TF2_BCD_arr   = []
        TF2_MJD_arr   = []
        TF2_arr       = []
        TF2err_arr    = []
        TF2_sta_index = []

        target_names = []
        V2_BCD_arr   = []
        V2_MJD_arr   = []
        V2_arr       = []
        V2err_arr    = []
        V2_sta_index = []

        target_names_CP = []
        CP_BCD_arr      = []
        CP_MJD_arr      = []
        CP_arr          = []
        CPerr_arr       = []
        CP_sta_index    = []

        for dic in list_of_dicts:
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            if sum(wlenRange_idx) > 0:
                category = dic['CATEGORY'].lower()
                if 'cal' in category:
                    try:
                        datay = np.array(dic['VIS2']['VIS2'])
                        datayerr = np.array(dic['VIS2']['VIS2ERR'])
                        datax = np.array(dic['VIS2']["TIME"])
                        n_rows = datay.shape[0]
                        # print datay.shape
                        for i in range(n_rows):
                            # print robust.mean(datay[i, wlenRange_idx])
                            if dic['BCD1NAME']   == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME']   == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            V2_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            V2err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V2_BCD_arr_cal.append([BCD1, BCD2])
                            V2_MJD_arr_cal.append(datax[i])
                            target_names_cal.append(dic['TARGET'])
                            sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                            V2_sta_index_cal.append(sta_index)
                    except:
                        print(dic['TARGET'], dic['DATEOBS'], "No CAL VIS2 data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            CP_BCD_arr_cal.append([BCD1, BCD2])
                            CP_MJD_arr_cal.append(datax[i])
                            target_names_CP_cal.append(dic['TARGET'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index_cal.append(sta_index)
                    except:
                        print(dic['TARGET'], dic['DATEOBS'], "No CAL CP data found.")
                    try:
                        datay = np.array(dic['TF2']['TF2'])
                        datayerr = np.array(dic['TF2']['TF2ERR'])
                        datax = np.array(dic['TF2']["TIME"])
                        n_rows = datay.shape[0]
                        # print datay.shape
                        for i in range(n_rows):
                            # print robust.mean(datay[i, wlenRange_idx])
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            TF2_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            TF2err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            TF2_BCD_arr.append([BCD1, BCD2])
                            TF2_MJD_arr.append(datax[i])
                            target_names_TF2.append(dic['TARGET'])
                            sta_index = np.sort(dic['TF2']['STA_INDEX'][i])
                            TF2_sta_index.append(sta_index)
                    except:
                        print(dic['TARGET'], dic['DATEOBS'], "No CAL TF2 data found.")
                if 'sci' in category:
                    try:
                        datay = np.array(dic['VIS2']['VIS2'])
                        datayerr = np.array(dic['VIS2']['VIS2ERR'])
                        datax = np.array(dic['VIS2']["TIME"])
                        n_rows = datay.shape[0]
                        #print (datay.shape)
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            #print i
                            #print datay[i, wlenRange_idx]
                            V2_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            V2err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V2_BCD_arr.append([BCD1, BCD2])
                            V2_MJD_arr.append(datax[i])
                            target_names.append(dic['TARGET'])
                            sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                            V2_sta_index.append(sta_index)
                    except:
                        print(dic['TARGET'], dic['DATEOBS'], "No SCI VIS2 data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            target_names_CP.append(dic['TARGET'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index.append(sta_index)
                            CP_BCD_arr.append([BCD1, BCD2])
                            CP_MJD_arr.append(datax[i])
                    except:
                        print(dic['TARGET'], dic['DATEOBS'], "No SCI CP data found.")
            else:
                print("Wavelength out of range.")

        sta_names = dic['STA_NAME']

        target_names_cal = np.array(target_names_cal)
        V2_BCD_arr_cal   = np.array(V2_BCD_arr_cal)
        V2_MJD_arr_cal   = np.array(V2_MJD_arr_cal)
        V2_arr_cal       = np.array(V2_arr_cal)
        V2err_arr_cal    = np.array(V2err_arr_cal)
        V2_sta_index_cal = np.array(V2_sta_index_cal)

        target_names_CP_cal = np.array(target_names_CP_cal)
        CP_BCD_arr_cal      = np.array(CP_BCD_arr_cal)
        CP_MJD_arr_cal      = np.array(CP_MJD_arr_cal)
        CP_arr_cal          = np.array(CP_arr_cal)
        CPerr_arr_cal       = np.array(CPerr_arr_cal)
        CP_sta_index_cal    = np.array(CP_sta_index_cal)

        target_names_TF2 = np.array(target_names_TF2)
        TF2_BCD_arr   = np.array(TF2_BCD_arr)
        TF2_MJD_arr   = np.array(TF2_MJD_arr)
        TF2_arr       = np.array(TF2_arr)
        TF2err_arr    = np.array(TF2err_arr)
        TF2_sta_index = np.array(TF2_sta_index)

        target_names = np.array(target_names)
        V2_BCD_arr   = np.array(V2_BCD_arr)
        V2_MJD_arr   = np.array(V2_MJD_arr)
        V2_arr       = np.array(V2_arr)
        V2err_arr    = np.array(V2err_arr)
        V2_sta_index = np.array(V2_sta_index)

        target_names_CP = np.array(target_names_CP)
        CP_BCD_arr      = np.array(CP_BCD_arr)
        CP_MJD_arr      = np.array(CP_MJD_arr)
        CP_arr          = np.array(CP_arr)
        CPerr_arr       = np.array(CPerr_arr)
        CP_sta_index    = np.array(CP_sta_index)

        print(V2_MJD_arr_cal)
        print(V2_MJD_arr)

        if len(V2_sta_index_cal) > 0:
            sta_indices = np.unique(V2_sta_index_cal, axis=0)
        elif len(V2_sta_index) > 0:
            sta_indices = np.unique(V2_sta_index, axis=0)
        else:
            print ("Data arrays empty. Quitting.")
            return
        n_max_config = np.nanmax([6, sta_indices.shape[0]])
        # print sta_indices.shape

        if len(V2_MJD_arr_cal) > 0 and len(V2_MJD_arr) > 0:
            MJD_range = [np.nanmin([np.nanmin(V2_MJD_arr_cal), np.nanmin(V2_MJD_arr)]),
                         np.nanmax([np.nanmax(V2_MJD_arr_cal), np.nanmax(V2_MJD_arr)])]
        elif len(V2_MJD_arr) > 0:
            MJD_range = [np.nanmin(V2_MJD_arr), np.nanmax(V2_MJD_arr)]
        elif len(V2_MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(V2_MJD_arr_cal), np.nanmax(V2_MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        print(MJD_range)
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0
        fig1, axs1 = plt.subplots(3, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        text_y = 1.55
        for i in range(n_max_config):
            if len(V2_sta_index_cal) > 0:
                idxst = np.all(V2_sta_index_cal == sta_indices[i], axis=1)
                if len(V2_arr_cal[idxst]) > 0:
                    if showvis == True:
                        label = 'V cal '
                    else:
                        label = 'V2 cal '
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V2_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(V2_arr_cal[cidxst]) > 0:
                            if plot_errorbars == True:
                                if showvis == True:
                                    axs1[i].errorbar(V2_MJD_arr_cal[cidxst], np.sqrt(np.abs(V2_arr_cal[cidxst])),
                                                     yerr=0.5 * V2err_arr_cal[cidxst] / np.sqrt(np.abs(V2_arr_cal[cidxst])),
                                                     fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(V2_MJD_arr_cal[cidxst], V2_arr_cal[cidxst],
                                                     yerr=V2err_arr_cal[cidxst],
                                                     fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis == True:
                                    axs1[i].errorbar(V2_MJD_arr_cal[cidxst], np.sqrt(np.abs(V2_arr_cal[cidxst])),
                                                     fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(V2_MJD_arr_cal[cidxst], V2_arr_cal[cidxst],
                                                     fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if V2_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                if (target_names_cal[idxst][j] in starFlux) :
                                    axs1[i].text(V2_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' ')+ \
                                                 " ("+np.str(starFlux[target_names_cal[idxst][j]])+"Jy)", rotation=90, \
                                                 va='bottom',fontsize=10)
                                else :
                                    axs1[i].text(V2_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' '), rotation=90, \
                                                 va='bottom',fontsize=10)

                                text_tag_flag = 0
                                prev_text_MJD = V2_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]
            if len(TF2_sta_index) > 0:
                if showvis == True:
                    label = 'TF '
                else:
                    label = 'TF2 '
                idxst = np.all(TF2_sta_index == sta_indices[i], axis=1)
                if len(TF2_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(TF2_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(TF2_arr[cidxst]) > 0:
                            if plot_errorbars == True:
                                if showvis == True:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     yerr=0.5 * TF2err_arr[cidxst] / np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                    z=np.polyfit(TF2_MJD_arr[cidxst]-np.min(TF2_MJD_arr[cidxst]), np.sqrt(np.abs(TF2_arr[cidxst])),3)
                                    p=np.poly1d(z)
                                    x=(np.max(TF2_MJD_arr[cidxst])-np.min(TF2_MJD_arr[cidxst]))*np.arange(100)/100.
                                    axs1[i].plot(np.min(TF2_MJD_arr[cidxst])+x,p(x),color=TF2_colors[j])
                                    val=np.sqrt(np.abs(TF2_arr[cidxst]))-p(TF2_MJD_arr[cidxst]-np.min(TF2_MJD_arr[cidxst]))
                                    valmed=np.median(val)
                                    print np.median(np.abs(val-valmed)),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
                                else:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], TF2_arr[cidxst], yerr=TF2err_arr[cidxst],
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis == True:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], TF2_arr[cidxst],
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
            if len(V2_sta_index) > 0:
                if showvis == True:
                    label = 'V sci '
                else:
                    label = 'V2 sci '
                idxst = np.all(V2_sta_index == sta_indices[i], axis=1)
                if len(V2_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V2_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(V2_arr[cidxst]) > 0:
                            if plot_errorbars == True:
                                if showvis == True:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], np.sqrt(V2_arr[cidxst]),
                                                     yerr=0.5 * V2err_arr[cidxst] / np.sqrt(V2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], V2_arr[cidxst], yerr=V2err_arr[cidxst],
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis == True:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], np.sqrt(V2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], V2_arr[cidxst],
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
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

            axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs1[i].set_ylim([-0.1, 1.5])
            if showvis == True:
                ylabel = '$V$'
            else:
                ylabel = '$V^2$'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel('$\mathrm{MJD}$')
        if showvis == True:
            plt.suptitle('$V\mathrm{\ vs.\ time}$')
        else:
            plt.suptitle('$V^2\mathrm{\ vs.\ time}$')
        fig1.subplots_adjust(hspace=0, wspace=0)
        for i in range(4):
            plt.setp(axs1[i].get_xticklabels(), visible=False)
            x_axis = axs1[i].axes.get_xaxis()
            x_axis.get_label().set_visible(False)
        for i in range(1, 6, 2):
            plt.setp(axs1[i].get_yticklabels(), visible=False)
            y_axis = axs1[i].axes.get_yaxis()
            y_axis.get_label().set_visible(False)
        # plt.tight_layout()
        if saveplots == True:
            if showvis == True:
                label = '_VIS_TF'
            else:
                label = '_VIS2_TF2'
            fig1.savefig(output_path + label + '.png', dpi=150)
            fig1.savefig(output_path + label + '.eps', format='eps', dpi=300)
            plt.close(fig1)

        if len(CP_sta_index_cal) > 0:
            CP_sta_indices = np.unique(CP_sta_index_cal, axis=0)
        else:
            CP_sta_indices = np.unique(CP_sta_index, axis=0)
        # print CP_sta_indices
        n_max_config = np.nanmax([4, CP_sta_indices.shape[0]])

        fig2, axs = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)
        axs = axs.ravel()
        text_y = 60
        for i in range(n_max_config):
            axs[i + 0].plot(MJD_range, [0.0, 0.0], '-', color='gray', lw=1.5)
            if len(CP_sta_index_cal) > 0:
                idxst = np.all(CP_sta_index_cal == CP_sta_indices[i], axis=1)
                if len(CP_arr_cal[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr_cal[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst], yerr=CPerr_arr_cal[cidxst],
                                                fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                label='CP cal ' + BCD_labels[j])
                                z=np.polyfit(CP_MJD_arr_cal[cidxst]-np.min(CP_MJD_arr_cal[cidxst]), CP_arr_cal[cidxst],1)
                                p=np.poly1d(z)
                                x=(np.max(CP_MJD_arr_cal[cidxst])-np.min(CP_MJD_arr_cal[cidxst]))*np.arange(100)/100.
                                axs[i + 0].plot(np.min(CP_MJD_arr_cal[cidxst])+x,p(x),color=V2_cal_colors[j])
                                val=np.sqrt(np.abs(CP_arr_cal[cidxst]))-p(CP_MJD_arr_cal[cidxst]-np.min(CP_MJD_arr_cal[cidxst]))
                                valmed=np.median(val)
                                print "CP",np.median(np.abs(val-valmed)),BCD_labels[j],sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0]+ ' - ' + sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst],
                                                    fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                    label='CP cal ' + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP_cal[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                if (target_names_CP_cal[idxst][j] in starFlux):
                                    axs[i + 0].text(CP_MJD_arr_cal[idxst][j], ymax * 1.05, \
                                                    target_names_CP_cal[idxst][j].replace('_', ' ')+ \
                                                    " ("+np.str(starFlux[target_names_CP_cal[idxst][j]])+"Jy)", rotation=90,
                                                    va='bottom',fontsize=8)
                                else:
                                     axs[i + 0].text(CP_MJD_arr_cal[idxst][j], ymax * 1.05, \
                                                    target_names_CP_cal[idxst][j].replace('_', ' '), rotation=90,
                                                    va='bottom',fontsize=8)

                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_CP_cal[idxst][j]
            if len(CP_sta_index) > 0:
                idxst = np.all(CP_sta_index == CP_sta_indices[i], axis=1)
                if len(CP_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst], yerr=CPerr_arr[cidxst],
                                                fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                label='CP sci ' + BCD_labels[j])
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst],
                                                    fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                    label='CP sci ' + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                axs[i + 0].text(CP_MJD_arr[idxst][j], ymax * 1.05,
                                                target_names_CP[idxst][j].replace('_', ' '),
                                                rotation=90,
                                                va='bottom')
                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr[idxst][j]
                                prev_target_name = target_names_CP[idxst][j]
            axlabel = sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs[i + 0].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                            transform=axs[i + 0].transAxes, bbox=props)
            if i == 0:
                leg = axs[i + 0].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs[i + 0].set_ylabel('$CP\,\left(^\circ\\right)$')
            axs[i + 0].set_xlabel('$\mathrm{MJD}$')
        plt.suptitle('$CP\mathrm{\ vs.\ time}$')
        fig2.subplots_adjust(hspace=0, wspace=0)
        for i in range(2):
            plt.setp(axs[i + 0].get_xticklabels(), visible=False)
            x_axis = axs[i + 0].axes.get_xaxis()
            x_axis.get_label().set_visible(False)
        for i in range(1, 4, 2):
            plt.setp(axs[i + 0].get_yticklabels(), visible=False)
            y_axis = axs[i + 0].axes.get_yaxis()
            y_axis.get_label().set_visible(False)
        # plt.tight_layout()
        if saveplots == True:
            fig2.savefig(output_path + '_CP' + '.png', dpi=150)
            fig2.savefig(output_path + '_CP' + '.eps', format='eps', dpi=300)
            plt.close(fig2)
        else:
            plt.show()
        print ("Plots READY.")

###############################################################################

def open_oi_dir(input_dir, verbose=True):
    oifits_file_list = glob.glob(input_dir + '/*fits*')

    N_files = len(oifits_file_list)
    list_of_dicts = []
    for file in oifits_file_list:
        if "LAMP" not in file:
            dic = open_oi(file)
            if dic:
                if verbose==True:
                    print(dic['TARGET'] + " " + dic['DATEOBS'] + " " + dic['BAND'] + " " + dic['DISP'] + " " + str(dic['DIT']) + " " + dic['CATEGORY'])
                list_of_dicts.append(dic)

    print("Done!")
    return list_of_dicts

###############################################################################
# dates = example format: ["2018-03-16"]
# bands = 'L','M','LM', 'N'
# spectral_resolutions: 'LOW','MED','HIGH'
# DIT_range: [min,max] (s)
# targets = []

def filter_oi_list(list_of_dicts, dates=[], bands=[], spectral_resolutions=[], DIT_range=[], targets=[], BCD1=[], BCD2=[], WLEN_range=[]):
    filtered_list_of_dicts = []
    if bands:
        # print  'old:',bands
        bands_new = []
        for i in range(len(bands)):
            if bands[i] == 'M':
                bands_new.append('LM')
            elif bands[i] == 'L':
                bands_new.append('LM')
            else:
                bands_new.append(bands[i])
                # print 'new: ', bands_new
    for dic in list_of_dicts:
        if dic:
            date = dic['DATEOBS'][0:10]
            #if dates:
            #    if date not in dates:
            #        continue
            if bands:
                if dic['BAND'] not in bands_new:
                    continue
            if BCD1:
                if dic['BCD1NAME'] not in BCD1:
                    continue
            if BCD2:
                if dic['BCD2NAME'] not in BCD2:
                    continue
            if spectral_resolutions:
                if dic['DISP'] not in spectral_resolutions:
                    continue
            if DIT_range:
                if not (dic['DIT'] >= DIT_range[0] and dic['DIT'] <= DIT_range[1]):
                    continue
            target = dic['TARGET']
            if targets:
                targets = [x.lower().replace("_", " ") for x in targets]
                target = target.lower().replace("_", " ")
                if target not in targets:
                    continue
            print("Selected: ", target, date, dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'])
            filtered_list_of_dicts.append(dic)

            # filter oifits files on wavelength
            if WLEN_range:
                wl1 = WLEN_range[0];
                wl2 = WLEN_range[1];
                dic['VIS2']['VIS2']      = dic['VIS2']['VIS2'][wl1:wl2]
                dic['VIS2']['VIS2ERR']   = dic['VIS2']['VIS2ERR'][wl1:wl2]
                dic['VIS']['VISAMP']     = dic['VIS']['VISAMP'][wl1:wl2]
                dic['VIS']['VISAMPERR']  = dic['VIS']['VISAMPERR'][wl1:wl2]
                dic['VIS']['DPHI']       = dic['VIS']['DPHI'][wl1:wl2]
                dic['VIS']['DPHIERR']    = dic['VIS']['DPHIERR'][wl1:wl2]
                dic['TF2']['TF2']        = dic['TF2']['TF2'][wl1:wl2]
                dic['TF2']['TF2ERR']     = dic['TF2']['TF2ERR'][wl1:wl2]
                dic['T3']['T3AMP']       = dic['T3']['T3AMP'][wl1:wl2]
                dic['T3']['T3AMPERR']    = dic['T3']['T3AMPERR'][wl1:wl2]
                dic['T3']['T3CLOS']      = dic['T3']['CLOS'][wl1:wl2]
                dic['T3']['T3CLOSERR']   = dic['T3']['CLOSERR'][wl1:wl2]
                dic['FLUX']['FLUX']      = dic['T3']['CLOS'][wl1:wl2]
                dic['FLUX']['T3CLOSERR'] = dic['T3']['CLOSERR'][wl1:wl2]

    return filtered_list_of_dicts

###############################################################################
# dates = example format: ["2018-03-16"]
# bands = 'L','M','LM', 'N'
# spectral_resolutions: 'LOW','MED','HIGH'
# DIT_range: [min,max] (s)
# targets = []
def filter_oi_list_night(list_of_dicts, dates='2000-01-01', bands=[], spectral_resolutions=[], DIT_range=[], targets=[]):
    t=Time(dates)
    mjd0=t.mjd
    filtered_list_of_dicts = []
    if bands:
        # print  'old:',bands
        bands_new = []
        for i in range(len(bands)):
            if bands[i] == 'M':
                bands_new.append('LM')
            elif bands[i] == 'L':
                bands_new.append('LM')
            else:
                bands_new.append(bands[i])
                # print 'new: ', bands_new
    for dic in list_of_dicts:
        if dic:
            date = dic['DATEOBS']
            tm   = Time(date)
            if (np.abs(tm.mjd-mjd0) > 0.5):
                continue
            if bands:
                if dic['BAND'] not in bands_new:
                    continue
            if spectral_resolutions:
                if dic['DISP'] not in spectral_resolutions:
                    continue
            if DIT_range:
                if not (dic['DIT'] >= DIT_range[0] and dic['DIT'] <= DIT_range[1]):
                    continue
            target = dic['TARGET']
            if targets:
                targets = [x.lower().replace("_", " ") for x in targets]
                target  = target.lower().replace("_", " ")
                if target not in targets:
                    continue
            print("Selected: ", target, date, dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'])
            filtered_list_of_dicts.append(dic)

    return filtered_list_of_dicts


###############################################################################
# Example code for TF and VIS plots.
# name_dir = r"D:\jvarga\Dokumentumok\MATISSE\data\OIFITS/"
# # outputdir = r"D:/jvarga/Dokumentumok/MATISSE/data/OIFITS/"
# list_of_dicts = open_oi_dir(name_dir+"2018-03-15")
# filtered_list_of_dicts = filter_oi_list(list_of_dicts,dates=["2018-03-15"],bands=['L'],spectral_resolutions=['LOW'],DIT_range=[0.045,0.055],targets=[])
# show_vis2_tf2_vs_time(filtered_list_of_dicts, wlenRange=[3.6,4.0], showvis=False, saveplots=False)

# dates=["2018-03-14","2018-03-15","2018-03-16","2018-03-11","2018-03-13"]
# bands=['L','M', 'N']
# sp_res = ['MED','LOW','HIGH']
# wlenRange = [[3.6,4.0],[4.6,4.8],[10.0,11.0]]
# DITs = np.array([0.2,0.02,0.05,0.1])
# dt = 0.001
# DITranges = np.transpose(np.array([DITs-dt,DITs+dt])).tolist()
# i=0
# j=0
# k=0
# l=0
# for i in range(len(dates)):
#     list_of_dicts = open_oi_dir(name_dir + dates[i])
#     for j in range(len(bands)):
#         for k in range(len(sp_res)):
#             for l in range(len(DITranges)):
#                 filtered_list_of_dicts = filter_oi_list(list_of_dicts,dates=[dates[i]],bands=[bands[j]],spectral_resolutions=[sp_res[k]],DIT_range=DITranges[l],targets=[])
#                 #show_oi_vs_time(filtered_list_of_dicts, [3.5, 3.95], key="VIS2", datatype='VIS2') #[3.5, 3.95] [10.2,10.9]
#                 print "Selected",len(filtered_list_of_dicts),"objects"
#                 fname = dates[i] + "_" + bands[j] + "_" + sp_res[k] + "_DIT" + ("%.2f"%(DITs[l])).replace('.','_')
#                 #show_oi_vs_wlen(filtered_list_of_dicts[3],datatype="VIS2")
#                 show_vis2_tf2_vs_time(filtered_list_of_dicts,wlenRange=wlenRange[j],showvis=False,saveplots=True,output_path=outputdir+fname) # wlenRange=[3.5, 3.95]);
#
# raise SystemExit

class oi_data_select_frame(wx.Frame):
    wl_min_def = np.array([1.175,1.5, 2.05,3.0,4.5, 7.8,3.0]) #um, J,H,K,L,M,N,LM
    wl_max_def = np.array([1.325,1.75,2.35,4.2,5.0,12.0,5.0 ]) #um, J,H,K,L,M,N,LM
    DIT       = 0.2
    DIT_range = 0.4
    wl_min = wl_min_def[6]
    wl_max = wl_max_def[6]
    date   = "2018-03-16"
    bands  = np.array(['L','M','LM','N'])
    spectral_resolutions = np.array(['LOW','MED','HIGH'])
    name_file       = ""
    name_dir        = ""
    target_selected = ""

    ###########################################################################

    def __init__(self, *args, **kwds):
        self.dic = {}
        self.list_of_dicts = [{}]
        self.filtered_list_of_dicts = [{}]

        # begin wxGlade: oi_data_select_frame.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((420, 559))

        self.statusbar = self.CreateStatusBar(1)
        self.statusbar.SetStatusText('')

        # Menu Bar
        #self.frame_menubar = wx.MenuBar()
        #wxglade_tmp_menu   = wx.Menu()
        #self.frame_menubar.Append(wxglade_tmp_menu, "File")
        #self.SetMenuBar(self.frame_menubar)
        # Menu Bar end
        self.panel = wx.Panel(self, wx.ID_ANY)
        self.btn_open_oifits = wx.Button(self.panel, 15, "Open OIFITS data")
        self.btn_load_oifits = wx.Button(self.panel, 14, "Load data")
        self.tb_path    = wx.TextCtrl(self.panel, wx.ID_ANY, "")
        self.tb_target  = wx.TextCtrl(self.panel, wx.ID_ANY, "")
        self.tb_date    = wx.TextCtrl(self.panel, wx.ID_ANY, self.date)
        self.tb_DIT     = wx.TextCtrl(self.panel, wx.ID_ANY, "%.2f"%(self.DIT))
        self.tb_DIT_range = wx.TextCtrl(self.panel, wx.ID_ANY, "%.3f"%(self.DIT_range))
        #self.cb_disp_LOW  = wx.CheckBox(self.panel, wx.ID_ANY, "LOW")
        #self.cb_disp_MED  = wx.CheckBox(self.panel, wx.ID_ANY, "MED")
        #self.cb_disp_HIGH = wx.CheckBox(self.panel, wx.ID_ANY, "HIGH")
        self.cb_disp_LOW  = wx.RadioButton(self.panel, wx.ID_ANY, "LOW", style=wx.RB_GROUP)
        self.cb_disp_MED  = wx.RadioButton(self.panel, wx.ID_ANY, "MED")
        self.cb_disp_HIGH = wx.RadioButton(self.panel, wx.ID_ANY, "HIGH")

        self.cb_pl_OI_t = wx.RadioButton(self.panel, wx.ID_ANY, "Plot OI vs. time", style=wx.RB_GROUP)
        self.cb_pl_OI_wl = wx.RadioButton(self.panel, wx.ID_ANY, "Plot OI vs. wl.")

        self.cb_b_J  = wx.CheckBox(self.panel, wx.ID_ANY, "J")
        self.cb_b_H  = wx.CheckBox(self.panel, wx.ID_ANY, "H")
        self.cb_b_K  = wx.CheckBox(self.panel, wx.ID_ANY, "K")
        self.cb_b_L  = wx.CheckBox(self.panel, wx.ID_ANY, "L")
        self.cb_b_M  = wx.CheckBox(self.panel, wx.ID_ANY, "M")
        self.cb_b_N  = wx.CheckBox(self.panel, wx.ID_ANY, "N")
        self.cb_b_LM = wx.CheckBox(self.panel, wx.ID_ANY, "LM")

        self.tb_wl_min = wx.TextCtrl(self.panel, wx.ID_ANY, "3.6")
        self.tb_wl_max = wx.TextCtrl(self.panel, wx.ID_ANY, "4.0")
        self.btn_def_wl  = wx.Button(self.panel, 13, "Default wl")

        self.cb_errorbar = wx.CheckBox(self.panel, wx.ID_ANY, "Plot errorbars")

        self.btn_CFXAMP = wx.Button(self.panel, 16, "CFXAMP")
        self.btn_CFX2   = wx.Button(self.panel, 17, "CFX2")
        self.btn_VISAMP = wx.Button(self.panel, 0, "VISAMP")
        self.btn_VIS2   = wx.Button(self.panel, 1, "VIS2")
        self.btn_VISPHI = wx.Button(self.panel, 2, "VISPHI")
        self.btn_VIS    = wx.Button(self.panel, 7, "VIS")
        self.btn_T3AMP  = wx.Button(self.panel, 4, "T3AMP")
        self.btn_TF2    = wx.Button(self.panel, 3, "TF2")
        self.btn_T3PHI  = wx.Button(self.panel, 6, "T3PHI")
        self.btn_TF     = wx.Button(self.panel, 8, "TF")
        self.btn_FLUX   = wx.Button(self.panel, 5, "FLUX")
        self.btn_VIS2_vs_freq = wx.Button(self.panel, 9, "VIS2 vs. spatial frequency")
        self.btn_VIS_vs_freq = wx.Button(self.panel, 10, "VIS vs. spatial frequency")
        self.btn_vis2_tf2_cp_vs_time = wx.Button(self.panel, 11, "VIS2, TF2, T3PHI vs. time")
        self.btn_vis_tf_cp_vs_time = wx.Button(self.panel, 12, "VIS, TF, T3PHI vs. time")

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    ###########################################################################

    def __set_properties(self):
        # begin wxGlade: oi_data_select_frame.__set_properties
        self.SetTitle("OIFITS plotter")
        self.tb_path.SetMinSize((134, 23))
        self.tb_wl_min.SetMinSize((60, 23))
        self.tb_wl_max.SetMinSize((60, 23))
        self.btn_VIS_vs_freq.SetMinSize((153, 26))
        self.btn_vis2_tf2_cp_vs_time.SetMinSize((153, 26))
        self.btn_vis_tf_cp_vs_time.SetMinSize((153, 26))
        # end wxGlade

        self.cb_disp_MED.SetValue(True)
        self.cb_b_LM.SetValue(True)
        self.cb_errorbar.SetValue(True)
        self.cb_pl_OI_wl.SetValue(True)
        self.Bind(wx.EVT_BUTTON, self.OnButtonClicked)

    ###########################################################################

    def __do_layout(self):
        # begin wxGlade: oi_data_select_frame.__do_layout
        sizer_1      = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.FlexGridSizer(1, 2, 0, 0)
        sizer_4      = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Band"), wx.VERTICAL)
        sizer_9      = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Plot special"), wx.VERTICAL)
        sizer_10     = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Special options"), wx.VERTICAL)
        sizer_5      = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_5 = wx.FlexGridSizer(0, 3, 0, 0)
        grid_sizer_4 = wx.FlexGridSizer(3, 3, 0, 0)
        sizer_7      = wx.BoxSizer(wx.VERTICAL)
        sizer_16     = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Plot"), wx.HORIZONTAL)
        grid_sizer_6 = wx.FlexGridSizer(6, 2, 0, 0)
        sizer_11     = wx.FlexGridSizer(1, 2, 0, 0)
        sizer_15     = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Plot options"), wx.VERTICAL)
        sizer_12     = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Spectral res."), wx.VERTICAL)
        sizer_13     = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Target, date and DIT"), wx.HORIZONTAL)
        grid_sizer_8 = wx.FlexGridSizer(4, 3, 0, 0)
        sizer_3      = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY,
                                                      "Input path"), wx.VERTICAL)
        grid_sizer_3 = wx.FlexGridSizer(0, 2, 0, 0)
        sizer_14     = wx.BoxSizer(wx.HORIZONTAL)
        sizer_14.Add(self.btn_open_oifits, 0, 0, 0)
        sizer_14.Add(self.btn_load_oifits, 0, 0, 0)
        sizer_3.Add(sizer_14, 1, wx.EXPAND, 0)
        label_1      = wx.StaticText(self.panel, wx.ID_ANY, "File/folder")
        grid_sizer_3.Add(label_1, 0, wx.ALL, 3)
        grid_sizer_3.Add(self.tb_path, 0, 0, 0)
        sizer_3.Add(grid_sizer_3, 1, wx.ALL | wx.EXPAND, 3)
        sizer_7.Add(sizer_3, 1, wx.EXPAND, 0)
        label_29     = wx.StaticText(self.panel, wx.ID_ANY, "Target")
        grid_sizer_8.Add(label_29, 0, 0, 0)
        grid_sizer_8.Add(self.tb_target, 0, 0, 0)
        grid_sizer_8.Add((0, 0), 0, 0, 0)
        label_23     = wx.StaticText(self.panel, wx.ID_ANY, "Date")
        grid_sizer_8.Add(label_23, 0, 0, 0)
        grid_sizer_8.Add(self.tb_date, 0, 0, 0)
        label_24     = wx.StaticText(self.panel, wx.ID_ANY, "")
        grid_sizer_8.Add(label_24, 0, 0, 0)
        label_25     = wx.StaticText(self.panel, wx.ID_ANY, "DIT")
        grid_sizer_8.Add(label_25, 0, 0, 0)
        grid_sizer_8.Add(self.tb_DIT, 0, 0, 0)
        label_26     = wx.StaticText(self.panel, wx.ID_ANY, "s")
        grid_sizer_8.Add(label_26, 0, 0, 0)
        label_27     = wx.StaticText(self.panel, wx.ID_ANY, "DITrange")
        grid_sizer_8.Add(label_27, 0, 0, 0)
        grid_sizer_8.Add(self.tb_DIT_range, 0, 0, 0)
        label_28     = wx.StaticText(self.panel, wx.ID_ANY, "s")
        grid_sizer_8.Add(label_28, 0, 0, 0)
        sizer_13.Add(grid_sizer_8, 1, wx.ALL | wx.EXPAND, 3)
        sizer_7.Add(sizer_13, 1, wx.EXPAND, 0)
        sizer_12.Add(self.cb_disp_LOW, 0, 0, 0)
        sizer_12.Add(self.cb_disp_MED, 0, 0, 0)
        sizer_12.Add(self.cb_disp_HIGH, 0, 0, 0)
        sizer_11.Add(sizer_12, 1, wx.ALL | wx.EXPAND, 3)
        sizer_15.Add(self.cb_pl_OI_t, 0, 0, 0)
        sizer_15.Add(self.cb_pl_OI_wl, 0, 0, 0)
        sizer_15.Add((0, 0), 0, 0, 0)
        sizer_11.Add(sizer_15, 1, wx.ALL | wx.EXPAND, 3)
        sizer_7.Add(sizer_11, 1, wx.EXPAND, 0)
        grid_sizer_6.Add(self.btn_VISAMP, 0, 0, 0)
        grid_sizer_6.Add(self.btn_VIS2, 0, 0, 0)
        grid_sizer_6.Add(self.btn_VISPHI, 0, 0, 0)
        grid_sizer_6.Add(self.btn_VIS, 0, 0, 0)
        grid_sizer_6.Add(self.btn_T3AMP, 0, 0, 0)
        grid_sizer_6.Add(self.btn_TF2, 0, 0, 0)
        grid_sizer_6.Add(self.btn_T3PHI, 0, 0, 0)
        grid_sizer_6.Add(self.btn_TF, 0, 0, 0)
        grid_sizer_6.Add(self.btn_CFX2, 0, 0, 0)
        grid_sizer_6.Add(self.btn_FLUX, 0, 0, 0)
        grid_sizer_6.Add(self.btn_CFXAMP, 0, 0, 0)
        grid_sizer_6.Add((0, 0), 0, 0, 0)
        sizer_16.Add(grid_sizer_6, 1, wx.EXPAND, 0)
        sizer_7.Add(sizer_16, 1, wx.EXPAND, 0)
        grid_sizer_1.Add(sizer_7, 1, wx.EXPAND, 0)
        grid_sizer_4.Add(self.cb_b_J, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_H, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_K, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_L, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_M, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_N, 0, 0, 0)
        grid_sizer_4.Add(self.cb_b_LM, 0, 0, 0)
        grid_sizer_4.Add((0, 0), 0, 0, 0)
        grid_sizer_4.Add((0, 0), 0, 0, 0)
        sizer_4.Add(grid_sizer_4, 1, wx.EXPAND, 0)
        label_7      = wx.StaticText(self.panel, wx.ID_ANY, "Wl min")
        grid_sizer_5.Add(label_7, 0, 0, 0)
        grid_sizer_5.Add(self.tb_wl_min, 0, 0, 0)
        label_9      = wx.StaticText(self.panel, wx.ID_ANY, "um")
        grid_sizer_5.Add(label_9, 0, 0, 0)
        label_8      = wx.StaticText(self.panel, wx.ID_ANY, "Wl max")
        grid_sizer_5.Add(label_8, 0, 0, 0)
        grid_sizer_5.Add(self.tb_wl_max, 0, 0, 0)
        label_10     = wx.StaticText(self.panel, wx.ID_ANY, "um")
        grid_sizer_5.Add(label_10, 0, 0, 0)
        sizer_5.Add(grid_sizer_5, 1, wx.EXPAND, 0)
        sizer_5.Add(self.btn_def_wl, 0, 0, 0)
        sizer_4.Add(sizer_5, 1, wx.EXPAND, 0)
        sizer_10.Add(self.cb_errorbar, 0, 0, 0)
        sizer_4.Add(sizer_10, 1, wx.EXPAND, 0)
        sizer_9.Add(self.btn_VIS2_vs_freq, 0, 0, 0)
        sizer_9.Add(self.btn_VIS_vs_freq, 0, 0, 0)
        sizer_9.Add(self.btn_vis2_tf2_cp_vs_time, 0, 0, 0)
        sizer_9.Add(self.btn_vis_tf_cp_vs_time, 0, 0, 0)
        sizer_4.Add(sizer_9, 1, wx.EXPAND, 0)
        grid_sizer_1.Add(sizer_4, 1, wx.EXPAND, 0)
        self.panel.SetSizer(grid_sizer_1)
        sizer_1.Add(self.panel, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

    ###########################################################################

    def OnButtonClicked(self, e):
        plot_t_flag = self.cb_pl_OI_t.GetValue()
        plot_wl_flag = self.cb_pl_OI_wl.GetValue()
        plot_errorbar_flag = self.cb_errorbar.GetValue()
        self.DIT = float(self.tb_DIT.GetValue())
        self.DIT_range = float(self.tb_DIT_range.GetValue())
        band_mask = [self.cb_b_L.GetValue(),self.cb_b_M.GetValue(),self.cb_b_LM.GetValue(),self.cb_b_N.GetValue()]
        selected_bands = self.bands[band_mask].tolist()
        spectral_resolution_mask = [self.cb_disp_LOW.GetValue(),self.cb_disp_MED.GetValue(),self.cb_disp_HIGH.GetValue()]
        selected_spectral_resolutions = self.spectral_resolutions[spectral_resolution_mask].tolist()
        self.wl_min = float(self.tb_wl_min.GetValue())
        self.wl_max = float(self.tb_wl_max.GetValue())
        self.target_selected = self.tb_target.GetValue()
        eID = e.GetId()
        if eID < 13 or eID > 15:
            #first check if list of dictionaries is empty
            if len(self.list_of_dicts) == 1:
                if not self.list_of_dicts[0]:
                    self.statusbar.SetStatusText('Load OIFITS data.')
                    print ('Load OIFITS data.')
                    self.name_file = self.tb_path.GetValue()
                    self.OI_load_data(self.name_file)

            if len(self.list_of_dicts) > 1:
                if not self.target_selected.strip():
                    self.filtered_list_of_dicts = filter_oi_list(self.list_of_dicts, dates=[self.date],
                                   bands=selected_bands,
                                   spectral_resolutions=selected_spectral_resolutions,
                                   DIT_range=[self.DIT - self.DIT_range,
                                              self.DIT + self.DIT_range],
                                   targets=[])
                else:
                    self.filtered_list_of_dicts = filter_oi_list(self.list_of_dicts, dates=[self.date], bands=selected_bands, spectral_resolutions=selected_spectral_resolutions, DIT_range=[self.DIT - self.DIT_range, self.DIT + self.DIT_range], targets=[self.target_selected])
            else:
                self.filtered_list_of_dicts = self.list_of_dicts

                # Display of visibility amplitude
            if eID == 0:
                self.statusbar.SetStatusText('Show VISAMP.')
                print ('Show VISAMP.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts, [self.wl_min,self.wl_max], key="VIS", datatype="CFLUX",showvis=False,plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS',
                                    datatype='CFLUX',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of squared visibility
            elif eID == 1:
                self.statusbar.SetStatusText('Show VIS2.')
                print ('Show VIS2.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts, [self.wl_min,self.wl_max], key="VIS2", datatype="VIS2",showvis=False,plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS2',
                                    datatype='VIS2',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of differential visibility
            elif eID == 2:
                self.statusbar.SetStatusText('Show VISPHI.')
                print ('Show VISPHI.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="VIS",
                                    datatype="DPHI",showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS',
                                    datatype='DPHI',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of transfer function
            elif eID == 3:
                self.statusbar.SetStatusText('Show TF2.')
                print ('Show TF2.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="TF2",
                                    datatype="TF2",showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='TF2',
                                    datatype='TF2',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of triple product amplitude
            elif eID == 4:
                self.statusbar.SetStatusText('Show T3AMP.')
                print ('Show T3AMP.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="T3",
                                    datatype="T3AMP",showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='T3',
                                    datatype='T3AMP',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of flux
            elif eID == 5:
                self.statusbar.SetStatusText('Show FLUX.')
                print ('Show FLUX.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="FLUX",
                                    datatype="FLUX",showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='FLUX',
                                    datatype='FLUX',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of closure phase
            elif eID == 6:
                self.statusbar.SetStatusText('Show T3PHI.')
                print ('Show T3PHI.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="T3",
                                    datatype="CLOS",showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],
                                    key='T3',datatype='CLOS',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of differential visibility
            elif eID == 7:
                self.statusbar.SetStatusText('Show VIS.')
                print ('Show VIS.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min,self.wl_max], key="VIS2",
                                    datatype="VIS2",showvis=True,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS2',
                                    datatype='VIS2',showvis=True,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of transfer function
            elif eID == 8:
                self.statusbar.SetStatusText('Show TF.')
                print ('Show TF.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts,
                                    [self.wl_min, self.wl_max], key="TF2",
                                    datatype="TF2",
                                    showvis=True,
                                    plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0], key='TF2',
                                    datatype='TF2', showvis=True,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of visibilities as a function of spatial frequency
            elif eID == 9:
                self.statusbar.SetStatusText('Show VIS2 vs. spatial frequency.')
                print ('Show VIS2 vs. spatial frequency.')
                show_oi_vs_freq(self.filtered_list_of_dicts[0],log=False,
                                showvis=False)
            elif eID == 10:
                self.statusbar.SetStatusText('Show VIS vs. spatial frequency.')
                print ('Show VIS vs. spatial frequency.')
                show_oi_vs_freq(self.filtered_list_of_dicts[0],log=False,
                                showvis=True)
            elif eID == 11:
                self.statusbar.SetStatusText('Show VIS2, TF2, T3PHI vs time.')
                print ('Show VIS2, TF2, T3PHI vs time.')
                show_vis2_tf2_vs_time(self.filtered_list_of_dicts,
                                      [self.wl_min,self.wl_max], showvis=False, saveplots=False, output_path="",plot_errorbars=plot_errorbar_flag)
            elif eID == 12:
                self.statusbar.SetStatusText('Show VIS, TF, T3PHI vs time.')
                print ('Show VIS2, TF2, T3PHI vs time.')
                show_vis2_tf2_vs_time(self.filtered_list_of_dicts,
                                      [self.wl_min,self.wl_max], showvis=True,
                                      saveplots=False, output_path="",
                                      plot_errorbars=plot_errorbar_flag)
                # Display of squared correlated flux
            elif eID == 17:
                self.statusbar.SetStatusText('Show CFX2.')
                print ('Show CFX2.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts, [self.wl_min,self.wl_max], key="VIS2", datatype="CFX2",showvis=False,plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS2',
                                    datatype='CFX2',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
                # Display of correlated flux
            elif eID == 16:
                self.statusbar.SetStatusText('Show CFXAMP.')
                print ('Show CFXAMP.')
                if plot_t_flag == True:
                    show_oi_vs_time(self.filtered_list_of_dicts, [self.wl_min,self.wl_max], key="VIS", datatype="CFLUX",showvis=False,plot_errorbars=plot_errorbar_flag)
                elif plot_wl_flag == True:
                    show_oi_vs_wlen(self.filtered_list_of_dicts[0],key='VIS',
                                    datatype='CFLUX',showvis=False,
                                    plot_errorbars=plot_errorbar_flag)
        else:
            if eID == 13:
                self.statusbar.SetStatusText('Default wavelength range.')
                print ('Default wavelength range.')
                self.tb_wl_min.Clear()
                self.tb_wl_max.Clear()
                if self.cb_b_J.GetValue() == True:
                    self.wl_min = self.wl_min_def[0] #um, J,H,K,L,M,N,LM
                    self.wl_max = self.wl_max_def[0]
                if self.cb_b_H.GetValue() == True:
                    self.wl_min = self.wl_min_def[1]  # um, J,H,K,L,M,N,LM
                    self.wl_max = self.wl_max_def[1]
                if self.cb_b_K.GetValue() == True:
                    self.wl_min = self.wl_min_def[2]  # um, J,H,K,L,M,N,LM
                    self.wl_max = self.wl_max_def[2]
                if self.cb_b_M.GetValue() == True:
                    self.wl_min = self.wl_min_def[4] #um, J,H,K,L,M,N,LM
                    self.wl_max = self.wl_max_def[4]
                if self.cb_b_L.GetValue() == True:
                    self.wl_min = self.wl_min_def[3]
                    self.wl_max = self.wl_max_def[3]
                if self.cb_b_N.GetValue() == True:
                    self.wl_min = self.wl_min_def[5]
                    self.wl_max = self.wl_max_def[5]
                if self.cb_b_LM.GetValue() == True:
                    self.wl_min = self.wl_min_def[6]
                    self.wl_max = self.wl_max_def[6]
                self.tb_wl_min.write("%.2f" % (self.wl_min))
                self.tb_wl_max.write("%.2f" % (self.wl_max))
                # e.Skip()
            elif eID == 15:
                self.statusbar.SetStatusText('Open OIFITS data.')
                print ('Open OIFITS data.')
                print("Running file selector...")
                openFileDialog = mat_FileDialog(None, 'Open a file', "lmk,")
                if openFileDialog.ShowModal() == wx.ID_OK:
                    self.name_file = openFileDialog.GetPaths()[0]
                    print(self.name_file)
                else:
                    self.name_file = ""
                openFileDialog.Destroy()
                self.tb_path.SetValue(self.name_file)
                self.OI_load_data(self.name_file)
            elif eID == 14:
                self.statusbar.SetStatusText('Load OIFITS data.')
                print ('Load OIFITS data.')
                print(self.tb_path.GetValue())
                if not self.tb_path.GetValue():
                    openFileDialog = mat_FileDialog(None, 'Open a file', "lmk,")
                    if openFileDialog.ShowModal() == wx.ID_OK:
                        self.name_file = openFileDialog.GetPaths()[0]
                        self.OI_load_data(self.name_file)
                else:
                    self.name_file = self.tb_path.GetValue()
                    self.OI_load_data(self.name_file)

    ###########################################################################

    def OI_load_data(self,path):
        if os.path.isfile(path):
            self.dic = {}
            print("Reading file " + path + "...")
            self.dic = open_oi(path)
            self.list_of_dicts = [self.dic]
            # update the values in the form
            if self.dic['BAND'] == 'LM':
                self.cb_b_L.SetValue(False)
                self.cb_b_M.SetValue(False)
                self.cb_b_N.SetValue(False)
                self.cb_b_LM.SetValue(True)
                self.wl_min = self.wl_min_def[6] # um, J,H,K,L,M,N,LM
                self.wl_max = self.wl_max_def[6]
            elif self.dic['BAND'] == 'N':
                self.cb_b_L.SetValue(False)
                self.cb_b_M.SetValue(False)
                self.cb_b_N.SetValue(True)
                self.cb_b_LM.SetValue(False)
                self.wl_min = self.wl_min_def[5]
                self.wl_max = self.wl_max_def[5]

            if self.dic['DISP'] == 'LOW':
                self.cb_disp_LOW.SetValue(True)
                self.cb_disp_MED.SetValue(False)
                self.cb_disp_HIGH.SetValue(False)
            elif self.dic['DISP'] == 'MED':
                self.cb_disp_LOW.SetValue(False)
                self.cb_disp_MED.SetValue(True)
                self.cb_disp_HIGH.SetValue(False)
            elif self.dic['DISP'] == 'HIGH':
                self.cb_disp_LOW.SetValue(False)
                self.cb_disp_MED.SetValue(False)
                self.cb_disp_HIGH.SetValue(True)

            self.tb_DIT.SetValue("%.3f" % (self.dic['DIT']))
            self.tb_wl_min.Clear()
            self.tb_wl_max.Clear()
            self.tb_wl_min.write("%.2f" % (self.wl_min))
            self.tb_wl_max.write("%.2f" % (self.wl_max))
        elif os.path.isdir(path):
            name_dir = path
            self.dic = {}
            self.list_of_dicts = [{}]
            self.list_of_dicts = open_oi_dir(name_dir)
        else:
            "File/directory not found."
        self.date = self.list_of_dicts[0]['DATEOBS'][0:10]
        self.tb_date.SetValue(self.date)

###############################################################################

class OI_plotter(wx.App):
    def OnInit(self):
        self.frame = oi_data_select_frame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True

###############################################################################

if __name__ == '__main__':
    listArg = sys.argv
    name_file = []
    typePlot = "VIS2"
    for elt in listArg:
        if ('--help' in elt):
            print("Usage: mat_show_....py [--dir=start directory]")
            sys.exit(0)
        elif ('--typePlot' in elt):
            typePlot = elt.split("=")[1]
            # elif len(listArg) == 2:
            #     name_file = sys.argv[1]
            #     print(name_file)

    if typePlot == "VIS2":
        table = "VIS2"
    elif typePlot == "CLOS":
        table = "T3"
    elif typePlot == "TF2":
        typePlot = "TF2"
        table = "TF2"
    elif typePlot == "FLUX":
        table = "FLUX"
    elif typePlot == "DPHI":
        table = "VIS"

    app = OI_plotter(0)
    app.MainLoop()

    # app = wx.App()
    # if not name_file:
    #     print("No input name given, running file selector...")
    #     openFileDialog = mat_FileDialog(None, 'Open a file', "lmk,")
    #     if openFileDialog.ShowModal() == wx.ID_OK:
    #         name_file = openFileDialog.GetPaths()[0]
    #         print(name_file)
    #     openFileDialog.Destroy()
    # app.MainLoop()
    # app.Destroy()
    #
    # dic = {};
    # if os.path.isfile(name_file):
    #     print("Reading file " + name_file + "...")
    #     dic = open_oi(name_file)
    #     print("Plotting data " + name_file + "...")
    #     show_oi_vs_freq(dic)
    #     print("Plotting data " + name_file + "...")
    #     plt.figure()
    #     show_oi_vs_wlen(dic['FLUX'], dic['WLEN'], datatype='FLUX')
    #     print("Plotting data " + name_file + "...")
    #     plt.figure()
    #     wlen = dic['WLEN']
    #     print(wlen)
    # elif os.path.isdir(name_file):
    #     name_dir = name_file
    #     list_of_dicts = open_oi_dir(name_dir)
    #     filtered_list_of_dicts = filter_oi_list(list_of_dicts)
    #     show_vis2_tf2_vs_time(filtered_list_of_dicts, [3.5, 3.95])
    #     # show_oi_vs_time(filtered_list_of_dicts ,[3.5,3.95],key=table, datatype=typePlot)
