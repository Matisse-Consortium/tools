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
  2018-03-26: new GUI interface ready: oi_data_select_frame (jvarga)
  2018-04-04: updated GUI and extended functionality: input file/folder textbox, filter for target name,
              more bands (JHK) available (for e.g. AMBER data), plot with or without errorbars, plot V or V2 (jvarga)

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
import robust


###############################################################################
# showvis: if True, plot visibilities (V) instead of V^2: V is calculated from V^2 (not from the VISAMP table)
def show_seeing(list_of_dicts, saveplots=False, output_path="",show=True):
    # check if list is not empty:
    if list_of_dicts:
        seeing_arr = []
        tau0_arr = []
        mjd_arr =[]
 
        for dic in list_of_dicts:
            seeing_arr.append(dic['SEEING'])
            tau0_arr.append(dic['TAU0']*1.E3)
            mjd_arr.append(np.array(dic['VIS2']["TIME"])[0])
        fig=plt.figure(figsize=(7, 7))
        axs2 = fig.add_subplot(2, 1, 2)
        axs1 = fig.add_subplot(2, 1, 1, sharex=axs2 )
        plt.setp(axs1.get_xaxis().get_offset_text(), visible=False)
        axs1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        
        axs2.plot(mjd_arr,tau0_arr,'s')
        axs2.set_ylabel("Tau0 (ms)")
        axs2.set_xlabel('$\mathrm{MJD}$')
      
        axs1.plot(mjd_arr,seeing_arr,'d')
        axs1.set_ylabel("Seeing")
 
        axs1.set_ylim(0.,np.max(seeing_arr)+0.5)
        axs2.set_ylim(0.,np.max(tau0_arr)+1.)
        
        for dic in list_of_dicts:
            target_names=dic['TARGET']
            axs1.text(np.array(dic['VIS2']["TIME"])[0], np.max(seeing_arr)+0.5+0.1,target_names.replace('_', ' '), rotation=90, va='bottom',fontsize=8)
 
        
        plt.subplots_adjust(hspace=0)

        if saveplots == True:
            label = '_SEEING'
            fig.savefig(output_path + label + '.png', dpi=150)
            fig.savefig(output_path + label + '.eps', format='eps', dpi=300)
            plt.close(fig)
        if show == True:
            plt.show()

    return 0
