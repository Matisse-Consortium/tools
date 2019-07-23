#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
  $Id$

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur

  Created on Wed Apr 12 12:04:01 2017
  @author: ame

  This software offers a GUI to the esorex MATISSE recipes.

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

import wx
from mat_generic_gui import displayGui
from mat_fileDialog import identifyFile


fileTypes = ["HOT_DARK", "CALIB_SRC_RAW", "BADPIX", "OBS_FLATFIELD",
             "NONLINEARITY", "SHIFT_MAP", "KAPPA_MATRIX"]
checkPresent = [1,1,1,1,1,1,0]
             
GuiTitle="mat_raw_estimates"

class mat_raw_estimates_gui(displayGui):
    def __init__(self):
          super(mat_raw_estimates_gui, self).__init__(GuiTitle,".",fileTypes, checkPresent)
       
if __name__ == "__main__":
    app   = wx.App(False)
    frame = mat_raw_estimates_gui()
    frame.Show()
    app.MainLoop()
    app.Destroy()
