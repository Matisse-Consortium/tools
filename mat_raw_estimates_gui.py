# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 12:04:01 2017

@author: ame
"""

import wx
from mat_generic_gui import displayGui


fileTypes = ["HOT_DARK", "CALIB_SRC_RAW", "BADPIX", "OBS_FLATFIELD",
             "NONLINEARITY", "SHIFT_MAP", "KAPPA_MATRIX"]
checkPresent = [1,1,1,1,1,1,0]
             
GuiTitle="MAT_RAW_ESTIMATE"

class mat_raw_estimates_gui(displayGui):
    def __init__(self):
          super(mat_raw_estimates_gui, self).__init__(GuiTitle,".",fileTypes, checkPresent)
       


if __name__ == "__main__":
    app = wx.App(False)
    frame = mat_raw_estimates_gui()
    frame.Show()
    app.MainLoop()
    app.Destroy()

 