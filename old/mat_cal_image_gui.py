#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:18:07 2017

@author: fmillour
"""

import wx
from mat_generic_gui import displayGui

fileTypes = ["HOT_DARK", "CALIB_SRC_RAW", "BADPIX", "OBS_FLATFIELD",
             "SHIFT_MAP", "NONLINEARITY"]
             
checkPresent = [1,1,1,1,1,1]
             
GuiTitle="mat_est_kappa"

class mat_cal_image_gui(displayGui):
    def __init__(self):
          super(mat_cal_image_gui, self).__init__(GuiTitle,".",fileTypes, checkPresent)
       


if __name__ == "__main__":
    app = wx.App(False)
    frame = mat_cal_image_gui()
    frame.Show()
    app.MainLoop()
    app.Destroy()
