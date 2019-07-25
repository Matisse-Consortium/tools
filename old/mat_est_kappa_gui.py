#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
  Created on Wed Apr 12 12:04:01 2017

@author: ame
"""

import wx
from mat_generic_gui import displayGui


fileTypes = ["KAPPA_HOTDARK", "KAPPA_SRC", "BADPIX", "OBS_FLATFIELD",
             "NONLINEARITY", "SHIFT_MAP"]
checkPresent = [1,1,1,1,1,1]
             
GuiTitle="mat_est_kappa"

class mat_est_kappa_gui(displayGui):
    def __init__(self):
          super(mat_est_kappa_gui, self).__init__(GuiTitle,".",fileTypes, checkPresent)
       


if __name__ == "__main__":
    app = wx.App(False)
    frame = mat_est_kappa_gui()
    frame.Show()
    app.MainLoop()
    app.Destroy()
