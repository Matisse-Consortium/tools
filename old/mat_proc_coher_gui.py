#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:18:07 2017

@author: fmillour
"""
import wx
import os
from mat_generic_gui import displayGui


fileTypes = ["OBJ_CORR_FLUX", "PHOT_BEAMS", "OI_OPDWVPO"]
checkPresent = [1,1,0]
             
GuiTitle="mat_proc_coher"

class mat_proc_coher_gui(displayGui):
    def __init__(self):
          super(mat_proc_coher_gui, self).__init__(GuiTitle,".",fileTypes, checkPresent)
       


if __name__ == "__main__":
    app = wx.App(False)
    frame = mat_proc_coher_gui()
    frame.Show()
    app.MainLoop()
    app.Destroy()
