# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 13:50:15 2018

@author: fmillour
"""

from   astropy.io import fits as fits
from mat_show_oifits import open_oi,open_oi_dir,filter_oi_list
import numpy as np
import wx
import os
import time
from mat_fileDialog import mat_FileDialog


class changeKeys(wx.Frame):
    def __init__(self, parent, title):
        super(changeKeys, self).__init__(parent, title = title, size = (400,350))
        panel = wx.Panel(self)

        # Structure the window vertically
        vbox  = wx.BoxSizer(wx.VERTICAL)
        #vbox.AddSpacer(10)

        # Exit button
        self.exitButton    = wx.Button(panel,label="Exit")
        vbox.Add(self.exitButton,   border=10,flag=wx.ALIGN_CENTER)
        self.exitButton.Bind(wx.EVT_BUTTON, self.exitClicked)

        # Load data button
        self.loadDirButton = wx.Button(panel,label="Load files")
        vbox.Add(self.loadDirButton,border=10,flag=wx.EXPAND|wx.ALIGN_CENTER)
        self.loadDirButton.Bind(wx.EVT_BUTTON, self.loadFiles)
        panel.SetSizer(vbox)

        self.dir=os.getcwd()
        self.Centre()
        self.Show()
        self.Fit()

    def exitClicked(self,event):
        print("clicked itou!")
        self.Destroy()

    def loadFiles(self,event):
        print("clicked!")
        openFileDialog = mat_FileDialog(None, 'Select a directory',self.dir)
        print("marche ?")
        if openFileDialog.ShowModal() == wx.ID_OK:
            print(openFileDialog.GetPaths())
        print("marche pas ?")
        # Get list of files, target names, and data types
        self.dic = open_oi_dir(dir)
        print("ou bien ")

    def getUniques(dic):


        listofstuff = []
        for idat in datas:
            obj    = idat['TARGET']
            typ    = idat['HDR']['ESO PRO CATG']
            caract = obj+';'+typ
            print(obj+';'+typ)
            listofstuff = np.append(listofstuff, caract)

        redlist = set(listofstuff)

    def changeOiFitsToCalib(filename):
        # Change value of
        fits.setval(filename, 'ESO PRO CATG', value='CALIB_RAW_INT')

    def changeOiFitsToScience(filename):
        # Change value of
        fits.setval(filename, 'ESO PRO CATG', value='CALIB_RAW_INT')



