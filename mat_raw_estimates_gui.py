#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
  $Id:  $
 
  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 
Created on Wed Apr  5 10:18:07 2017

  $Author: fmillour $
  $Date:  $
  $Revision:  $
  $Name:  $
"""

import wx
import os
import pwd
import errno

fileTypes = ["HOT_DARK", "CALIB_SRC_RAW", "BADPIX", "OBS_FLATFIELD",
             "NONLINEARITY", "SHIFT_MAP", "KAPPA_MATRIX"]

mainPath = ""
sofFile  = ""
runDir   = ""
row = 0;
col = 0;
guiTitle = "mat_raw_estimates";
waitforupdate=0
flagnosave=0
 
########################################################################
class displayGui(wx.Frame):
    #----------------------------------------------------------------------
    def __init__(self):
        global row, col, mainPath, guiTitle
        
        wx.Frame.__init__(self, None, wx.ID_ANY, guiTitle)
        panel = wx.Panel(self, wx.ID_ANY)
 
        sizer = wx.GridBagSizer(3, len(fileTypes))
                         
        # Define the main buttons with the file selection here
               # input files part 
        self.btnDark = tributton(panel, sizer, fileTypes[0], 1)
        self.btnRaw  = tributton(panel, sizer, fileTypes[1], 1)
        self.btnFlat = tributton(panel, sizer, fileTypes[2], 1)
        self.btnShift= tributton(panel, sizer, fileTypes[3], 1)
        self.btnBPM  = tributton(panel, sizer, fileTypes[4], 1)
        self.btnNLM  = tributton(panel, sizer, fileTypes[5], 1)
        self.btnKappa  = tributton(panel, sizer, fileTypes[6], 0)
        self.btns = [self.btnDark, self.btnRaw, self.btnFlat, self.btnShift, self.btnBPM, self.btnNLM, self.btnKappa]
        
        # Buttons for interactivity are set here
                
               # SOF file part 
        col = 0
        row += 1
        btnOpenSOF = wx.Button(panel, label='Open SOF')
        sizer.Add(btnOpenSOF, (row, col))
        btnOpenSOF.Bind(wx.EVT_BUTTON, self.OnOpenSOF)
        
        col = 1
        self.btnDescSOF  = wx.StaticText(panel, label="SOF file", name="dirRun")
        sizer.Add(self.btnDescSOF, (row, col))
                
        col = 2
        self.btnTextSOF = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)
        sizer.Add(self.btnTextSOF, (row, col))
        self.btnTextSOF.Bind(wx.EVT_TEXT, self.OnTextSOF)
                
        col = 4
        btnSOF     = wx.Button(panel, label='Generate SOF')
        sizer.Add(btnSOF, (row, col))
        btnSOF.Bind(wx.EVT_BUTTON, self.OnGenSOF)
                
               # Run part
        col = 0
        row += 1
        btnMagic = wx.Button(panel, label='auto set')
        sizer.Add(btnMagic, (row, col))
        btnMagic.Bind(wx.EVT_BUTTON, self.OnOpenMagic)
        
        col = 1
        self.btnDescRun  = wx.StaticText(panel, label="output dir", name="dirRun")
        sizer.Add(self.btnDescRun, (row, col))
                
        col = 2
        self.btnTextRun = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)
        sizer.Add(self.btnTextRun, (row, col))
        self.btnTextRun.Bind(wx.EVT_TEXT, self.OnTextRun)
                
        col = 4
        btnRun = wx.Button(panel, label='Run!')
        sizer.Add(btnRun, (row, col))
        btnRun.Bind(wx.EVT_BUTTON, self.OnRunRex)
                
        col = 0
        row += 1
        btnOK = wx.Button(panel, label='Cancel', size=(60, -1))
        sizer.Add(btnOK, (row, col))
        btnOK.Bind(wx.EVT_BUTTON, self.OnClose)
         
        # Set simple sizer for a nice border
        border = wx.BoxSizer()
        border.Add(sizer, 5, wx.ALL | wx.EXPAND, 5)

        # Use the sizers
        panel.SetSizerAndFit(border)  
        self.SetSizerAndFit(sizer)
        
       # panel.SetSizer(sizer)
        
    def OnTextSOF(self,e):
        global sofFile
       # print("yeah man!")     
        sofFile  = self.btnTextSOF.GetValue()   
        
    def OnTextRun(self,e):
        global runDir
       # print("yeah man!")     
        runDir  = self.btnTextRun.GetValue()   
        
    def OnOpenMagic(self,e):
        global sofFile, runDir
        print("Hocus Pocus, abracadabra et tutti quanti !")
            
            # Test default directory
        for idx, btn in enumerate(self.btns):
            if btn.filelist:
                path = os.path.dirname(btn.filelist[0])
                if btn.typeFile == fileTypes[0] or btn.typeFile == fileTypes[1]:
                    runDir = path+"../Results/";
                    
            # Otherwise use SOF file directory
        if sofFile != "":
            runDir = os.path.dirname(sofFile)+"/../Results/";
            
        print("Set the output path to"+runDir)
        # Expand environment variables
        runDir = os.path.expandvars(runDir)
        self.btnTextRun.SetValue(runDir)
        
    def OnOpenSOF(self,e):
        global mainPath, sofFile, waitforupdate
        waitforupdate = 1;
        dlg = wx.FileDialog(
        None,
        "Choose a SOF file",
        mainPath,
        "",
        "SOF files (*.sof*)|*.sof*|" \
         "All files (*.*)|*.*"
        )
        if dlg.ShowModal() == wx.ID_OK:
            sofFile = dlg.GetPath()
            print "You chosed the following file(s):"
            print sofFile
            
        self.btnTextSOF.SetValue(sofFile)
        
        for idx, btn in enumerate(self.btns):
            btn.filelist = [];
        
        f = open(sofFile, 'r')
        for line in f:
            line = line.strip()
            if line:
                columns  = line.split()
                filename = columns[0]
                filetype = columns[1]
                for idx, btn in enumerate(self.btns):
                    if filetype == fileTypes[idx]:
                        btn.filelist.append(filename)
           # print(filename, filetype)
        f.close()
        
        for idx, btn in enumerate(self.btns):
            if btn.filelist:
                path = os.path.dirname(btn.filelist[0])
                if btn.typeFile == fileTypes[0] or btn.typeFile == fileTypes[1]:
                    mainPath = path
                    print("Set the main path to"+mainPath)
                    # Expand environment variables
                    mainPath = os.path.expandvars(mainPath)
                btn.dirtxt.SetValue(path)
                # Get file list
                mystring = [os.path.basename(fil) for fil in btn.filelist]
                text = mystring[0];
                for texti in mystring[1:]:
                    text += ", "+texti
               #     print("What's in the file")
               #     print(text)
                btn.filetxt.SetValue(text)
            
        waitforupdate = 0;
       # wx.MessageBox('Nothing there! Come again later ;-)', 'Error', wx.OK | wx.ICON_ERROR)
        
    def OnGenSOF(self,e):
        global mainPath, sofFile, guiTitle, flagnosave
        if not flagnosave:
                    
                # First things first: check values have been filled
            for idx, btn in enumerate(self.btns):
                if not btn.filelist:
                    if btn.checkPresent:
                        wx.MessageBox('Please set first '+fileTypes[idx], 'Error', wx.OK | wx.ICON_ERROR)
                        return -1
                
            print("writing SOF file "+sofFile)
        
            if sofFile == "":
                sofFile = mainPath+'/../'+pwd.getpwuid( os.getuid() )[ 0 ]+'/sof/'+ guiTitle+'.sof'
            self.make_sure_path_exists(sofFile)
        
            # Get file list
            text = "";
            for idx, btn in enumerate(self.btns):
                for texti in btn.filelist:
                    text += texti+" "+fileTypes[idx]+"\n"
            f = open(sofFile, 'w')
            f.writelines(text)
            f.close()
        else:
            wx.MessageBox('Please correct first red boxes', 'Error', wx.OK | wx.ICON_ERROR)
                            
    def OnRunRex(self,e):
        global guiTitle
        global runDir
        if runDir!="":
            self.make_sure_path_exists(runDir)
            os.chdir(runDir)
        else:
            wx.MessageBox('Please set an output path first', 'Error', wx.OK | wx.ICON_ERROR)
            return -1
            
        command = "esorex "+guiTitle+" "+sofFile+"&";
        print(command)
        os.system(command)
        print("Ouaf ! Ouaf !")
        wx.MessageBox("Process launched! This may take a while though. Why not drink a cup of tea?", 'Important information', wx.ICON_EXCLAMATION)
        #dlg.Destroy()
                        
    def OnClose(self,e):
        self.Close() 
                
    def make_sure_path_exists(self,path):
        directory = os.path.dirname(path)
        try:
            os.makedirs(directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
  
########################################################################
class tributton(wx.Frame):
    #-------------------
    def __init__(self, panel, sizer, typeFile, checkPresent):
        # Placing variables for GUI are global
        global row, col, mainPath
        self.checkPresent = checkPresent
        # output filelist for the recipe is init here
        self.filelist = [];
        # Put typefile as internal global variable for the class
        self.typeFile = typeFile
        
        # Select main directory
        self.txt     = wx.StaticText(panel, label=typeFile, name=typeFile)
        
        self.dirdsc  = wx.StaticText(panel, label="dir", name="dir")
        self.dirtxt  = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)
        self.dirmod  = wx.StaticText(panel, label="", name="dirmod")
        
        self.filedsc = wx.StaticText(panel, label="files", name="files")
        self.filetxt = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)
        self.filemod = wx.StaticText(panel, label="", name="filemod")
        
        self.btn     = wx.Button    (panel, label="Select", name=typeFile)
        
        self.dirtxt.Bind( wx.EVT_TEXT,   self.onDirTxt)
        self.filetxt.Bind(wx.EVT_TEXT,   self.onFileTxt)
        self.btn.Bind(    wx.EVT_BUTTON, self.onButton)
        #sizer.Add(btn, 0, wx.ALL, 5)
        
        row += 1
        col  = 0
        sizer.Add(self.txt, (row, col))
        col += 1
        sizer.Add(self.dirdsc, (row, col))
        col += 1
        sizer.Add(self.dirtxt, (row, col))
        col += 1
        sizer.Add(self.dirmod, (row, col))
        row += 1
        col  = 1
        sizer.Add(self.filedsc, (row, col))
        col += 1
        sizer.Add(self.filetxt,(row, col))
        col += 1
        sizer.Add(self.filemod, (row, col))
        col += 1
        sizer.Add(self.btn, (row, col))
                        
    #----------------------------------------------------------------------
    def onDirTxt(self, event):
        global mainPath, sofFile, waitforupdate, flagnosave
        if not waitforupdate:
            self.dirmod.SetLabel("*")
        
            text = self.dirtxt.GetValue()
            if not os.path.isdir(text):
                self.dirtxt.SetBackgroundColour(wx.RED)
            if self.typeFile == fileTypes[1]:
                mainPath = os.path.dirname(self.filelist[0])
                print("Set the main path to"+mainPath)
                flagnosave = 1;
            else:
                self.dirtxt.SetBackgroundColour(wx.WHITE)
                self.filelist = self.getFilelistFromText()
                flagnosave = 0
            
                    
    #----------------------------------------------------------------------
    def onFileTxt(self, event):
        global mainPath, sofFile, waitforupdate, flagnosave
        if not waitforupdate:
            self.filemod.SetLabel("*")
        
            filelist = self.getFilelistFromText()
           # print(filelist)
       #     print(txtfile)
            if not all(os.path.isfile(x) for x in filelist):
                self.filetxt.SetBackgroundColour(wx.RED)
                flagnosave = 1;
            else:
                self.filetxt.SetBackgroundColour(wx.WHITE)
                self.filelist = self.getFilelistFromText()
                flagnosave = 0
                            
    def getFilelistFromText(self):
        dirtxt  = self.dirtxt.GetValue()
        filtxt  = self.filetxt.GetValue()
        files = filtxt.split(", ")
        filelist = []
        for fle in files:
            filelist.append(dirtxt + "/" + fle)
        return filelist
        
        
    #----------------------------------------------------------------------
    def onButton(self, event):
        global mainPath, sofFile, waitforupdate
        """
        This method is fired when its corresponding button is pressed
        """
        waitforupdate = 1;
        #button = event.GetEventObject()
        
        dlg = wx.FileDialog(
        None,
        "Choose a "+self.typeFile,
        mainPath,
        "",
        "FITS files (*.fits*; *.oifits*)|*.fits*;*.oifits*|" \
         "FITS files (*.fits*)|*.fits*|" \
         "OIFITS files (*.oifits*)|*.oifits*|" \
         "All files (*.*)|*.*",
        wx.MULTIPLE
        )
        if dlg.ShowModal() == wx.ID_OK:
            self.filelist = dlg.GetPaths()
            
            print "You chosed the following file(s):"
            print self.filelist
        
            path = os.path.dirname(self.filelist[0])
            if self.typeFile == fileTypes[0] or self.typeFile == fileTypes[1]:
                mainPath = path
                print("Set the main path to"+mainPath)
                
            self.dirtxt.SetValue(path)
            # Get file list
            mystring = [os.path.basename(fil) for fil in self.filelist]
            text = mystring[0];
            for texti in mystring[1:]:
                text += ", "+texti
            #self.filetxt.SetValue(os.path.basename(text))
            self.filetxt.SetValue(text)
            
            print("selected files")
            print(text)
            print "You selected a "+self.typeFile
            dlg.Destroy()
            
        self.filemod.SetLabel("")
        self.filetxt.SetBackgroundColour(wx.WHITE)
        self.dirmod.SetLabel("")
        self.dirtxt.SetBackgroundColour(wx.WHITE)
        waitforupdate = 0;
                            
#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    app = wx.App(False)
    frame = displayGui()
    frame.Show()
    app.MainLoop()
