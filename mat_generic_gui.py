#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
  $Id: mat_raw_estimates_gui.py 61 2017-04-12 08:09:56Z fmillour $
 
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
  $Date: 2017-04-12 10:09:56 +0200 (mer., 12 avr. 2017) $
  $Revision: 61 $
  $Name:  $
"""

import wx
import os
import pwd
import errno
from mat_fileDialog import mat_FileDialog 

fileTypes = ["HOT_DARK", "CALIB_SRC_RAW", "BADPIX", "OBS_FLATFIELD","NONLINEARITY", "SHIFT_MAP", "KAPPA_MATRIX"]
checkPresent = [1,1,1,1,1,1,0]
             

                

########################################################################
class displayGui(wx.Frame):
    #----------------------------------------------------------------------
    def __init__(self,guiTitle="SOF GUI",mainpath=".",fileTypes=[],checkPresent=[]):
        
        self.waitforupdate=0    
        self.flagnosave=0
        self.sofFile= ""
        self.runDir= ""
        self.guiTitle=guiTitle
        self.mainPath=mainpath
        self.fileTypes=fileTypes
   
        
        wx.Frame.__init__(self, None, wx.ID_ANY, guiTitle)
        panel = wx.Panel(self, wx.ID_ANY)
             
        vbox=wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
                         
        # Define the main buttons with the file selection here
        # input files part 
        for typi in self.fileTypes:
            self.__dict__[typi] = tributton(self,panel, typi)
            vbox.Add(self.__dict__[typi],border=20,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        self.btns = [self.__dict__[typi] for typi in fileTypes]
        vbox.AddSpacer(30)     
        
        grid=wx.GridBagSizer(2,4)        
        btnOpenSOF       = wx.Button    (panel, label='Open SOF',size=(100, -1))
        self.btnDescSOF  = wx.StaticText(panel, label="SOF file", name="dirRun",size=(80, -1))   
        self.btnTextSOF  = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)
        btnSOF           = wx.Button    (panel, label='Generate SOF',size=(110, -1))
        btnMagic         = wx.Button    (panel, label='auto set',size=(100, -1))
        self.btnDescRun  = wx.StaticText(panel, label="output dir", name="dirRun",size=(80, -1))               
        self.btnTextRun  = wx.TextCtrl  (panel, size=(400, -1),style=wx.TE_RICH)            
        btnRun           = wx.Button    (panel, label='Run!',size=(110, -1))               
        grid.Add(btnOpenSOF, (0, 0))    
        grid.Add(self.btnDescSOF, (0, 1))
        grid.Add(self.btnTextSOF, (0, 2))     
        grid.Add(btnSOF, (0, 4)) 
        grid.Add(btnMagic, (1, 0))
        grid.Add(self.btnDescRun, (1, 1))    
        grid.Add(self.btnTextRun, (1, 2))        
        grid.Add(btnRun, (1, 4))           
        vbox.Add(grid,border=10,flag=wx.LEFT|wx.RIGHT|wx.EXPAND) 
        
    
        vbox.AddSpacer(30) 
        
        btnCancel = wx.Button(panel,wx.ID_CANCEL, label='Cancel', size=(60, -1))
        vbox.Add(btnCancel,flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)

        vbox.AddSpacer(10)       
       
        panel.SetSizerAndFit(vbox)  
        self.SetSizerAndFit(vbox)     
        
        btnOpenSOF.Bind(wx.EVT_BUTTON, self.OnOpenSOF)
        self.btnTextSOF.Bind(wx.EVT_TEXT, self.OnTextSOF)
        btnSOF.Bind(wx.EVT_BUTTON, self.OnGenSOF)
        btnMagic.Bind(wx.EVT_BUTTON, self.OnOpenMagic)
        btnRun.Bind(wx.EVT_BUTTON, self.OnRunRex)
        self.btnTextRun.Bind(wx.EVT_TEXT, self.OnTextRun)
        btnCancel.Bind(wx.EVT_BUTTON, self.OnClose)
        
    def OnTextSOF(self,e):
        print("yeah man!")     
        self.sofFile  = self.btnTextSOF.GetValue()   
        
    def OnTextRun(self,e):
        print("yeah man!")     
        self.runDir  = self.btnTextRun.GetValue()   
        
    def OnOpenMagic(self,e):
        print("Hocus Pocus, abracadabra et tutti quanti !")
            
            # Test default directory
        for idx, btn in enumerate(self.btns):
            if btn.filelist:
                path = os.path.dirname(btn.filelist[0])
                if btn.typeFile == fileTypes[0] or btn.typeFile == self.fileTypes[1]:
                    self.runDir = path+"../Results/";
                    
            # Otherwise use SOF file directory
        if self.sofFile != "":
            self.runDir = os.path.dirname(self.sofFile)+"/../Results_"+os.path.basename(os.path.splitext(self.sofFile)[0])+"/";
            
        print("Set the output path to"+self.runDir)
        # Expand environment variables
        self.runDir = os.path.expandvars(self.runDir)
        self.btnTextRun.SetValue(self.runDir)
        
    def OnOpenSOF(self,e):
        #global mainPath, sofFile, waitforupdate
        self.waitforupdate = 1;
        
        #Init all text boxes
        self.btnTextSOF.SetValue("")
        self.btnTextRun.SetValue("")
            # Test default directory
        for btn in self.btns:
            btn.dirtxt.SetValue("")
            btn.filetxt.SetValue("")
        dlg = wx.FileDialog(
        None,
        "Choose a SOF file",
        self.mainPath,
        "",
        "SOF files (*.sof*)|*.sof*|" \
         "All files (*.*)|*.*"
        )
        if dlg.ShowModal() == wx.ID_OK:
            self.sofFile = dlg.GetPath()
            print "You chosed the following file(s):"
            print self.sofFile
            
        self.btnTextSOF.SetValue(self.sofFile)
        
        for idx, btn in enumerate(self.btns):
            btn.filelist = [];
        
        f = open(self.sofFile, 'r')
        for line in f:
            line = line.strip()
            if line:
                columns  = line.split()
                filename = columns[0]
                filetype = columns[1]
                for idx, btn in enumerate(self.btns):
                    if filetype == self.fileTypes[idx]:
                        btn.filelist.append(filename)
           # print(filename, filetype)
        f.close()
        
        for idx, btn in enumerate(self.btns):
            if btn.filelist:
                path = os.path.dirname(btn.filelist[0])
                if btn.typeFile == self.fileTypes[0] or btn.typeFile == self.fileTypes[1]:
                    self.mainPath = path
                    print("Set the main path to"+self.mainPath)
                    # Expand environment variables
                    self.mainPath = os.path.expandvars(self.mainPath)
                btn.dirtxt.SetValue(path)
                # Get file list
                mystring = [os.path.basename(fil) for fil in btn.filelist]
                text = mystring[0];
                for texti in mystring[1:]:
                    text += ", "+texti
               #     print("What's in the file")
               #     print(text)
                btn.filetxt.SetValue(text)
            
        self.waitforupdate = 0;
       # wx.MessageBox('Nothing there! Come again later ;-)', 'Error', wx.OK | wx.ICON_ERROR)
        
    def OnGenSOF(self,e):
        #global mainPath, sofFile, guiTitle, flagnosave
        if not self.flagnosave:
                    
                # First things first: check values have been filled
            for idx, btn in enumerate(self.btns):
                if not btn.filelist:
                    if btn.checkPresent:
                        wx.MessageBox('Please set first '+self.fileTypes[idx], 'Error', wx.OK | wx.ICON_ERROR)
                        return -1
                
        
            if self.sofFile == "":
                self.sofFile = self.mainPath+'/../'+pwd.getpwuid( os.getuid() )[ 0 ]+'/sof/'+ self.guiTitle+'.sof' 
                #self.sofFile = self.mainPath+'/../'+'/sof/'+ self.guiTitle+'.sof'  => For Windows!!!

            print("writing SOF file "+self.sofFile)
            self.make_sure_path_exists(self.sofFile)
        
            # Get file list
            text = "";
            for idx, btn in enumerate(self.btns):
                for texti in btn.filelist:
                    text += texti+" "+self.fileTypes[idx]+"\n"
            f = open(self.sofFile, 'w')
            f.writelines(text)
            f.close()
        else:
            wx.MessageBox('Please correct first red boxes', 'Error', wx.OK | wx.ICON_ERROR)
                            
    def OnRunRex(self,e):
        global guiTitle
        global runDir
        if self.runDir!="":
            self.make_sure_path_exists(self.runDir)
            os.chdir(self.runDir)
        else:
            wx.MessageBox('Please set an output path first', 'Error', wx.OK | wx.ICON_ERROR)
            return -1
            
        command = "esorex "+self.guiTitle+" "+self.sofFile+"&";
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
class tributton(wx.GridBagSizer):
    #-------------------
    def __init__(self,parent,panel, typeFile, checkPresent=1):
        super(tributton, self).__init__(2,4)
        self.parent=parent
        # Placing variables for GUI are global
       
        self.checkPresent = checkPresent
        # output filelist for the recipe is init here
        self.filelist = [];
        # Put typefile as internal global variable for the class
        self.typeFile = typeFile
        
        # Select main directory
        self.txt     = wx.StaticText (panel, label=typeFile, name=typeFile,size=(120, -1))       
        self.dirdsc  = wx.StaticText (panel, label="dir", name="dir",size=(40, -1))
        self.dirtxt  = wx.TextCtrl   (panel, size=(400, -1),style=wx.TE_RICH)
        self.dirmod  = wx.StaticText (panel, label="", name="dirmod",size=(7, -1))        
        self.filedsc = wx.StaticText (panel, label="files", name="files",size=(40, -1))
        self.filetxt = wx.TextCtrl   (panel, size=(400, -1),style=wx.TE_RICH)
        self.filemod = wx.StaticText (panel, label="", name="filemod",size=(7, -1))      
        self.btn     = wx.Button     (panel, label="Select", name=typeFile)        
        self.dirtxt.Bind( wx.EVT_TEXT,   self.onDirTxt)
        self.filetxt.Bind(wx.EVT_TEXT,   self.onFileTxt)
        self.btn.Bind(    wx.EVT_BUTTON, self.onButton)
        #sizer.Add(btn, 0, wx.ALL, 5)     
         
        self.Add(self.txt,     (0, 0))
        self.Add(self.dirdsc,  (0, 1))
        self.Add(self.dirtxt,  (0, 2))
        self.Add(self.dirmod,  (0, 3))
        self.Add(self.filedsc, (1, 1))
        self.Add(self.filetxt, (1, 2))
        self.Add(self.filemod, (1, 3))
        self.Add(self.btn,     (1, 4))
                        
    #----------------------------------------------------------------------
    def onDirTxt(self, event):
        if not self.parent.waitforupdate:
            self.dirmod.SetLabel("*")
        
            text = self.dirtxt.GetValue()
            if not os.path.isdir(text):
                self.dirtxt.SetBackgroundColour(wx.RED)
            if self.typeFile == self.parent.fileTypes[1]:
                self.parent.mainPath = os.path.dirname(self.filelist[0])
                print("Set the main path to"+self.parent.mainPath)
                self.parent.flagnosave = 1;
            else:
                self.dirtxt.SetBackgroundColour(wx.WHITE)
                self.filelist = self.getFilelistFromText()
                self.parent.flagnosave = 0
            
                    
    #----------------------------------------------------------------------
    def onFileTxt(self, event):
        if not self.parent.waitforupdate:
            self.filemod.SetLabel("*")
        
            filelist = self.getFilelistFromText()
           # print(filelist)
       #     print(txtfile)
            if not all(os.path.isfile(x) for x in filelist):
                self.filetxt.SetBackgroundColour(wx.RED)
                self.parent.flagnosave = 1;
            else:
                self.filetxt.SetBackgroundColour(wx.WHITE)
                self.filelist = self.getFilelistFromText()
                self.parent.flagnosave = 0
                            
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
        """
        This method is fired when its corresponding button is pressed
        """
        self.parent.waitforupdate = 1;
        
        dlg=mat_FileDialog(None,"Choose a "+self.typeFile, self.parent.mainPath)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.filelist = dlg.GetPaths()
            
            print "You chosed the following file(s):"
            print self.filelist
        
            path = os.path.dirname(self.filelist[0])
            if self.typeFile == self.parent.fileTypes[0] or self.typeFile == self.parent.fileTypes[1]:
                self.parent.mainPath = path
                print("Set the main path to"+self.parent.mainPath)
                
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
        self.parent.waitforupdate = 0;                
                
                
#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    app = wx.App(False)
    frame = displayGui("Test","./",fileTypes, checkPresent)
    frame.Show()
    app.MainLoop()
    app.Destroy()