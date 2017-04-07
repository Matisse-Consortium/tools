#!/usr/bin/python  
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:18:07 2017

@author: fmillour
"""
import wx
import os

fileTypes = ["OBJ_CORR_FLUX", "PHOT_BEAMS", "OI_OPDWVPO"]

mainPath = ""
sofFile  = ""
row = 0;
col = 0;
guiTitle = "mat_proc_coher";
 
########################################################################
class displayGui(wx.Frame):
    #----------------------------------------------------------------------
    def __init__(self):
        global row, col, mainPath, guiTitle
        
        wx.Frame.__init__(self, None, wx.ID_ANY, guiTitle)
        panel = wx.Panel(self, wx.ID_ANY)
 
        sizer = wx.GridBagSizer(3, len(fileTypes))
                         
        # Define the main buttons with the file selection here
        self.btnCFlux = tributton(panel, sizer, fileTypes[0], 1)
        self.btnPBeams= tributton(panel, sizer, fileTypes[1], 1)
        self.btnOPD   = tributton(panel, sizer, fileTypes[2], 0)
        self.btns = [self.btnCFlux, self.btnPBeams, self.btnOPD]
        
        # Buttons for interactivity are set here
        btnOK = wx.Button(panel, label='Cancel', size=(60, -1))
        col  = 0
        row += 1
        sizer.Add(btnOK, (row, col))
        btnOK.Bind(wx.EVT_BUTTON, self.OnClose)
                
        btnOpenSOF = wx.Button(panel, label='Open SOF')
        col  = 1
        sizer.Add(btnOpenSOF, (row, col))
        btnOpenSOF.Bind(wx.EVT_BUTTON, self.OnOpenSOF)
                
        btnSOF = wx.Button(panel, label='Generate SOF')
        col  = 2
        sizer.Add(btnSOF, (row, col))
        btnSOF.Bind(wx.EVT_BUTTON, self.OnGenSOF)
         
        # Set simple sizer for a nice border
        border = wx.BoxSizer()
        border.Add(sizer, 5, wx.ALL | wx.EXPAND, 5)

        # Use the sizers
        panel.SetSizerAndFit(border)  
        self.SetSizerAndFit(sizer)
        
       # panel.SetSizer(sizer)
        
    def OnOpenSOF(self,e):
        global mainPath, sofFile
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
        
        for idx, btn in enumerate(self.btns):
            btn.filelist = [];
        
        f = open(sofFile, 'r')
        for line in f:
            line     = line.strip()
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
                mainPath = os.path.dirname(btn.filelist[0])
                btn.dirtxt.SetValue(mainPath)
                # Get file list
                mystring = [os.path.basename(fil) for fil in btn.filelist]
                text     = mystring[0];
                for texti in mystring[1:]:
                    text += ", "+texti
                btn.name.SetValue(os.path.basename(text))
            
       # wx.MessageBox('Nothing there! Come again later ;-)', 'Error', wx.OK | wx.ICON_ERROR)
        
    def OnGenSOF(self,e):
        global mainPath, sofFile, guiTitle
        # First things first: check values have been filled
        for idx, btn in enumerate(self.btns):
            if not btn.filelist:
                if btn.checkPresent:
                    wx.MessageBox('Please set first '+fileTypes[idx], 'Error', wx.OK | wx.ICON_ERROR)
                    return -1
        
        # Get file list
        text = "";
        for idx, btn in enumerate(self.btns):
            for texti in btn.filelist:
                text += texti+" "+fileTypes[idx]+"\n"
                
        if sofFile == "":
            sofFile = mainPath+'/'+ guiTitle+'.sof'
        print("writing SOF file "+sofFile)
        f = open(sofFile, 'w')
        f.writelines(text)
        f.close()
    
        
    def OnClose(self,e):
        self.Close(True)    
  
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
        self.dirtxt  =   wx.TextCtrl(panel, size=(400, -1))
        self.txt     = wx.StaticText(panel, label=typeFile, name=typeFile)
        self.dirdsc  = wx.StaticText(panel, label="dir",    name="dir")
        self.filedsc = wx.StaticText(panel, label="files",  name="files")
        self.name    =   wx.TextCtrl(panel, size=(400, -1))
        self.btn     =     wx.Button(panel, label="Select", name=typeFile)
        
        self.btn.Bind(wx.EVT_BUTTON, self.onButton)
        #sizer.Add(btn, 0, wx.ALL, 5)
        
        row += 1
        col  = 0
        sizer.Add(self.txt,     (row, col))
        col += 1
        sizer.Add(self.dirdsc,  (row, col))
        col += 1
        sizer.Add(self.dirtxt,  (row, col))
        row += 1
        col  = 1
        sizer.Add(self.filedsc, (row, col))
        col += 1
        sizer.Add(self.name,    (row, col))
        col += 1
        sizer.Add(self.btn,     (row, col))
        
    #----------------------------------------------------------------------
    def onButton(self, event):
        global mainPath, sofFile
        """
        This method is fired when its corresponding button is pressed
        """
        button = event.GetEventObject()
        dlg    = wx.FileDialog(
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
        
            mainPath = os.path.dirname(self.filelist[0])
            if sofFile == "":
                sofFile = mainPath+'/'+ guiTitle+'.sof'
            self.dirtxt.SetValue(mainPath)
            # Get file list
            mystring = [os.path.basename(fil) for fil in self.filelist]
            text     = mystring[0];
            for texti in mystring[1:]:
                text += ", "+texti
            self.name.SetValue(os.path.basename(text))
            dlg.Destroy()
            print "You selected a "+self.typeFile
                            
#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    app   = wx.App(False)
    frame = displayGui()
    frame.Show()
    app.MainLoop()