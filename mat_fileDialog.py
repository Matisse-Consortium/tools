#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
author ame
'''
import wx
import os
from ObjectListView import ObjectListView, ColumnDefn
from astropy.io import fits
import subprocess

fvpath = "C:/fv/bin/fv.exe"
 
 
class fileViewerKeyword:
    def __init__(self,headerkeyword=None,name=None,checkheader=None,checkheader_cases=None,function=None,source=None):
        self.headerkeyword=headerkeyword
        self.name=name
        self.checkheader=checkheader
        self.checkheader_cases=checkheader_cases
        self.function=function
        self.source=source

    def evaluate(self,header=None,filename=None):
        res=""
        if (self.function and self.source=="header"):
            res=eval("{0}(header)".format(self.function)) 
        elif (self.function and self.source=="filename"):
            res=eval("{0}(filename)".format(self.function))            
        elif self.headerkeyword and header:
            if self.checkheader:
               check=header[self.checkheader] 
               for i in range(len(self.checkheader_cases)):
                   if check==self.checkheader_cases[i]:
                       res=header[self.headerkeyword[i]]
            else:
               res=header[self.headerkeyword]           
        return res

    def __str__(self):
        
        res="NAME={0}".format(self.name)
        if  self.headerkeyword:
            res+=" HEADERKEYWORD={0}".format(self.headerkeyword)
        if  self.checkheader:
            res+=" CHECKHEADER={0}".format(self.checkheader)
        if  self.checkheader_cases:
            res+=" CHECKHEADER8CASE={0}".format(self.checkheader_cases)     
        if  self.function:
            res+=" FUNCTION={0}".format(self.function)  
        if  self.source:
            res+=" SOURCE={0}".format(self.source)    
        return res


def matisseType(header):   
    res=""
    catg=None
    typ=None
    tech=None
    try:
        catg=header['HIERARCH ESO PRO CATG']       
    except:
        try:
            catg=header['HIERARCH ESO DPR CATG']
            typ=header['HIERARCH ESO DPR TYPE']
            tech=header['HIERARCH ESO DPR TECH']
        except:
            pass        
    if catg=="CALIB" and typ=="DARK,DETCAL" and tech=="IMAGE":
        res="DARK"
    elif catg=="CALIB" and typ=="FLAT,DETCAL" and tech=="IMAGE":
        res="FLAT"
    elif catg=="CALIB" and typ=="DARK" and tech=="SPECTRUM":
        res="OBSDARK"   
    elif catg=="CALIB" and typ=="FLAT" and tech=="SPECTRUM":
        res="OBSFLAT"    
    elif catg=="CALIB" and typ=="DARK,WAVE" and tech=="IMAGE":
        res="DISTOR_HOTDARK" 
    elif catg=="CALIB" and typ=="SOURCE,WAVE" and tech=="IMAGE":
        res="DISTOR_IMAGES"     
    elif catg=="CALIB" and typ=="SOURCE,LAMP" and tech=="SPECTRUM":
        res="SPECTRA_HOTDARK" 
    elif catg=="CALIB" and typ=="SOURCE,WAVE" and tech=="SPECTRUM":
        res="SPECTRA_IMAGES"
    elif catg=="CALIB" and typ=="DARK,FLUX" and tech=="IMAGE":
        res="KAPPA_HOTDARK"
    elif catg=="CALIB" and typ=="SOURCE,FLUX" and tech=="IMAGE":
        res="KAPPA_SRC" 
    elif catg=="CALIB" and typ=="SKY" and tech=="IMAGE":
        res="KAPPA_SKY"
    elif catg=="CALIB" and typ=="OBJECT,FLUX" and tech=="IMAGE":
        res="KAPPA_OBJ"
    elif catg=="SCIENCE" and typ=="OBJECT" and tech=="IMAGE":
        res="TARGET_RAW"        
    elif catg=="CALIB" and typ=="OBJECT" and tech=="IMAGE":
        res="CALIB_RAW"         
    elif catg=="CALIB" and typ=="DARK,IMB" and tech=="IMAGE":
        res="IM_COLD"        
    elif catg=="CALIB" and typ=="FLAT,IME" and tech=="IMAGE":
        res="IM_FLAT"
    elif catg=="CALIB" and typ=="DARK,IME" and tech=="IMAGE":
        res="IM_DARK"
    elif catg=="CALIB" and typ=="DARK,FLAT" and tech=="IMAGE":
        res="IM_PERIODIC" 
    elif catg=="TECHNICAL" and typ=="DARK,FOCUS" and tech=="IMAGE":
        res="REF_HOTDARK" 
    elif catg=="TECHNICAL" and typ=="SOURCE,FOCUS" and tech=="IMAGE":
        res="IM_REF"
    elif catg=="CALIB" and typ=="DARK" and tech=="INTERFEROMETRY":
        res="HOT_DARK"
    elif catg=="CALIB" and typ=="SOURCE" and tech=="INTERFEROMETRY":
        res="CALIB_SRC_RAW"     
    else:
        res=catg        
    return res
 
 
  
keywords=[]
keywords.append(fileViewerKeyword(function="matisseType",source="header",name="DoCatg"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET CHIP NAME",name="Detector"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET NDIT",name="NDIT"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET SEQ1 DIT",name="DIT"))
keywords.append(fileViewerKeyword(headerkeyword=["HIERARCH ESO INS PIL NAME","HIERARCH ESO INS PIN NAME"],checkheader="HIERARCH ESO DET CHIP NAME",checkheader_cases=["HAWAII-2RG","AQUARIUS"],name="Mode"))
keywords.append(fileViewerKeyword(headerkeyword=["HIERARCH ESO INS DIL NAME","HIERARCH ESO INS DIN NAME"],checkheader="HIERARCH ESO DET CHIP NAME",checkheader_cases=["HAWAII-2RG","AQUARIUS"],name="Resolution"))


 
 
class matisseFile(object):
    def __init__(self,filename,folder):
        self.filename=filename
        self.folder=folder
        
        self.header=None
        self.isFits=True
        self.isMatisse=True
        
        
        if filename.endswith(".fits"):
            try:
                self.header=fits.getheader( self.folder+"/"+self.filename)                       
            except:
                self.isFits=False
                self.isMatisse=False
        else:
            self.isFits=False
            self.isMatisse=False    
                
        for keywordi in keywords:
            try:
                self.__dict__[keywordi.name]=keywordi.evaluate(header=self.header,filename=folder+"/"+filename)
                
            except:
                self.__dict__[keywordi.name]=""
                self.isMatisse=False
                
matisseColor={
"DARK":wx.Colour(255,255,0),
"FLAT":wx.Colour(150,0,0),
"OBSDARK":wx.Colour(150,0,150),
"OBSFLAT":wx.Colour(0,0,150),
"DISTOR_HOTDARK":wx.Colour(0,0,200),
"DISTOR_IMAGES":wx.Colour(0,220,0),
"SPECTRA_HOTDARK":wx.Colour(150,150,0),
"SPECTRA_IMAGES":wx.Colour(0,150,0),
"SPECTRA_HOTDARK":wx.Colour(200,0,0),
"SPECTRA_IMAGES":wx.Colour(0,0,50),
"KAPPA_HOTDARK":wx.Colour(50,0,0),
"KAPPA_SRC" :wx.Colour(50,50,0),
"KAPPA_SKY":wx.Colour(150,0,50),
"KAPPA_OBJ":wx.Colour(150,50,0),
"TARGET_RAW":wx.Colour(0,50,150),
"CALIB_RAW":wx.Colour(150,0,0),
"IM_COLD":wx.Colour(0,0,150),
"IM_FLAT":wx.Colour(0,200,0),
"IM_PERIODIC":wx.Colour(0,200,50), 
"REF_HOTDARK":wx.Colour(150,50,100), 
"IM_REF":wx.Colour(150,0,50),
"HOT_DARK":wx.Colour(50,50,50),
"CALIB_SRC_RAW":wx.Colour(0,50,0) 
}


    
class mat_FileDialog(wx.Dialog):
  
    def __init__(self, parent, title):
        super(mat_FileDialog, self).__init__(parent, title=title, size=(1200, 600))
            
        self.InitUI()
        self.Centre()
        self.Show()  
        self.path=""
        
    def InitUI(self):
    
        panel = wx.Panel(self)

        font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
        font.SetPointSize(9)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)

        splitter=wx.SplitterWindow(panel)    
        self.dirTree = wx.GenericDirCtrl(splitter, style=wx.DIRCTRL_DIR_ONLY|wx.DIRCTRL_EDIT_LABELS)
        self.fileList = ObjectListView(splitter,wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        
        cols=[ColumnDefn("Filename","left",300,"filename",minimumWidth=250)]
        cols.extend([ColumnDefn(keywordi.name,"left",75,keywordi.name,minimumWidth=75) for keywordi in keywords])  
        
        self.fileList.SetColumns(cols)       
        self.fileList.rowFormatter=self.setRowColor    
        self.fileList.AutoSizeColumns()
        
        splitter.SplitVertically(self.dirTree,self.fileList)
        splitter.SetMinimumPaneSize(400)
        vbox.Add(splitter, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)        
        vbox.AddSpacer(10)
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)  
        self.addDirButton=wx.Button(panel,label="+")
        hbox2.Add(self.addDirButton, proportion=0.1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        self.pathText=wx.TextCtrl(panel,style=wx.TE_READONLY)
        hbox2.Add(self.pathText, proportion=1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        self.filterBox=wx.ComboBox(panel,choices=["All Files","Fits Files","Matisse Files"],style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.filterBox.SetValue("All Files")
        hbox2.Add(self.filterBox, proportion=0.2,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)        
        self.cancelButton = wx.Button(panel,wx.ID_CANCEL)
        hbox2.Add(self.cancelButton,proportion=0.1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        self.okButton = wx.Button(panel,wx.ID_OK )
        hbox2.Add(self.okButton,proportion=0.1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        vbox.Add(hbox2, proportion=0.2, flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)        
        vbox.AddSpacer(10)
       
        panel.SetSizer(vbox)
       
        
        tree = self.dirTree.GetTreeCtrl()
        self.Bind(wx.wx.EVT_TREE_SEL_CHANGED, self.dirChanged,tree)
        self.Bind(wx.EVT_BUTTON, self.okClicked, self.okButton)
        self.Bind(wx.EVT_BUTTON, self.cancelClicked, self.cancelButton) 
        self.Bind(wx.EVT_COMBOBOX, self.filterChanged, self.filterBox) 
        self.Bind(wx.EVT_LIST_ITEM_SELECTED,self.fileSelected,self.fileList)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK,self.fileListRightClicked,self.fileList)
        self.Bind(wx.EVT_BUTTON, self.addDirClicked, self.addDirButton)
        
    def setRowColor(self,listItem, data):
       
        if data.DoCatg:
            try:
                col=matisseColor[data.DoCatg]
            except:
                col=wx.RED
        else:
            col=wx.Colour(180,180,180)
        listItem.SetTextColour(col)
        c=listItem.GetId() % 2
        listItem.SetBackgroundColour(wx.Colour(250+3*c,250+3*c,250+3*c))
       
    def dirChanged(self,treeEvent=None):
       self.pathText.SetValue("")
       newDir= self.dirTree.GetPath()    
       files=os.listdir(newDir)     
       matisseFileList=[]
       filt=self.filterBox.GetValue()  
            
       for filei in files:        
           if os.path.isfile(newDir+"/"+filei):
               matFile=matisseFile(filei,newDir)
               if filt=="All Files":
                   Append=True
               elif filt=="Fits Files" and matFile.isFits:
                   Append=True
               elif filt=="Matisse Files" and matFile.isMatisse:
                   Append=True
               else:
                   Append=False
               if Append:                               
                   matisseFileList.append(matFile)
       self.fileList.SetObjects(matisseFileList)
       self.fileList.AutoSizeColumns()
       
       for icol in range(len(keywords)+1):
           self.fileList.SetColumnWidth(icol,wx.LIST_AUTOSIZE)
           wc=self.fileList.GetColumnWidth(icol)
           if wc<self.fileList.columns[icol].minimumWidth:
               self.fileList.SetColumnWidth(icol,self.fileList.columns[icol].minimumWidth)

       
    def fileSelected(self,event):
        nfiles=self.fileList.GetSelectedItemCount()
        itemNum=-1
        selectedItems=[]
        for ifile in range(nfiles):
            itemNum=self.fileList.GetNextSelected(itemNum)
            selectedItems.append(self.fileList.GetItem(itemNum).GetText())
        self.path=selectedItems
        txt=""
        for itemi in selectedItems:
            txt=txt+itemi+" "
        self.pathText.SetValue(txt)
        
    def fileListRightClicked(self,event):
        menu = wx.Menu()
        menu.Append( 0, "Show Header" )
        wx.EVT_MENU( menu, 0, self.showHeader)
        menu.Append( 1, "Show IMAGING_DETECTOR")        
        wx.EVT_MENU( menu, 1, self.showImagingDetector)
        menu.Append( 2, "Show IMAGING_DATA")
        wx.EVT_MENU( menu, 2, self.showImagingData)
        menu.Append( 3, "Open with fv" )
        wx.EVT_MENU( menu, 3, self.openWithFv)
        
        self.fileList.PopupMenu( menu, event.GetPoint())
    def okClicked(self,event):
        self.EndModal(wx.ID_OK)
        
        
    def cancelClicked(self,event):
         self.EndModal(wx.ID_CANCEL)
         
    def addDirClicked(self,event):
        path=self.dirTree.GetPath()
        os.makedirs(path+"/"+"newDir")
        self.dirTree.ReCreateTree()
        self.dirTree.ExpandPath(path+"/"+"newDir")
    
    def filterChanged(self,event):
        self.dirChanged()
  
    def GetPath(self):
        return [self.dirTree.GetPath()+'/'+pathi for pathi in self.path]
        
    def showHeader(self,event):
        print "show header {0}".format(self.GetPath())
        
    def showImagingDetector(self,event):
        print "show IMAGING_DETECTOR {0}".format(self.GetPath())
        
    def showImagingData(self,event):
        print "show IMAGING_DATA {0}".format(self.GetPath())
        
    def openWithFv(self,event):
        print "Open with fv {0}".format(self.GetPath())
        for filei in self.path:
            if filei.endswith('.fits'):
                subprocess.Popen([fvpath,self.dirTree.GetPath()+'/'+filei])       
    
if __name__ == '__main__':
  
    app = wx.App()
    openFileDialog=mat_FileDialog(None, title='Open file')
    if openFileDialog.ShowModal()==wx.ID_OK:
        print openFileDialog.GetPath()
    openFileDialog.Destroy()
    app.MainLoop()
    app.Destroy()
    
    
    #GenericDirCtrl
    #ListView
