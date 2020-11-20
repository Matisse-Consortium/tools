#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 19:31:29 2020

@author: ame
"""

import wx
import os
from ObjectListView import ObjectListView, ColumnDefn
import argparse
import pickle
import sys
import numpy as np
#from mat_reducPkl import mat_dataReduction

_status=["NOT STARTED  ","PROCESSING...","PROCESSED  ","POSTPROCESSING...","COMPLETED    "]


class mat_dataReduction:

    def __init__(self):
        self.files=[]
        self.dirCalib=None
        self.dirResult=None
        self.paramL={}
        self.paramN={}
        self.postCmd=[]
        self.maxIter=4
        self.nbCore=1
        self.skipN=False
        self.skipL=False
        self.status=0
        self.filename=None
        self.logfile=None
    """    
    def save(self,filename=None):
        if filename==None:
            if self.filename!=None:
                filename= self.filename
            else:
                print("Error: filename not given")
                return
        f=open(filename,"wb")
        f.write(pickle.dumps(self.__dict__, protocol=2))
        f.close()
        self.filename=filename
    """
    def load(self,filename):
        f=open(filename,"rb")
        dataPickle = f.read()
        f.close()
        self.__dict__ = pickle.loads(dataPickle)
        f.close()
        self.filename=filename
        
    def getShortFilename(self):
        return self.filename.replace("\\","/").split("/")[-1]
    
    def getStatusName(self):
        return _status[self.status]

        

class mat_reductionGUI(wx.Dialog):

    def __init__(self, parent, dirPkl):
        super(mat_reductionGUI, self).__init__(parent, title="MATISSE reduction tool",
                                               style=wx.DEFAULT_FRAME_STYLE, size=(600, 900))
        
        
        self.dir=os.path.abspath(dirPkl)
        self.pklfiles=[os.path.join( self.dir,fi) for fi in os.listdir( self.dir) if ".pkl" in fi]
        
        self.reduc=[]
        for fi in self.pklfiles:
            r=mat_dataReduction()    
            r.load(fi)
            #r.status=np.random.randint(0,5)
            self.reduc.append(r)
   
        self.InitUI()
        self.Centre()
        

    def InitUI(self):
        
        panel = wx.Panel(self)
        font  = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
        font.SetPointSize(12)
        
        vbox   = wx.BoxSizer(wx.VERTICAL)
        hbox   = wx.BoxSizer(wx.HORIZONTAL)
        
        
        self.reducList = ObjectListView(panel,wx.ID_ANY, style=wx.LC_REPORT)
        cols=[ ColumnDefn("Filename",   "left",200,"getShortFilename",   minimumWidth=20),
               ColumnDefn("dirResult","left",300,"dirResult",minimumWidth=20),
               ColumnDefn("status",   "left",100, "getStatusName",  minimumWidth=20)]
        self.reducList.SetColumns(cols)
        self.reducList.rowFormatter=self.rowFormatter
        self.reducList.AutoSizeColumns()
        self.reducList.SortBy(0, ascending=False)
        self.reducList.SetFont(font)
        self.reducList.SetObjects(self.reduc)

        
        vbox.Add(self.reducList,proportion=10,flag=wx.LEFT|wx.RIGHT|wx.EXPAND|wx.TOP, border=10)
        
        
        self.but1    = wx.Button(panel,label="but1",style=wx.BU_EXACTFIT)
        self.but2    = wx.Button(panel,label="but2",style=wx.BU_EXACTFIT)
        hbox.Add(self.but1,proportion=1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND|wx.TOP, border=10)
        hbox.Add(self.but2,proportion=1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND|wx.TOP, border=10)        
        vbox.Add(hbox,proportion=1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND|wx.TOP, border=10)
       
        
        panel.SetSizer(vbox)
        
        self.Bind(wx.EVT_CLOSE, self.byebye)


    def rowFormatter(self,listItem, data):
        listItem.SetTextColour(wx.BLACK)
        if data.status == 0:
            listItem.SetBackgroundColour(wx.Colour(100,100,255))
        elif data.status<4:
            listItem.SetBackgroundColour(wx.Colour(255,255,100))
        else:
            listItem.SetBackgroundColour(wx.Colour(100,255,100))
        

    def byebye(self,event):
        self.Destroy()
        return 0
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='mat_reductionGUI.py')
    
    parser.add_argument('dir', default="",  \
    help='The path to the directory containing your raw data.')


    
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_reductionGUI.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        print("\n     Example : python mat_reductionGUI.py .")
        sys.exit(0)
    
        
    try :
        app = wx.App()
    except:
        pass
    
    dir0=os.path.abspath(args.dir)
    os.chdir(dir0)
    
    gui = mat_reductionGUI(None,dir0)
    gui.Show()
    app.MainLoop()


