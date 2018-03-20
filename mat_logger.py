#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
  $Id:  $

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur

  Created on 2018-03-15
  @author: ame

  This software is a computer program whose purpose is to produce a smart log
  for the MATISSE instrument.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
"""

# Import necessary files
from libAutoPipeline import matisseType
from mat_fileDialog import mat_FileDialog
from mat_fileDialog import identifyFile
import wx
import os
from ObjectListView import ObjectListView, ColumnDefn
from astropy.io import fits
from fitsheaderviewer import fitsHeaderViewer
import distutils.spawn
import sys
import pickle
import socket



# Set useful paths
fvpath    = distutils.spawn.find_executable("fv")
iconspath = os.path.join(os.path.dirname(__file__),"icons")

#various dir0 for debug 
if socket.gethostname()=="pc-amepro":
    dir0="D:\\Documents\\Travail\\MATISSE\\"
else:
    #on wmtpipeline
    dir0="/data/RawDataMatisse/"




def findHeaderKeyword(h,key):
    try :
        res=h[key]
    except:
        print ("keyword {0} not found".format(key))
        res=""
    return res

    
###############################################################################

class mat_logData():
      def __init__(self,tplstart,tplid,target,progid,nbFiles,nexp,comment,firstFileData,status):
        self.tplstart    = tplstart
        self.tplid       = tplid
        self.target      = target
        self.progid      = progid
        self.nbFiles     = nbFiles
        self.nexp        = nexp
        self.comment     = comment
        self.listOfFiles = [firstFileData]
        self.status      = status
    
###############################################################################
        
class mat_fileData():
      def __init__(self,filename,h):
        self.filename    = filename
        self.doCatg      = matisseType(h)
        self.date    = findHeaderKeyword(h,'DATE')
        self.expno= "{0}/{1}".format(findHeaderKeyword(h,'HIERARCH ESO TPL EXPNO'),
                                   findHeaderKeyword(h,'HIERARCH ESO TPL NEXP'))
 
    
###############################################################################

class mat_logger(wx.Dialog):

    def __init__(self, parent, date):
        super(mat_logger, self).__init__(parent, title="MATISSE Log of {0}".format(date), size=(1200, 900))
        self.date=date
        self.logfilename=dir0+"/"+self.date+'/mat_logDialog.pkl'
        self.tplList=[]       
        self.tplListObj=[]
        self.fileList=[]
        
        self.InitUI()
        self.Centre()
        self.Show()
        self.path=""
            
#------------------------------------------------------------------------------

    def InitUI(self):

        
        panel = wx.Panel(self)

        font  = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
            
        font.SetPointSize(9)

        hbox   = wx.BoxSizer(wx.HORIZONTAL)   
        vbox   = wx.BoxSizer(wx.VERTICAL)   
        vbox2  = wx.BoxSizer(wx.VERTICAL)   
       
        self.tplListWidget = ObjectListView(panel,wx.ID_ANY, style=wx.LC_REPORT)   
        cols=[ ColumnDefn("TPL START","left",160,"tplstart",minimumWidth=20),
               ColumnDefn("TPL ID","left",140,"tplid",minimumWidth=20),
               ColumnDefn("TARGET","left",80,"target",minimumWidth=20),
               #ColumnDefn("PROG. ID","left",70,"progid",minimumWidth=20),                        
               ColumnDefn("Files","left",30,"nbFiles",minimumWidth=20),                         
               ColumnDefn("NEXP","left",30,"nexp",minimumWidth=20)]    
        self.tplListWidget.SetColumns(cols)
        self.tplListWidget.rowFormatter=self.setRowColor
        self.tplListWidget.AutoSizeColumns()
        self.tplListWidget.SortBy(0, ascending=True)
                
              
        self.commentTxtCtrl= wx.TextCtrl(panel,style=wx.TE_MULTILINE) 
        
        self.fileListWidget =  ObjectListView(panel,wx.ID_ANY, style=wx.LC_REPORT) 
        cols2=[ColumnDefn("Date","left",200,"date",minimumWidth=20),
               ColumnDefn("File name","left",340,"filename",minimumWidth=20),
               ColumnDefn("DO CATG","left",100,"doCatg",minimumWidth=20),
               ColumnDefn("EXPNO","left",50,"expno",minimumWidth=20)]
        self.fileListWidget.SetColumns(cols2)
        self.fileListWidget.AutoSizeColumns()
        self.fileListWidget.SortBy(0, ascending=True)
    
        self.updateBut=wx.Button(panel,label="Update",style=wx.BU_EXACTFIT)
        self.updateBut.Bind(wx.EVT_BUTTON, self.updateClicked)
        
       
        
        vbox2.Add(self.fileListWidget, proportion=3, flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)
        vbox2.Add(self.commentTxtCtrl, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)
        
        hbox.Add(self.tplListWidget, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND,border=10)  
        hbox.Add(vbox2, proportion=2,flag=wx.LEFT|wx.RIGHT|wx.EXPAND);
   

        vbox.Add(hbox,proportion=1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)
        vbox.Add(self.updateBut,proportion=0.1,flag=wx.LEFT|wx.RIGHT|wx.EXPAND)



        panel.SetSizer(vbox)

         
        self.Bind(wx.EVT_LIST_ITEM_SELECTED,self.tplSelected,self.tplListWidget)
        self.Bind(wx.EVT_LIST_ITEM_SELECTED,self.updateClicked,self.updateBut)
        self.Bind(wx.EVT_TEXT,self.commentChanged,self.commentTxtCtrl)
        
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK,self.fileListRightClicked,self.fileListWidget)
        
        
        
        if os.path.isfile(self.logfilename):
            print ("log file exist for {0} ... loading".format(self.date))
            pik=open(self.logfilename, 'rb')
            self.tplList    = pickle.load(pik) 
            self.tplListObj = pickle.load(pik)
            self.fileList   = pickle.load(pik)       
        else:
            print ("No log file for {0} ...".format(self.date))
                         
        self.getInfosFromNight()
        
#------------------------------------------------------------------------------                
    def fileListRightClicked(self,event):
        menu = wx.Menu()
        m1=menu.Append( 0, "Show Header" )
        menu.Bind(wx.EVT_MENU,self.showHeader,m1)
        #wx.EVT_MENU( menu, 0, self.showHeader)
        m2=menu.Append( 1, "Show RAW DATA")
        menu.Bind(wx.EVT_MENU,self.showRawData,m2)
        #wx.EVT_MENU( menu, 1, self.showRawData)
        self.fileListWidget.PopupMenu( menu, event.GetPoint())
        
#------------------------------------------------------------------------------
        
    def showHeader(self,event):
        print("Show header")
        itemNum=self.fileListWidget.GetNextSelected(-1)            
        idx=self.fileListWidget.GetItem(itemNum).GetData()
        l=self.fileListWidget.GetObjects()
        filename=l[idx].filename

        fitsHeaderViewer(dir0+"/"+self.date+'/'+filename)       
                 
#------------------------------------------------------------------------------
                
    def showRawData(self,event):
        print ("Show IMAGING_DATA")    
#------------------------------------------------------------------------------
                
    def saveBackup(self,event):
        print ("Saving pkl file")    
        
#------------------------------------------------------------------------------

    def saveData(self):
        pik     = open(self.logfilename, 'wb')
        pickle.dump(self.tplList, pik, pickle.HIGHEST_PROTOCOL)
        pickle.dump(self.tplListObj, pik, pickle.HIGHEST_PROTOCOL)
        pickle.dump(self.fileList, pik, pickle.HIGHEST_PROTOCOL)  
#------------------------------------------------------------------------------
        
    def updateClicked(self,event):
        self.getInfosFromNight()
        self.tplListWidget.SortBy(0, ascending=True)
#------------------------------------------------------------------------------
        
        
    def commentChanged(self,event):
        itemNum=self.tplListWidget.GetNextSelected(-1)
        if (itemNum!=-1):
            txt=self.tplListWidget.GetItemText(itemNum)
            idx=self.tplList.index(txt)   
            self.tplListObj[idx].comment=self.commentTxtCtrl.GetValue()
            self.saveData()
#------------------------------------------------------------------------------
           
    def tplSelected(self,event):
        #nfiles=self.tplListWidget.GetSelectedItemCount()
        itemNum=self.tplListWidget.GetNextSelected(-1)
        
        
        txt=self.tplListWidget.GetItemText(itemNum)
        idx=self.tplList.index(txt)
        self.fileListWidget.SetObjects(self.tplListObj[idx].listOfFiles) 
        self.commentTxtCtrl.SetValue(self.tplListObj[idx].comment)
        
#------------------------------------------------------------------------------

    def getInfosFromNight(self):
                  
        files=os.listdir(dir0+self.date)  
        for filei in files:
            if not(filei in self.fileList):
                if filei.endswith(".fits"):
                    try :
                        h=fits.getheader( dir0+"/"+self.date+"/"+filei)
                        try : 
                            tplstart=findHeaderKeyword(h,'HIERARCH ESO TPL START')
                            #tplstart=h['HIERARCH ESO TPL START']
                            print ("{0} ==> {1}".format(filei,tplstart))
                            if (tplstart in self.tplList):           
                                i=self.tplList.index(tplstart)
                                print("not new")
                                self.tplListObj[i].nbFiles+=1
                                #if self.tplListObj[i].nexp<findHeaderKeyword(h,['HIERARCH ESO TPL NEXP']):
                                #    self.tplListObj[i].nexp=findHeaderKeyword(h,['HIERARCH ESO TPL NEXP'])
                                if self.tplListObj[i].nexp<h['HIERARCH ESO TPL NEXP']:
                                    self.tplListObj[i].nexp=h['HIERARCH ESO TPL NEXP']
                                self.tplListObj[i].listOfFiles.append(mat_fileData(filei,h))
                                self.fileList.append(filei)      
                            else:
                                print("new")
                                self.tplList.append(tplstart)    
                                try :
                                    target= findHeaderKeyword(h,'HIERARCH ESO OBS TARG NAME')                                    
                                except:
                                    target=""
                                try :
                                    progid= findHeaderKeyword(h,'HIERARCH ESO OBS PROG ID')
                                except:
                                    progid=""
                             
                                    
                                self.tplListObj.append(mat_logData(tplstart,
                                    findHeaderKeyword(h,'HIERARCH ESO TPL ID'),target,progid,1,
                                    findHeaderKeyword(h,'HIERARCH ESO TPL NEXP'),"No Comment",
                                    mat_fileData(filei,h),"Started"))
                                self.fileList.append(filei)  
                        except:
                            print ("skipping file {0}".format(filei)    )                                                
                    except:    
                        print( "skipping file {0} : not valid fits file".format(filei))
                else:
                    print ("skipping file {0} : not fits file".format(filei))
        self.tplListWidget.SetObjects(self.tplListObj)
        self.saveData()        
#------------------------------------------------------------------------------


    # Colour of text
    def setRowColor(self,listItem, data):

        if data.tplid == "MATISSE_hyb_obs":
                txtcol=wx.Colour(116,196,147)
        elif data.tplid == "MATISSE_img_acq":
                txtcol=wx.Colour(116,147,196)
        else:
                txtcol=wx.Colour(155,155,155)      
         
            
        listItem.SetBackgroundColour(txtcol)
    
###############################################################################
    
###############################################################################
   
if __name__ == '__main__':
    
    # Scan the command line arguments
    listArg = sys.argv
    name_file = []
    for elt in listArg:
        if ('--help' in elt):
            print( "Usage: mat_logger.py <name of night directory>")
            sys.exit(0)
        elif len(listArg) == 2:
            name_file = sys.argv[1]
            print(name_file)
            
    # If no argument is given, then open the file dialog to select directory
    if not name_file:
        app2 = wx.App()
        print("No input directory given, running file selector...")
        openFileDialog = mat_FileDialog(None, 'Open a night directory',"lmk,")
        if openFileDialog.ShowModal() == wx.ID_OK:
            name_file = openFileDialog.GetPaths()[0]
            print( name_file)
        openFileDialog.Destroy()
        app2.MainLoop()
        app2.Destroy()
    
    try :
        app = wx.App()
    except:
        pass
    
    openLogger = mat_logger(None,name_file)
    if openLogger.ShowModal() == wx.ID_OK:
        print ("OK")
    openLogger.Destroy()
    app.MainLoop()
    app.Destroy()
