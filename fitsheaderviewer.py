# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:28:51 2016

@author: ame
"""

import wx
from ObjectListView import ObjectListView, ColumnDefn 
from astropy.io import fits

class headerKeyword(object):
    def __init__(self,key,value,comment):
        self.keyword=key
        self.value=value
        self.comment=comment

class fitsHeaderViewer(wx.Frame):
    def __init__(self,headerOrFilename):
        super(fitsHeaderViewer, self).__init__(None,title="FitsHeaderViewer",size=(500, 400))  
       
        self.tableView=ObjectListView(self,wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER)     
        self.tableView.rowFormatter=self.setRowColor  
        cols=[ColumnDefn("Keyword","left",150,"keyword",minimumWidth=100),ColumnDefn("Value","left",150,"value",minimumWidth=100),ColumnDefn("Comment","left",150,"comment",minimumWidth=100)]
        
        self.tableView.SetColumns(cols)  

        if type(headerOrFilename)==type("") or  type(headerOrFilename)==type(u""):
            header=fits.getheader(headerOrFilename)
        else:
            header=headerOrFilename 
            
        self.keywords=[]
        for keyi in header.keys():
            self.keywords.append(headerKeyword(keyi,header[keyi],header.comments[keyi]))
        
  
        self.tableView.SetObjects( self.keywords)
        

        
        self.Centre()
        self.Show()
 
    def setRowColor(self,listItem, data):
       
        c=listItem.GetId() % 2
        listItem.SetBackgroundColour(wx.Colour(250+3*c,250+3*c,250+3*c))
       
if __name__ == '__main__':
    dir0="D:\\Documents\\Travail\\Etoiles\\AMBER_BE\\alpha ara\data\\"
    filename="PRODUCT_Alpha_Arae_1.89-2.59micron_2007-07-28T05_51_13.5574.fits"
    app = wx.App()
    f=fitsHeaderViewer(dir0+filename)
    f.Show()
    app.MainLoop()

         