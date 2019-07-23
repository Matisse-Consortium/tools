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
    def __init__(self,parent,headerOrFilename):
        
        if type(headerOrFilename)==type("") or  type(headerOrFilename)==type(u""):
            header=fits.getheader(headerOrFilename)
            filename=headerOrFilename
        else:
            header=headerOrFilename 
            filename = "Unkwnown"
        wx.Frame.__init__(self, parent, title=filename,size=(800,700))
        self.tableView=ObjectListView(self,wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER)     
        self.tableView.rowFormatter=self.setRowColor  
        cols=[ColumnDefn("Keyword","left",250,"keyword",minimumWidth=100),ColumnDefn("Value","left",250,"value",minimumWidth=100),ColumnDefn("Comment","left",290,"comment",minimumWidth=100)]
        
        self.tableView.SetColumns(cols)  

       
            
        self.keywords=[]
        for keyi in header.keys():
            self.keywords.append(headerKeyword(keyi,header[keyi],header.comments[keyi]))
        
  
        self.tableView.SetObjects( self.keywords)
 
    def setRowColor(self,listItem, data):
       
        c=listItem.GetId() % 2
        listItem.SetBackgroundColour(wx.Colour(250+3*c,250+3*c,250+3*c))
  
