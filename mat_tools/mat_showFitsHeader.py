# -*- coding: utf-8 -*-
"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Wed Dec 21 15:28:51 2016
@author: ame

Fits header display

This software is a computer program whose purpose is to show oifits
files from the MATISSE instrument.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software. 

You can use, modify and/ or redistribute the software under the
terms of the CeCILL license as circulated by CEA, CNRS and INRIA at
the following URL "http://www.cecill.info". You have a copy of the
licence in the LICENCE.md file.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""

import wx
from ObjectListView import ObjectListView, ColumnDefn
from astropy.io import fits

class headerKeyword(object):
    def __init__(self,key,value,comment):
        self.keyword=key
        self.value=value
        self.comment=comment

class mat_showFitsHeader(wx.Frame):
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

