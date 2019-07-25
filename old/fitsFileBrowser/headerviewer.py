# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:28:51 2016

@author: ame
"""

from PySide import QtCore, QtGui
import sys
from astropy.io import fits



class headerModel(QtCore.QAbstractTableModel):
    def __init__(self,header):
        super(headerModel, self).__init__(None) 
        self._data=[]
        
        for key in header.keys():
            if key!='':
                self._data.append([key,header[key],header.comments[key]])

        
    def rowCount(self,parent=QtCore.QModelIndex()):
        return len(self._data)
        
    def columnCount(self,parent=QtCore.QModelIndex()):
        return 3
        
    def data(self,index,role = QtCore.Qt.DisplayRole):
      
        if role == QtCore.Qt.DisplayRole: 
            return self._data[index.row()][index.column()]
        elif role == QtCore.Qt.EditRole:
            return self._data[index.row()][index.column()]              
        elif role==QtCore.Qt.TextColorRole: 
            return QtGui.QColor(QtCore.Qt.black)             
        else:
            return None
            
    def headerData(self, col, orientation, role=QtCore.Qt.DisplayRole):
      if orientation==QtCore.Qt.Horizontal:

          if role==QtCore.Qt.DisplayRole:
            if col==0:
                return 'KEYNAME'
            if col==1:
                return 'VALUE'
            if col==2:
                return 'COMMENTS'
          else:
              return None
               
        


class fitsHeaderViewer(QtGui.QWidget):
    def __init__(self,header):
        super(fitsHeaderViewer, self).__init__(None)  

        self.headerModel=headerModel(header)
        self.tableView=QtGui.QTableView()
        self.tableView.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows);
        self.tableView.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.tableView.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch )
        self.tableView.horizontalHeader().setStyleSheet("::section{background:rgb(250,250,250);}")
        self.tableView.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.tableView.setModel(self.headerModel)
        self.layout=QtGui.QVBoxLayout()
        self.layout.addWidget(self.tableView)
        self.setLayout(self.layout)
        self.resize(600,800)  
        
if __name__ == "__main__":  
    try:
        app = QtGui.QApplication(sys.argv) 
        
    except:
        app = QtGui.QApplication.instance()
    
    dir="D:/Documents/Travail/MATISSE/testplan/2016-12-20 franges/"
    filename=dir+"LIGHT_BSN1234.fits"
    header=fits.getheader(filename)
    w=fitsHeaderViewer(header)
    w.show()
   
    app.exec_()                