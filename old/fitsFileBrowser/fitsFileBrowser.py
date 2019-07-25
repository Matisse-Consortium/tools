# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 13:37:33 2016

@author: ame
"""

from PySide import QtCore, QtGui
import sys
import os
import datetime
from astropy.io import fits
from ui_fitsFileBrowser import *
import subprocess
import headerviewer  



class sofFileParameter:
    def __init__(self, docatg=None,minFile=1,maxFile=1):
        self.docatg=docatg
        self.minFile=minFile
        self.maxFile=maxFile
        

    

class fileViewerKeywordWidget(QtGui.QWidget):
    def __init__(self,keyword=None):
        super(fileViewerKeywordWidget, self).__init__(None)  
        
        
        if not(keyword):
            keyword=fileViewerKeyword()
        self.layout=QtGui.QHBoxLayout() 
          
        self.nameL=QtGui.QLabel("NAME")           
        self.nameW=QtGui.QTextEdit(keyword.name)     
        self.layout.addWidget(self.nameL)
        self.layout.addWidget(self.nameW)
        
        self.headerkeywordL=QtGui.QLabel("HEADERKEYWORD")           
        self.headerkeywordW=QtGui.QListWidget()     
        if type(keyword.headerkeyword)==type([]):    
            for item in keyword.headerkeyword:
                self.headerkeywordW.addItem(item)
        else:
            self.headerkeywordW.addItem(keyword.headerkeyword)
        self.layout.addWidget(self.headerkeywordL)
        self.layout.addWidget(self.headerkeywordW)            
            
        self.checkheaderL=QtGui.QLabel("CHECKHEADER")           
        self.checkheaderW=QtGui.QListWidget()     
        if type(keyword.checkheader)==type([]):    
            for item in keyword.checkheader:
                self.checkheaderW.addItem(item)
        else:
            self.checkheaderW.addItem(keyword.checkheader)                                 
        self.layout.addWidget(self.checkheaderL)
        self.layout.addWidget(self.checkheaderW)

        self.checkheader_casesL=QtGui.QLabel("CHECKHEADER_CASES")           
        self.checkheader_casesW=QtGui.QListWidget()     
        if type(keyword.checkheader_cases)==type([]):    
            for item in keyword.checkheader_cases:
                self.checkheader_casesW.addItem(item)
        else:
            self.checkheader_casesW.addItem(keyword.checkheader_cases)                                 
        self.layout.addWidget(self.checkheader_casesL)
        self.layout.addWidget(self.checkheader_casesW)  

        self.functionL=QtGui.QLabel("FUNCTION")           
        self.functionW=QtGui.QTextEdit(keyword.function)                                   
        self.layout.addWidget(self.functionL)
        self.layout.addWidget(self.functionW)  
        
        self.sourceL=QtGui.QLabel("SOURCE")           
        self.sourceW=QtGui.QTextEdit(keyword.source)                                   
        self.layout.addWidget(self.sourceL)
        self.layout.addWidget(self.sourceW)      
        
        self.setLayout(self.layout)
        
class EditfileViewerKeywordsDialog(QtGui.QDialog):
    def __init__(self,keywords=None):
        super(EditfileViewerKeywordsDialog, self).__init__(None)  
        self.layout=QtGui.QVBoxLayout()
        self.setLayout(self.layout)
        self.keywordsWidgetList=[]
        for keyword in keywords:
            print keyword
            w=fileViewerKeywordWidget(keyword)
            self.keywordsWidgetList.append(w)
            self.layout.addWidget(w)
        self.resize(1800,len(keywords)*20)        

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



class fileModel(QtCore.QAbstractTableModel):
    def __init__(self,dir=None):
        super(fileModel, self).__init__(None) 
        self._data=[]
        self.isMatisse=[]
        if dir:
            for filei in os.listdir(dir):
                isMatissei=True
                dati=[filei] 
                if os.path.isfile(os.path.join(dir, filei)):
                    if filei.endswith(".fits"):
                        header=fits.getheader(dir+"/"+filei) 
                        for keywordi in keywords:
                            try:
                                dati.append(keywordi.evaluate(header=header,filename=dir+"/"+filei))    
                            except:
                                dati.append('')
                                isMatissei=False
                    else:
                        isMatissei=False
                        for keywordi in keywords:
                            dati.append('')
                    dati.append(os.path.getmtime(dir+"/"+filei))
                    self._data.append(dati)
                    self.isMatisse.append(isMatissei)
                     
       
       
    def rowCount(self,parent=QtCore.QModelIndex()):
        return len(self._data)
        
    def columnCount(self,parent=QtCore.QModelIndex()):
        return 2+len(keywords)
        
    def data(self,index,role = QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole: 
            if index.column()<self.columnCount()-1:
                return self._data[index.row()][index.column()]    
            if index.column()==self.columnCount()-1:
                return datetime.datetime.fromtimestamp(self._data[index.row()][index.column()]).strftime("%Y-%m-%d %H:%M:%S")
                
        elif role == QtCore.Qt.EditRole:
            return self._data[index.row()][index.column()]
        elif role==QtCore.Qt.TextColorRole: 
            if self.isMatisse[index.row()]:
                return QtGui.QColor(QtCore.Qt.black) 
            else :
                return QtGui.QColor(QtCore.Qt.gray)               
        else:
            return None
            
    def headerData(self, col, orientation, role=QtCore.Qt.DisplayRole):
      if orientation==QtCore.Qt.Horizontal:

          if role==QtCore.Qt.DisplayRole:
              if col==0:
                return 'Filename'
              elif col<self.columnCount()-1:
                return keywords[col-1].name
              elif col==self.columnCount()-1:
                  return 'Last Modification'
          else:
              return None
       


class filterProxyModel(QtGui.QSortFilterProxyModel):
    def __init__(self,filt=True):
        super(filterProxyModel, self).__init__(None) 
        #self.filterData(filt)
        self.fiterOn=filt
    def filterAcceptsRow(self, row, parent):
        model = self.sourceModel()              
        data=model._data[row][0]
        return (model.isMatisse[row] or  not(self.fiterOn))

                
        return self.filterFits.exactMatch(data)
        
    def filterData(self,filt=True):
        self.fiterOn=filt
        self.invalidate()
            
        
        
class FitsFileBrowser(QtGui.QMainWindow, Ui_FitsFileBrowser):

    def __init__(self,cwd=None):
        
        super(FitsFileBrowser, self).__init__(None) 
        self.setupUi(self)
        self.show()
        self.dirModel = QtGui.QFileSystemModel()         
        self.dirModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.Dirs)
        self.dirView.setModel(self.dirModel)       
        self.dirView.hideColumn(1)
        self.dirView.hideColumn(2)      
        self.dirView.hideColumn(3)         
        select=self.dirView.selectionModel()        
        select.selectionChanged.connect(self.dirChanged)      
        
        self.splitter.setSizes([50,300])        
        self.fileView.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows);
        self.fileView.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.fileView.horizontalHeader().setResizeMode( 0, QtGui.QHeaderView.Stretch )
        self.fileView.horizontalHeader().setStyleSheet("::section{background:rgb(250,250,250);}")
        self.fileView.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.filterModel=filterProxyModel(False)
        
        self.dirModel.setRootPath(cwd)
        self.dirView.setCurrentIndex(self.dirModel.index(cwd))
        self.dirView.scrollTo(self.dirView.selectionModel().selectedIndexes()[0])
        self.dirChanged(None)
        
        self.actionFilter.triggered.connect(self.isFilterOn)
        self.actionEditFilter.triggered.connect(self.editFilter)   
        self.actionViewHeader.triggered.connect(self.viewHeader)   
        self.actionOpenFile.triggered.connect(self.OpenFile)           
        
        self.fitsHeaderViewer=[None]*100
        


    def isFilterOn(self):
        if self.actionFilter.isChecked():
            self.filterModel.filterData(True)
        else:
            self.filterModel.filterData(False) 
            
    def editFilter(self):
        w=EditfileViewerKeywordsDialog(keywords)
        w.exec_()
        
    def viewHeader(self):
        for index in self.fileView.selectionModel().selectedIndexes():
            if index.column()==0:
                try:
                    fileName ="{0}/{1}".format(self.dirModel.rootPath(),index.data())
                    header=fits.getheader(fileName)         
                    w=headerviewer.fitsHeaderViewer(header)
                    w.show()   
                    self.fitsHeaderViewer[index.row()]=w
                except:
                    pass
        
        
       
    def OpenFile(self):
       for index in self.fileView.selectionModel().selectedIndexes():
          if index.column()==0:
              if self.fileModel.isMatisse[index.row()]:
                  fileName ="{0}/{1}".format(self.dirModel.rootPath(),index.data())
                  programName = "C:/fv/bin/fv.exe"
                  #programName ="notepad.exe"
                  subprocess.Popen([programName,fileName])
       
       
       
       
       
    def dirChanged(self,selec):
        
        if selec:
            dirTxt=self.dirModel.filePath(selec.indexes()[0])
            #print dirTxt
            self.dirModel.setRootPath(dirTxt)
        dirTxt=self.dirModel.rootPath()
        self.fileModel=fileModel(dirTxt)
        self.filterModel.setSourceModel(self.fileModel)
        self.fileView.setModel(self.filterModel)







keywords=[]
keywords.append(fileViewerKeyword(function="matisseType",source="header",name="DoCatg"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET CHIP NAME",name="Detector"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET NDIT",name="NDIT"))
keywords.append(fileViewerKeyword(headerkeyword="HIERARCH ESO DET SEQ1 DIT",name="DIT"))
keywords.append(fileViewerKeyword(headerkeyword=["HIERARCH ESO INS PIL NAME","HIERARCH ESO INS PIN NAME"],checkheader="HIERARCH ESO DET CHIP NAME",checkheader_cases=["HAWAII-2RG","AQUARIUS"],name="Mode"))
keywords.append(fileViewerKeyword(headerkeyword=["HIERARCH ESO INS DIL NAME","HIERARCH ESO INS DIN NAME"],checkheader="HIERARCH ESO DET CHIP NAME",checkheader_cases=["HAWAII-2RG","AQUARIUS"],name="Resolution"))




if __name__ == "__main__":  
    try:
        app = QtGui.QApplication(sys.argv) 
        
    except:
        app = QtGui.QApplication.instance()
    
    os.chdir("D:\\Documents\\Travail\\MATISSE\\datatest")
    cwd = os.getcwd()
    w=FitsFileBrowser(cwd)
    w.show()
   
    app.exec_()