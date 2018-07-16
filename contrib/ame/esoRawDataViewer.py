# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 16:54:17 2018

@author: ame
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 06 18:25:27 2016

@author: ame
"""



import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import astropy.io.fits as fits

class matisseRawDataViewer():
    def __init__(self,filename):
        self.data=fits.open(filename)
        self.nframe=self.data['IMAGING_DATA'].header['NAXIS2']
        self.nregion=self.data['IMAGING_DATA'].header['NREGION']
        self.fig=[]
        self.ax=[]
        det=self.data[0].header['HIERARCH ESO DET NAME']
        
        if (det=="MATISSE-N"):
            detsize=1024
        else:
            detsize=2048
        
        self.fig,self.ax = plt.subplots() 
        plt.subplots_adjust(left=0.0, bottom=0.1,right=1,top=1)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)  
        axcolor = 'lightgoldenrodyellow'
        axslider =  plt.axes([0.1, 0.05, 0.8, 0.05], facecolor=axcolor)
        self.slider = Slider(axslider, 'Frame', 0, self.nframe-1, valinit=0, valfmt='%0.0f')
        self.slider.on_changed(self.showFrame)
        self. showFrame(0)
             
        
        self.ax.set_ylim((0,detsize))
        self.ax.set_xlim((0,detsize)) 
        plt.show()
        
    def showFrame(self,val):
        ymin,ymax=self.ax.get_ylim()
        xmin,xmax=self.ax.get_xlim()
        self.ax.cla()
        self.ax.set_ylim((ymin,ymax))
        self.ax.set_xlim((xmin,xmax)) 
        iframe = int(self.slider.val)
        for iregion in range(self.nregion):
            left=self.data['IMAGING_DETECTOR'].data['CORNER'][iregion][0]
            right=self.data['IMAGING_DETECTOR'].data['CORNER'][iregion][0]+ self.data['IMAGING_DETECTOR'].data['NAXIS'][iregion][0]              
            top=self.data['IMAGING_DETECTOR'].data['CORNER'][iregion][1]
            bottom=self.data['IMAGING_DETECTOR'].data['CORNER'][iregion][1]+ self.data['IMAGING_DETECTOR'].data['NAXIS'][iregion][1]                  
            self.ax.imshow(self.data['IMAGING_DATA'].data['DATA{0}'.format(iregion+1)][iframe,:,:],extent=[left,right,top,bottom])

arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "matisseRawDataViewer script to visualize raw data"
    print "Usage : filename [-options]"
    print "options :"
else:   
    filename=arg[1]
           

        
            
d = matisseRawDataViewer(filename) 

    