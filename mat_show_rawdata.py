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
import numpy as np
from PIL import Image
import astropy.io.fits as fits

class mat_show_rawdata():
    def __init__(self,filename):
        self.data=fits.open(filename)
        self.nframe=self.data['IMAGING_DATA'].header['NAXIS2']
        self.nregion=self.data['IMAGING_DATA'].header['NREGION']
        self.tartyp = self.data['IMAGING_DATA'].data.field('TARTYP')

 

        self.fig = plt.figure() 
        self.ax =  plt.axes([0.1, 0.26, 0.8, 0.7] ,facecolor="white")
        
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)  
        
        
        
        img = Image.new('RGB', (self.nframe, 1), (0, 255, 0))
        for i,t in enumerate(self.tartyp):
            if t=='T':
                img.putpixel((i,0), (0,255,0))
            elif t=='U':
                img.putpixel((i,0), (255,0,0))
            elif t=='S':
                img.putpixel((i,0), (0,255,255))
        
        self.tartype_axe =  plt.axes([0.1, 0.10, 0.8, 0.05] )
        plt.xticks([]), plt.yticks([])
        plt.tight_layout()
        plt.imshow(img,extent=[0,255,0,10])
        self.tartype_axe.set_ylim((0,10))
        self.tartype_axe.set_xlim((0,255)) 
        
        axslider =  plt.axes([0.1, 0.05, 0.8, 0.05] ,facecolor="lightgrey")
        self.slider = Slider(axslider, 'Frame', 0, self.nframe-1, valinit=0, valfmt='%0.0f' , color="yellow")
        self.slider.on_changed(self.showFrame)
        
        axsliderMax =  plt.axes([0.1, 0.20, 0.8, 0.02] ,facecolor="lightgrey")
        self.sliderMax = Slider(axsliderMax, 'Max Cut', 0, 50000, valinit=50000, valfmt='%0.0f',  color="skyblue")
        self.sliderMax.on_changed(self.updateFrame)
   
   
        axsliderMin =  plt.axes([0.1, 0.23, 0.8, 0.02]  ,facecolor="lightgrey")
        self.sliderMin = Slider(axsliderMin, 'Min Cut', 0, 50000, valinit=0, valfmt='%0.0f', color="deepskyblue")
        self.sliderMin.on_changed(self.updateFrame)     
        
        
        self. showFrame(True)
        plt.show()
        
    def updateFrame(self,val):
        
        self.showFrame(False)
        
    def  showFrame(self,auto = False):
        ymin,ymax=self.ax.get_ylim()
        xmin,xmax=self.ax.get_xlim()
        self.ax.cla()
        self.ax.set_ylim((ymin,ymax))
        self.ax.set_xlim((xmin,xmax)) 
        iframe = int(self.slider.val)
        maxi = self.sliderMax.val
        mini = self.sliderMin.val
        if (self.nregion % 3) == 0:
            nrow  = self.nregion / 3
            ypos = 5;
            for iy in range(3):
                xpos = 5; 
                for ix in range(nrow):
                    iregion = iy*nrow + ix+1
                   
                    img = self.data['IMAGING_DATA'].data['DATA{0}'.format(iregion)][iframe,:,:]
                    
                    left = xpos
                    right = xpos + np.shape(img)[1]
                    top = ypos 
                    bottom = ypos + np.shape(img)[0]
                    
                    if iy == 1 and ix == 1 and auto == True:
                        print ("min = {0} max = {1}".format(self.sliderMax.val,self.sliderMin.val))
                        image_histogram, bins = np.histogram(img.flatten(),
                             int(np.sqrt(len(img.flatten()))), normed=False)
                        maximg = image_histogram - np.median(image_histogram) > 15 
                        maxi = np.max((bins[1:])[maximg])
                        mini = np.median(img.flatten())
                        maxi = mini + (maxi - mini)/3
                        self.sliderMax.set_val(maxi)
                        self.sliderMin.set_val(mini)                                     
                        print ("min = {0} max = {1}".format(self.sliderMin.val,self.sliderMax.val))
                       
                    
                    if maxi<=mini:
                        maxi = mini+1
                        self.sliderMax.set_val(maxi)
                    self.ax.imshow(img,extent=[left,right,top,bottom],
                      vmin=mini,vmax=maxi, interpolation = 'nearest',
                                                        cmap = 'afmhot')
                    print ("iregion = {0} : [{1},{2},{3},{4}]".format(iregion,left,right,top,bottom))
                    
                    
                    
                   
                    xpos += np.shape(img)[1]+5
                ypos += np.shape(img)[0]+5                
            self.ax.set_ylim((0,ypos))
            self.ax.set_xlim((0,xpos)) 
                
if __name__ == '__main__':
    arg=sys.argv

    if (arg[1]=="--help" or arg[1]== "-h"):
        print "mat_show_rawdata script to visualize raw data"
        print "Usage : filename [-options]"
        print "options :"
    else:   
        filename=arg[1]
               

            
                
    d = mat_show_rawdata(filename) 

    
