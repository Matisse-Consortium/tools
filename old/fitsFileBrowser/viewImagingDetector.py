# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 11:49:31 2016

@author: ame
"""


from matplotlib import pyplot as plt
from matplotlib import patches
from astropy.io import fits

class detectorRegion:
    def __init__(self,x0,y0,dx,dy,name):
        self.x0=x0
        self.y0=y0
        self.x1=x0+dx
        self.y1=y0+dy
        self.dx=dx
        self.dy=dy
        self.name=name
        self.xc=x0+dx/2
        self.yc=y0+dy/2
        
    def __str__(self):
        return "size {0}x{1}  : pixels ({2},{3})->({4},{5})".format(self.dx,self.dy,self.x0,self.y0,self.x1,self.y1)


if __name__ == "__main__":  
    dir="D:/Documents/Travail/MATISSE/testplan/2016-12-20 franges/"
    filename=dir+"LIGHT_BSN1234.fits"
    
    data=fits.getdata(filename,'IMAGING_DETECTOR')
    
    nregion=len(data)
    
    regions=[]
    
    fig=plt.figure()
    ax = fig.add_subplot(111)
    plt.ylim((0,1024))
    plt.xlim((0,1024))    
    for datai in data:       
        region=detectorRegion(datai.field('CORNER')[0],datai.field('CORNER')[1],datai.field('NAXIS')[0],datai.field('NAXIS')[1],datai.field('REGNAME').split('\x00')[0])
        regions.append(region)
        p=patches.Rectangle((region.x0,region.y0),region.dx,region.dy)
        ax.add_patch(p)
        ax.text(region.xc,region.yc,region.name, horizontalalignment='center',verticalalignment='center',fontsize=8)
       
    