# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 14:32:33 2018

@author: ame
"""
import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np



arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_image_chopping_plot script to visualize chopping on acquisition images"
    print "Usage : filename [-options]"
    print "options :"
    print "--skip nframe (or -s)"
    print "--beam ibeam (or -b)"
else:   
    filename=arg[1]
    
    
    
    nskip=0
    ibeam=1
    narg=len(arg)
    print arg
    for i in range(2,narg):
        if (arg[i] == '--skip' or arg[i] == '-s'):
            nskip=int(arg[i+1])
        elif (arg[i]=='--beam' or arg[i] == '-b'):
            ibeam=int(arg[i+1])
    
    iregion=['9','10','12','13']  
    print "Analyzing Chopping in Acq file={0} on beam {1} ('DATA{2}') and skipping {3} frames".format(filename,ibeam,iregion[ibeam-1],nskip)
    d=fits.open(filename)
    
    tartype=d['IMAGING_DATA'].data['TARTYP']
    tt=(tartype=='T')*1+(tartype=='U')*0.5
       
    
   
    phot=d['IMAGING_DATA'].data['DATA{0}'.format(iregion[ibeam-1])]*1.0
    print np.shape(phot)
    
    start=0
    tt=tt[nskip:]    
    phot=phot[nskip:]     
    
    
    colors=[]
    symsize=[]
    for i in range(len(tt)):
        if (tt[i]==0):
            colors.append('red')
            symsize.append(20)
        elif (tt[i]==0.5):
            colors.append('#00FF00')
            symsize.append(20)
        else:
            colors.append('blue')
            symsize.append(20)
     
    det=d[0].header['HIERARCH ESO DET NAME']
    print "det={0}".format(det)
     
    if (det=="MATISSE-N"):
        xc=40
        yc=109
        dx=25
        dy=6
       
    else:
        xc=73
        yc=11
        dx=15
        dy=5
    
    
    dx2=20
    
    photIntegr=np.sum(np.sum(phot[:,yc-dy/2-1:yc-1+dy/2,xc-dx/2-1:xc-1+dx/2],axis=1),axis=1)

    maxi=np.max(photIntegr)
    mini=np.min(photIntegr)
    if (maxi!=mini):
        photIntegr=(photIntegr-mini)/(maxi-mini)
   
    x=np.arange(len(photIntegr))
    #plt.scatter(x,tt,marker="s",color=colors,s=symsize*100)   
    plt.scatter(x,photIntegr,marker="o",s=symsize,color=colors)
    plt.ylim((0,1))


    plt.show()
