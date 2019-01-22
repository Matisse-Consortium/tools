# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 14:32:33 2018

@author: ame
"""
import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_image_chopping script to visualize acquisition Image with choppping"
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
        dx=13
        dy=3
       
    else:
        xc=73
        yc=11
        dx=13
        dy=3
    
    
    dx2=20
    nadd=0
    nrem=0
    l=np.shape(phot)[0]
    imageadd=np.zeros([np.shape(phot)[1],np.shape(phot)[2]])
    imagerem=np.zeros([np.shape(phot)[1],np.shape(phot)[2]])
    print "Integrating the {0} images".format(l)
    for i in range(l):
      if tt[i]==0:
          nrem+=1
          imagerem+=phot[i,:,:]
      elif tt[i]==1:
          imageadd+=phot[i,:,:]
          nadd+=1
    
    image=imageadd/nadd-imagerem/nrem
    plt.imshow(image)
    
   
    
    
    def g2d ((x,y), amplitude, xo, yo, sigma_x, sigma_y, offset):
        xo = float(xo)
        yo = float(yo)    
        a = 1./(2*sigma_x**2)   
        b = 1./(2*sigma_y**2)
        res=offset + amplitude*np.exp( - (a*((x-xo)**2) + b*((y-yo)**2)))
        return res.ravel()
      
    initial_guess = (3,xc,yc,dx/2.35,dy/2.35,0)
    x = np.linspace(0, np.shape(image)[1]-1, np.shape(image)[1])
    y = np.linspace(0, np.shape(image)[0]-1, np.shape(image)[0])
    x,y = np.meshgrid(x, y)
    
    print np.shape(image)
    print np.shape(x)
    print np.shape(y)
    popt, pcov = opt.curve_fit(g2d, (x,y), image.reshape( np.shape(image)[1]*np.shape(image)[0]), p0 = initial_guess)
    xMes=popt[1]
    yMes=popt[2] 
    wMes=popt[3]*2.35
    hMes=popt[4]*2.35
    mMes=popt[0]    
    print "min={0} max={1} nobj={2} nsky={3}".format(np.min(image),np.max(image),nadd,nrem)
    plt.scatter(xMes,yMes,marker="*",s=15,color='red')
    print "xc={0} yc={1} w={2} h={3} ampl={4}".format(xMes,yMes,wMes,hMes,mMes)
    plt.show()
