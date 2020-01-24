# -*- coding: utf-8 -*-
"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Fri May 11 05:35:18 2018

Acquisition display

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



import os,time,stat
from datetime import datetime
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import ndimage
from mpl_toolkits import mplot3d
import fnmatch
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
import matplotlib.patches as patches
from scipy import ndimage
import scipy.optimize as opt




def gaussian(center_x, center_y, width_x, width_y,height,background):
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp( -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)+background

def file_is_empty(path):
    return os.stat(path).st_size==0

def mat_showAcq(filename):
    sky=[]
    target=[]
    undy=[]
    compteursky=0
    compteurtarget=0
    compteurundy=0

    hdu=fits.open(filename)
    TorS=hdu['IMAGING_DATA'].data['TARTYP']
    detecteur=hdu[0].header['HIERARCH ESO DET CHIP NAME']
    print(detecteur)
    a=1
    b=1
    if detecteur == 'HAWAII-2RG':
        DLx=19.68
        DLy=4.92
        tolx=4
        toly=1

    else:
        DLx=31.5
        DLy=7.87
        tolx=7
        toly=1.5
    for j in range(np.size(TorS)):
        if TorS[j] == 'S':
            sky=np.append(sky,int(j))
            compteursky=compteursky+1
        elif TorS[j] == 'T':
            target=np.append(target,int(j))
            compteurtarget=compteurtarget+1
        else:
            undy=np.append(undy,int(j))
            compteurundy=compteurundy+1

    sky=list(sky)
    target=list(target)
    undy=list(undy)
    sky=[int(sky[m]) for m in range(compteursky) ]
    target=[int(target[mm]) for mm in range(compteurtarget)]
    undy=[int(undy[mmm]) for mmm in range(compteurundy)]
    csky=np.where(np.diff(sky)!=1)
    ctarget=np.where(np.diff(target)!=1)
    csky=list(np.append(-1,csky[0]))
    ctarget=list(np.append(-1,ctarget[0]))


    if len(csky)!=len(ctarget):
        del sky[-1:]
    if ctarget[-1:]!=target[-1:]:
        ctarget=np.append(ctarget,len(np.diff(target)))



    if detecteur == 'HAWAII-2RG':
        m9x=(75.10-1)*a
        m9y=(11.85-1)*b
        m10x=(74.81-1)*a
        m10y=(11.51-1)*b
        m12x=(72.88-1)*a
        m12y=(11.55-1)*b
        m13x=(72.36-1)*a
        m13y=(11.62-1)*b
    else:
        m9x=(42.28-1)*a
        m9y=(109.07-1)*b
        m10x=(36.31-1)*a
        m10y=(109.06-1)*b
        m12x=(43.19-1)*a
        m12y=(108.99-1)*b
        m13x=(38.26-1)*a
        m13y=(109.02-1)*b




    for i in range(len(csky)-1):

        vecsky=sky[csky[i]+1:csky[i+1]+1]
        vectarget=target[ctarget[i]+1:ctarget[i+1]+1]
        frame09=hdu['IMAGING_DATA'].data['DATA9'][vectarget,:,:].astype(float)
        frame10=hdu['IMAGING_DATA'].data['DATA10'][vectarget,:,:].astype(float)
        frame12=hdu['IMAGING_DATA'].data['DATA12'][vectarget,:,:].astype(float)
        frame13=hdu['IMAGING_DATA'].data['DATA13'][vectarget,:,:].astype(float)
        sky09=hdu['IMAGING_DATA'].data['DATA9'][vecsky,:,:].astype(float)
        sky10=hdu['IMAGING_DATA'].data['DATA10'][vecsky,:,:].astype(float)
        sky12=hdu['IMAGING_DATA'].data['DATA12'][vecsky,:,:].astype(float)
        sky13=hdu['IMAGING_DATA'].data['DATA13'][vecsky,:,:].astype(float)


        img9=np.mean(frame09,axis=0)-np.mean(sky09,axis=0)
        img10=np.mean(frame10,axis=0)-np.mean(sky10,axis=0)
        img12=np.mean(frame12,axis=0)-np.mean(sky12,axis=0)
        img13=np.mean(frame13,axis=0)-np.mean(sky13,axis=0)


        if detecteur=='HAWAII-2RG':
            params9 = (11,72,1,4,img9.max(),1e4)
            params10 = (11,72,1,4,img10.max(),1e4)
            params12 = (12,72,1,4,img12.max(),1e4)
            params13 = (11,72,1,4,img13.max(),1e4)
        else:
            params9 = (109,40,1,4,img9.max(),1e4)
            params10 = (109,40,1,4,img10.max(),1e4)
            params12 = (109,40,1,4,img12.max(),1e4)
            params13 = (109,40,1,4,img13.max(),1e4)

        errorfunction9 = lambda p: np.ravel(gaussian(*p)(*np.indices(img9.shape)) - img9)
        params9, success9 = opt.leastsq(errorfunction9, params9)
        g9 = gaussian(*params9)
        im9=g9(*np.indices(np.shape(img9)))

        errorfunction10 = lambda p: np.ravel(gaussian(*p)(*np.indices(img10.shape)) - img10)
        params10, success10 = opt.leastsq(errorfunction10, params10)
        g10 = gaussian(*params10)
        im10=g10(*np.indices(np.shape(img10)))

        errorfunction12 = lambda p: np.ravel(gaussian(*p)(*np.indices(img12.shape)) - img12)
        params12, success12 = opt.leastsq(errorfunction12, params12)
        g12 = gaussian(*params12)
        im12=g12(*np.indices(np.shape(img12)))

        errorfunction13 = lambda p: np.ravel(gaussian(*p)(*np.indices(img13.shape)) - img13)
        params13, success13 = opt.leastsq(errorfunction13, params13)
        g13 = gaussian(*params13)
        im13=g13(*np.indices(np.shape(img13)))

        b9x=params9[1]*a
        b9y=params9[0]*b
        b10x=params10[1]*a
        b10y=params10[0]*b
        b12x=params12[1]*a
        b12y= params12[0]*b
        b13x=params13[1]*a
        b13y=params13[0]*b

        print('barycentre au cycle:'+str(i))
        print(b9x,b9y)
        print(b10x,b10y)
        print(b12x,b12y)
        print(b13x,b13y)


        plt.figure(1,figsize=(15,7))
        plt.subplot(2,2,1)

        plt.scatter(b9x*a,b9y*b ,marker='o',color="blue")
        plt.text(b9x*a,b9y*b , str(i), color="red", fontsize=12)
        axes=plt.gca()
        axes.add_artist(patches.Ellipse((m9x,m9y), DLx*a, DLy*b , 0,edgecolor = 'magenta', fill=False, zorder = 2))
        axes.add_artist(patches.Ellipse((m9x,m9y), tolx*a, toly*b , 0,edgecolor = 'orange', fill=False, zorder = 2))
        plt.xlim(m9x-DLx*a/2-2*a,m9x+DLx*a/2+2*a)
        plt.ylim(m9y-DLy*b/2-2*b ,m9y+DLy*b/2+2*b)
        plt.ylabel('IP3,T2')


        plt.subplot(2,2,2)
        plt.scatter(b10x*a,b10y*b ,marker='o',color="blue")
        plt.text(b10x*a,b10y*b , str(i), color="red", fontsize=12)
        axes=plt.gca()
        axes.add_artist(patches.Ellipse((m10x, m10y), DLx*a, DLy*b , 0,edgecolor = 'magenta', fill=False, zorder = 2))
        axes.add_artist(patches.Ellipse((m10x,m10y), tolx*a, toly*b , 0,edgecolor = 'orange', fill=False, zorder = 2))
        plt.xlim(m10x-DLx*a/2-2*a,m10x+DLx*a/2+2*a)
        plt.ylim(m10y-DLy*b/2-2*b,m10y+DLy*b/2+2*b)
        plt.ylabel(' IP1,T1')

        plt.subplot(2,2,3)
        plt.scatter(b12x*a,b12y*b ,marker='o',color="blue")
        plt.text(b12x*a,b12y*b , str(i), color="red", fontsize=12)
        axes=plt.gca()
        axes.add_artist(patches.Ellipse((m12x, m12y), DLx*a,DLy*b , 0,edgecolor = 'magenta', fill=False, zorder = 2))
        axes.add_artist(patches.Ellipse((m12x,m12y), tolx*a, toly*b , 0,edgecolor = 'orange', fill=False, zorder = 2))
        plt.xlim(m12x-DLx*a/2-2*a,m12x+DLx*a/2+2*a)
        plt.ylim(m12y-DLy*b/2-2*b,m12y+DLy*b/2+2*b)
        plt.ylabel(' IP5,T3')

        plt.subplot(2,2,4)
        plt.scatter(b13x*a,b13y*b ,marker='o',color="blue")
        plt.text(b13x*a,b13y*b , str(i), color="red", fontsize=12)
        axes=plt.gca()
        axes.add_artist(patches.Ellipse((m13x, m13y), DLx*a, DLy*b , 0,edgecolor = 'magenta', fill=False, zorder = 2))
        axes.add_artist(patches.Ellipse((m13x,m13y), tolx*a, toly*b , 0,edgecolor = 'orange', fill=False, zorder = 2))
        plt.xlim(m13x-DLx*a/2-2*a,m13x+DLx*a/2+2*a)
        plt.ylim(m13y-DLy*b/2-2*b,m13y+DLy*b/2+2*b)
        plt.ylabel(' IP7,T4')







    plt.show()

if  __name__== '__main__' :


    arg=sys.argv
    filename=arg[1]
    mat_showAcq(filename)







