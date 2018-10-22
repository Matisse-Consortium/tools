# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:41:41 2018

@author: fmillour
"""

import numpy as np
from   matplotlib import pyplot as plt
from   astropy.io import fits as fits
from mat_show_oifits import open_oi

def open_hdr(oi_file):
    try:
        hdu = fits.open(oi_file)
    except IOError:
        print(("Unable to read fits file: " + oi_file))
        return {} 

    hdr = hdu[0].header
    
    return hdr

def subangle(angle1, angle2):
    delta_theta = (angle1 > angle2) * (360 - angle2 - angle1) + (angle2 > angle1) * (angle2 - angle1);
    return delta_theta

files = [u'C:/WR104/2018-07-17T064751_WR104_IR-N_IN.fits', u'C:/WR104/2018-07-17T064751_WR104_IR-LM_IN.fits']
files = [u'C:/WR104/2018-07-17T073446_94Aqr_IR-LM_OUT.fits', u'C:/WR104/2018-07-17T081931_94Aqr_IR-LM_OUT.fits', u'C:/WR104/2018-07-17T085626_94Aqr_IR-LM_OUT.fits', u'C:/WR104/2018-07-17T092831_94Aqr_IR-LM_OUT.fits', u'C:/WR104/2018-07-17T095736_94Aqr_IR-LM_OUT.fits']

BX = []
BY = []
BXA = []
BYA = []

plt.figure(1)
plt.clf();
plt.subplot(2,2,1)

for file in files:
    res = open_hdr(file)

    base=0;
    for i in np.arange(1,5):
        for j in np.arange(i+1,5):
            tel1 = i;
            tel2 = j;
            base+=1;
            
            blen = (res['HIERARCH ESO ISS PBL'+str(tel1)+str(tel2)+' START']+res['HIERARCH ESO ISS PBL'+str(tel1)+str(tel2)+' END'])/2
            bang = (res['HIERARCH ESO ISS PBLA'+str(tel1)+str(tel2)+' START']+res['HIERARCH ESO ISS PBLA'+str(tel1)+str(tel2)+' END'])/2
            
            bx = blen * np.sin(bang * np.pi / 180.)
            by = blen * np.cos(bang * np.pi / 180.)
            
            BX = np.append(BX, bx)
            BY = np.append(BY, by)
            
            BXA = np.append(BXA, bx)
            BYA = np.append(BYA, by)
            BXA = np.append(BXA, -bx)
            BYA = np.append(BYA, -by)
                
            print(base,i,j,tel1,tel2,blen, bang,bx,by)
            
            plt.axis('equal')
            plt.plot(bx,by, marker='o', markersize=6, color="red")
            plt.plot(-bx,-by, marker='o', markersize=3, color="red")
       
plt.plot(0,0, marker='+', markersize=3, color="black")     
        
BLDST = [];
BLANG = [];
for idx1,bas1 in enumerate(BX):
    for idx2 in np.arange(idx1,len(BX)):
        #if idx1 != idx2:
            #print(idx1,idx2)
        BLDST = np.append(BLDST, np.linalg.norm([BX[idx1]-BX[idx2],BY[idx1]-BY[idx2]]))
        BLANG = np.append(BLANG, 180/np.pi*np.arctan2(BX[idx1]-BX[idx2],BY[idx1]-BY[idx2]))
           #        print(bdist)

plt.subplot(2,2,3)
n, bins, patches = plt.hist(x=np.linalg.norm([BX,BY],axis=0), bins=int(np.max(BX))+1, range=[0,int(np.max(BX))+1], color='#0504aa',
                            alpha=0.7, rwidth=0.85)

plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value (m)')
plt.ylabel('Frequency')
plt.title('Histogram of bases')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 1)



plt.subplot(2,2,4)
n, bins, patches = plt.hist(x=180/np.pi*np.arctan2(BXA,BYA), bins=90, color='#0504aa',
                            alpha=0.7, rwidth=0.85)

plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value (degrees)')
plt.ylabel('Frequency')
plt.title('Histogram of angles')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 1)