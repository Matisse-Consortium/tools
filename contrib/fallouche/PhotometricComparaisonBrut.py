# -*- coding: utf-8 -*-
"""
Created on Fri May 11 05:35:18 2018

@author: 
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

tplstart=sys.argv[1]
detec=sys.argv[2]
direct=sys.argv[3]
if direct[-1:] != '/' :
    direct=direct+'/'

cmd='python tplstartBrut.py ' +tplstart+' '+detec +' ' +direct
os.system(cmd)
file=np.load(direct+'pathAndfolderBrut.npz')

pathF=str(file['arr_0'])
nonchopper=file['arr_1']
chopper=file['arr_2']

cmd='rm -r '+direct+'pathAndfolderBrut.npz'
os.system(cmd)


moyenneC1=[]
moyenneC2=[]
moyenneC3=[]
moyenneC4=[]

moyenneNC1=[]
moyenneNC2=[]
moyenneNC3=[]
moyenneNC4=[]
thedateC=[]
thedateNC=[]
shut=[]
thebcdc=[]
thebcdnc=[]
framec=[]
framenc=[]
couleur=['blue','orange','red','green']
for i in range(len(chopper)):
    print(pathF+chopper[i])
    hduC=fits.open(pathF+chopper[i])
    dateC=hduC[0].header['DATE-OBS'][11:]
    thedateC=np.append(thedateC,dateC)
    bcd1C=hduC[0].header['HIERARCH ESO INS BCD1 NAME']
    bcd2C=hduC[0].header['HIERARCH ESO INS BCD2 NAME']
    bcdc=[bcd1C,bcd2C]
    detect=hduC[0].header['HIERARCH ESO DET CHIP NAME']
    if detect=='AQUARIUS':
        disp=hduC[0].header['HIERARCH ESO INS DIN NAME']
    if detect=='HAWAII-2RG':
        disp=hduC[0].header['HIERARCH ESO INS DIL NAME']
    varC1=hduC['IMAGING_DATA'].data['DATA9'].astype(float)
    varC2=hduC['IMAGING_DATA'].data['DATA10'].astype(float)
    varC3=hduC['IMAGING_DATA'].data['DATA12'].astype(float)
    varC4=hduC['IMAGING_DATA'].data['DATA13'].astype(float)
    target=hduC[0].header['HIERARCH ESO OBS TARG NAME']
    
    shutter1=hduC[0].header['HIERARCH ESO INS BSN1 ST']
    shutter2=hduC[0].header['HIERARCH ESO INS BSN2 ST']
    shutter3=hduC[0].header['HIERARCH ESO INS BSN3 ST']
    shutter4=hduC[0].header['HIERARCH ESO INS BSN4 ST']
    shutter=[shutter1,shutter2,shutter3,shutter4]
    shut=np.append(shut,np.where(shutter)[0][0])  
    thebcdc=np.append(thebcdc,[bcdc])
    nbframeschopper=np.shape(varC1)
    start=hduC[0].header['HIERARCH ESO TPL START']
   
    moyenneC1=np.append(moyenneC1,np.mean(varC1,axis=(1,2)))
    moyenneC2=np.append(moyenneC2,np.mean(varC2,axis=(1,2)))
    moyenneC3=np.append(moyenneC3,np.mean(varC3,axis=(1,2)))
    moyenneC4=np.append(moyenneC4,np.mean(varC4,axis=(1,2)))
    framec=np.append(framec,nbframeschopper[0])

    
            

for i in range(len(nonchopper)):
    hduNC=fits.open(pathF+'/'+nonchopper[i])
    dateNC=hduNC[0].header['DATE-OBS'][11:]
    thedateNC=np.append(thedateNC,dateNC)
    bcd1NC=hduNC[0].header['HIERARCH ESO INS BCD1 NAME']
    bcd2NC=hduNC[0].header['HIERARCH ESO INS BCD2 NAME']
    bcdNC=[bcd1NC,bcd2NC]
    detect=hduC[0].header['HIERARCH ESO DET CHIP NAME']
    if detect=='AQUARIUS':
        disp=hduC[0].header['HIERARCH ESO INS DIN NAME']
    if detect=='HAWAII-2RG':
        disp=hduC[0].header['HIERARCH ESO INS DIL NAME']
    varNC1=hduNC['IMAGING_DATA'].data['DATA9'].astype(float)
    varNC2=hduNC['IMAGING_DATA'].data['DATA10'].astype(float)
    varNC3=hduNC['IMAGING_DATA'].data['DATA12'].astype(float)
    varNC4=hduNC['IMAGING_DATA'].data['DATA13'].astype(float)
    thebcdnc=np.append(thebcdnc,[bcdNC])
    nbframesnonchopper=np.shape(varNC1)
    
    moyenneNC1=np.append(moyenneNC1,np.mean(varNC1,axis=(1,2)))
    moyenneNC2=np.append(moyenneNC2,np.mean(varNC2,axis=(1,2)))
    moyenneNC3=np.append(moyenneNC3,np.mean(varNC3,axis=(1,2)))
    moyenneNC4=np.append(moyenneNC4,np.mean(varNC4,axis=(1,2)))
    framenc=np.append(framenc,nbframesnonchopper[0])

   
        
print(shut)
   
lesmax=[np.max(moyenneC1),np.max(moyenneC2),np.max(moyenneC3),np.max(moyenneC4)]
lesmin=[np.min(moyenneC1),np.min(moyenneC2),np.min(moyenneC3),np.min(moyenneC4)]
lemin=np.min(lesmin)
lemax=np.max(lesmax)
cc1=int(np.around((np.max(moyenneC1)-np.min(moyenneC1))))
cc2=int(np.around((np.max(moyenneC2)-np.min(moyenneC2))))
cc3=int(np.around((np.max(moyenneC3)-np.min(moyenneC3))))
cc4=int(np.around((np.max(moyenneC4)-np.min(moyenneC4))))
ccc=[cc1,cc2,cc3,cc4]
cc=np.max(ccc)
c=4*cc
print(cc)
print(c)
print(lemin)
dicoSout={'BSN1':1,'BSN2':2,'BSN3':4,'BSN4':3}
dicoSin ={'BSN1':2,'BSN2':1,'BSN3':3,'BSN4':4}
s=['BSL1','BSL2','BSL4','BSL3']

thebcdc=np.reshape(thebcdc,(len(chopper),2))
thebcdnc=np.reshape(thebcdnc,(len(nonchopper),2))
plt.figure(1,figsize=(15,10))   
plt.title('donnees'+target+'_'+start+' choppees @'+disp+' in juska la '+str(np.sum(framec[0:4]))+' frames au dela bcd out  \n Integration sur toute la fenetre')
plt.plot(moyenneC1,color=couleur[0],label=('+0x'+str(cc)))
plt.plot(moyenneC2+cc,color=couleur[1],label=('+1x'+str(cc)))
plt.plot(moyenneC3+2*cc,color=couleur[2],label=('+2x'+str(cc)))
plt.plot(moyenneC4+3*cc,color=couleur[3],label=('+3x'+str(cc)))

for i in range(len(chopper)):
    plt.text(np.sum(framec[0:i]),lemin,thedateC[i])
    plt.text(np.sum(framec[0:i]),3*cc+np.max(moyenneC4),'BSN'+str(int(shut[i])+1)+'\n'+str(thebcdc[i,:]),fontsize=14)
    plt.plot(np.ones(c)*np.sum(framec[0:i]),np.arange(c)+lemin,linestyle='-.',color='black')
    axes=plt.gca()
    if thebcdc[i,0]=='OUT':
        axes.add_artist(patches.Rectangle((np.sum(framec[0:i]),cc*(dicoSout['BSN'+str(int(shut[i])+1)]-1)+lemin),(framec[i]),lemax-lemin,fill=True,color='grey',alpha=0.3,zorder = 2))
    if thebcdc[i,0]=='IN':
        print(cc*(dicoSin['BSN'+str(int(shut[i])+1)]-1))
        axes.add_artist(patches.Rectangle((np.sum(framec[0:i]),cc*(dicoSin['BSN'+str(int(shut[i])+1)]-1)+lemin),(framec[i]),lemax-lemin,fill=True,color='grey',alpha=0.3,zorder = 2))
    
for i in range(4):
    plt.text(len(moyenneC1),i*(cc)+(np.mean(moyenneC1)),s[i],rotation=90,color=couleur[i],fontsize=16)
plt.xlabel('number of frames')
plt.ylabel('averaged intensity in the photo channel per frame')
plt.legend()
#plt.savefig('/data/users/fal/RUN1C/'+target.replace(' ','')+'_Choppe_'+disp+'_'+start+'.pdf')


plt.figure(2,figsize=(15,10))
plt.title('donnees'+target+'_'+start+' non choppes @'+disp+' Integration  sur toute la fenetre')
plt.plot(moyenneNC1,color=couleur[0],label=('+0x'+str(cc)))
plt.plot(moyenneNC2+cc,color=couleur[1],label=('+1x'+str(cc)))
plt.plot(moyenneNC3+2*cc,color=couleur[2],label=('+2x'+str(cc)))
plt.plot(moyenneNC4+3*cc,color=couleur[3],label=('+3x'+str(cc)))

for i in range(len(nonchopper)):
    plt.text(np.sum(framenc[0:i]),lemin,thedateNC[i])
    plt.text(np.sum(framenc[0:i]),3*cc+np.max(moyenneNC4),'BSN'+str(int(shut[i])+1)+'\n'+str(thebcdnc[i,:]),fontsize=14)
    plt.plot(np.ones(c)*np.sum(framenc[0:i]),np.arange(c)+lemin,linestyle='-.',color='black')

for i in range(4):
    plt.text(len(moyenneNC1),i*(cc)+(np.mean(moyenneNC1)),s[i],rotation=90,color=couleur[i],fontsize=16)
plt.xlabel('number of frames')
plt.ylabel('averaged intensity in the photo channel per frame')
plt.legend()
#plt.savefig('/data/users/fal/RUN1C/'+target.replace(' ','')+'_NonChoppe_'+disp+'_'+start+'.pdf')



plt.show()

