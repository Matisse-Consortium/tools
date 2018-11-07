import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import ndimage
from mpl_toolkits import mplot3d
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

sumframes=[]


tplstart=sys.argv[1]    ####### in the format hh:mm:ss    ex 05:26:21 
detec=sys.argv[2]      #########  detector, HAWAI for HAWAII-2RG and AQUA for the AQUARIUS
direct=sys.argv[3]     ######### a directory where a buffering file is created under the tplstartBrut routine, then read then destroyed
if direct[-1:] != '/' :
    direct=direct+'/'

cmd='python tplstartBrut.py ' +tplstart+' '+detec +' ' +direct   ###### creating the buffer file
os.system(cmd)
file=np.load(direct+'pathAndfolderBrut.npz')   ###### reading the buffer file
pathF=str(file['arr_0'])
chopper=file['arr_2']

cmd='rm -r '+direct+'pathAndfolderBrut.npz'
os.system(cmd)    ###### destroying the buffer file

if detect=='HAWAI':
    DATAnb='DATA11'
if detect=='AQUA':
    DATAnb='DATA5'

pathfile=pathF
#filename=files[1]
for i in range (len(chopper)):
    filename=chopper[i]
    print(filename)
    tartine=[]
    hdu=fits.open(pathfile+filename)
    ################### DATA11 because we only read the chopping status in the L band files, for reading the N band use DATA5
    frames=hdu['IMAGING_DATA'].data[DATAnb].astype(float)
    ######### useless keywords in order, the chopping frequency (which can be missing), the bcd1 status, etc etc
    #frek=hdu[0].header['HIERARCH ESO ISS CHOP FREQ']
    bcd=hdu[0].header['HIERARCH ESO INS BCD1 NAME']
    start=hdu[0].header['HIERARCH ESO TPL START ']
    starname=hdu[0].header['HIERARCH ESO OBS TARG NAME']
    TorS=hdu['IMAGING_DATA'].data['TARTYP']
    hdu.close()
    
    sumframes0=np.mean(frames,axis=(1,2))
   
    
    for j in range(np.size(TorS)):
        if TorS[j] == 'S':
            tartine=np.append(tartine,0)
        elif TorS[j] == 'T':
            tartine=np.append(tartine,1)        
        else:
            tartine=np.append(tartine,0.5)
            
            
            
    plt.subplot(len(chopper),1,i+1)
    plt.plot(sumframes0/np.max(sumframes0),marker='*',label='flux',linestyle='')
    plt.plot(tartine,linestyle='',marker='.',label='on sky')
    plt.title(filename)
    plt.ylabel('tplStart='+str(start))
plt.show()
    


    
