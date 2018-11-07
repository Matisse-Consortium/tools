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
pathfile='/data/RawDataMatisse/2018-09-23/'
#filename='MATIS.2018-07-17T04:02:50.339.fits'
files=['MATIS.2018-09-24T05:38:45.341.fits','MATIS.2018-09-24T05:40:35.533.fits','MATIS.2018-09-24T05:42:10.713.fits','MATIS.2018-09-24T05:43:45.943.fits','MATIS.2018-09-24T05:45:30.826.fits','MATIS.2018-09-24T05:47:05.597.fits','MATIS.2018-09-24T05:48:40.828.fits','MATIS.2018-09-24T05:50:16.927.fits']

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


print('pathf',pathF)
print('chopper',chopper)
print('nonchopper',nonchopper)


cmd='rm -r '+direct+'pathAndfolderBrut.npz'
os.system(cmd)

pathfile=pathF
#filename=files[1]
for i in range (len(chopper)):
    filename=chopper[i]
    print(filename)
    tartine=[]
    hdu=fits.open(pathfile+filename)
    frames=hdu['IMAGING_DATA'].data['DATA11'].astype(float)
    #frek=hdu[0].header['HIERARCH ESO ISS CHOP FREQ']
    bcd=hdu[0].header['HIERARCH ESO INS BCD1 NAME']
    start=hdu[0].header['HIERARCH ESO TPL START ']
    starname=hdu[0].header['HIERARCH ESO OBS TARG NAME']
    TorS=hdu['IMAGING_DATA'].data['TARTYP']
    hdu.close()
    sumframes0=np.mean(frames,axis=(1,2))
   
    
    print('I just closed the file')
    
    
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
    


    
