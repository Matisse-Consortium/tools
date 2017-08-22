# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 09:09:57 2017

@author: pbe
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit


def gaus(x,a,b,x0,sigma):
    return b+a*np.exp(-(x-x0)**2/(2*sigma**2))


coefHighL=[-4.03534524,54.94605543,-277.68461223,926.20906903]
c=np.zeros(3,dtype=np.float)
hdu = fits.open('CORRFLUX_CAL_MATISSE_GEN_LAMP_L101_0001_HIGH+M.fits')
ftreal=hdu['OBJ_CORR_FLUX'].data['CORRFLUXREAL1']
ftimag=hdu['OBJ_CORR_FLUX'].data['CORRFLUXIMAG1']
corner=hdu['IMAGING_DETECTOR'].data['CORNER'][0][1]
naxis=hdu['IMAGING_DETECTOR'].data['NAXIS'][0]
hdu.close()

dsp=(ftreal**2+ftimag**2).mean(axis=0)
dsp/=np.max(dsp)
boundary=[408,448]


#iBase=0
#iWlen=1200
#n=boundary[iBase][1]-boundary[iBase][0]
#x=np.arange(n)+boundary[iBase][0]
#y=dsp[iWlen,boundary[iBase][0]:boundary[iBase][1]]
#mean = boundary[iBase][0]+10                  
#sigma = 3
#print mean,sigma
#popt,pcov = curve_fit(gaus,x,y,p0=[0.1,0.,mean,sigma])
#print popt
#plt.plot(x,y)
#plt.plot(x,gaus(x,*popt))

longu=np.zeros((naxis[1]),dtype=np.float)
#coef=np.zeros(4,dtype=np.float)
#coef[0]=-4.03534524
#coef[1]=54.94605543
#coef[2]-277.68461223
for iWlen in range(naxis[1]):
    coef=coefHighL
    n=boundary[1]-boundary[0]
    x=np.arange(n)+boundary[0]
    y=dsp[iWlen,boundary[0]:boundary[1]]
    try: 
        popt,pcov = curve_fit(gaus,x,y,p0=[1.,0.,boundary[0]+10. ,3.])
        coef[3]=926.20906903-popt[2]
        
        longu[iWlen]=np.float(np.roots(coef)[2])
    except RuntimeError:
        longu[iWlen]=-1

longu_save=np.copy(longu)

x=np.arange(naxis[1])+corner
coeflamb=np.polyfit(x[650:1450],longu[650:1450],1)
p=np.poly1d(coeflamb)
residu=longu-p(x)
resstd=np.std(residu[950:1050])
print resstd

for iWlen in range(naxis[1]):
    ecart=np.abs(longu[iWlen]-p(x[iWlen]))
    if ecart > 3*resstd :
        longu[iWlen]=p(x[iWlen])
        
plt.figure(1)
plt.plot(longu_save)
coeflamb=np.polyfit(x[650:1450],longu[650:1450],1)
p=np.poly1d(coeflamb)

plt.plot(p(x))
