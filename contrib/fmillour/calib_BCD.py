# -*- coding: utf-8 -*-
"""
Created on Apr 11 2019

@author: fmillour
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import os
import numpy as np
from astropy.io import fits
import fnmatch
from numpy import asarray as ar

direc = "/home/fmillour/fmillour/DATA_WRsurvey/reduceData/reduceData_OIFITS/"

'''
inin   = direc+"2019-04-08T025515_WR48a_IR-LM_IN_IN.fits"
inout  = direc+"2019-04-08T025515_WR48a_IR-LM_IN_OUT.fits"
outin  = direc+"2019-04-08T025515_WR48a_IR-LM_OUT_IN.fits"
outout = direc+"2019-04-08T025515_WR48a_IR-LM_OUT_OUT.fits"

'''
inin   = direc+"2019-04-08T022901_tet_Cen_IR-LM_IN_IN.fits"
inout  = direc+"2019-04-08T022901_tet_Cen_IR-LM_IN_OUT.fits"
outin  = direc+"2019-04-08T022901_tet_Cen_IR-LM_OUT_IN.fits"
outout = direc+"2019-04-08T022901_tet_Cen_IR-LM_OUT_OUT.fits"

'''
direc = "/data/users/fmillour/DATA_March/2019-03-24_OIFITS/"
inin   = direc+"2019-03-24T100706_V921_Sco_IR-LM_IN_IN.fits"
inout  = direc+"2019-03-24T100706_V921_Sco_IR-LM_IN_OUT.fits"
outin  = direc+"2019-03-24T100706_V921_Sco_IR-LM_OUT_IN.fits"
outout = direc+"2019-03-24T100706_V921_Sco_IR-LM_OUT_OUT.fits"
'''




dinin = fits.open(inin)
dinout = fits.open(inout)
doutin = fits.open(outin)
doutout = fits.open(outout)

iit3p =   dinin['OI_T3'].data['T3PHI'];
iot3p =  dinout['OI_T3'].data['T3PHI'];
oit3p =  doutin['OI_T3'].data['T3PHI'];
oot3p = doutout['OI_T3'].data['T3PHI'];

iidp =   dinin['OI_VIS'].data['VISPHI'];
iodp =  dinout['OI_VIS'].data['VISPHI'];
oidp =  doutin['OI_VIS'].data['VISPHI'];
oodp = doutout['OI_VIS'].data['VISPHI'];

print(np.shape(iit3p))
print(np.shape(iot3p))
print(np.shape(oit3p))
print(np.shape(oot3p))
print(np.arange(6))

lim = 3;

#CP 123

plt.figure (10)
plt.plot(iit3p[0,:])
plt.plot(-oit3p[0,:])
plt.plot(iot3p[3,:])
plt.plot(-oot3p[3,:])
plt.ylim(-lim,lim)


plt.figure (0)
#clos1 = iit3p[0,:]-oit3p[0,:]
#clos2 = iot3p[1,:]-oot3p[1,:]
clos1 = (iit3p[0,:]-oit3p[0,:])/2.
clos2 = (iot3p[3,:]-oot3p[3,:])/2.
plt.plot(iit3p[0,:])
plt.plot(clos1)
plt.plot(clos2)
plt.plot((clos2+clos1)/2.)
plt.ylim(-lim,lim)

#CP 124

plt.figure (11)
plt.plot(iit3p[1,:])
plt.plot(oit3p[2,:])
plt.plot(-iot3p[1,:])
plt.plot(-oot3p[2,:])
plt.ylim(-lim,lim)


plt.figure (1)
#clos1 = iit3p[0,:]-oit3p[0,:]
#clos2 = iot3p[1,:]-oot3p[1,:]
clos1 = (iit3p[1,:]+oit3p[2,:])/2.
clos2 = (-iot3p[1,:]-oot3p[2,:])/2.
plt.plot(iit3p[1,:])
plt.plot(clos1)
plt.plot(clos2)
plt.plot((clos2+clos1)/2.)
plt.ylim(-lim,lim)

#CP 134

plt.figure (12)
plt.plot(iit3p[2,:])
plt.plot(oit3p[1,:])
plt.plot(-iot3p[2,:])
plt.plot(-oot3p[1,:])
plt.ylim(-lim,lim)


plt.figure (2)
#clos1 = iit3p[0,:]-oit3p[0,:]
#clos2 = iot3p[1,:]-oot3p[1,:]
clos1 = (iit3p[2,:]+oit3p[1,:])/2.
clos2 = (-iot3p[2,:]-oot3p[1,:])/2.
plt.plot(iit3p[2,:])
plt.plot(clos1)
plt.plot(clos2)
plt.plot((clos2+clos1)/2.)
plt.ylim(-lim,lim)

#CP 234

plt.figure (13)
plt.plot(iit3p[3,:])
plt.plot(oit3p[3,:])
plt.plot(-iot3p[0,:])
plt.plot(-oot3p[0,:])
plt.ylim(-lim,lim)


plt.figure (3)
#clos1 = iit3p[0,:]-oit3p[0,:]
#clos2 = iot3p[1,:]-oot3p[1,:]
clos1 = (iit3p[3,:]-oit3p[3,:])/2.
clos2 = (iot3p[0,:]-oot3p[0,:])/2.
plt.plot(iit3p[3,:])
plt.plot(clos1)
plt.plot(clos2)
plt.plot((clos2+clos1)/2.)
plt.ylim(-lim,lim)

plt.show()
