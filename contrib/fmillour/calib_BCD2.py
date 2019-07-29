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
from shutil import copyfile



lim = 180;

direc = "/mnt/c/DATA/DATA_WRsurvey/reduceData_OIFITS/"

'''
inin   = direc+"2019-04-08T025515_WR48a_IR-LM_IN_IN.fits"
inout  = direc+"2019-04-08T025515_WR48a_IR-LM_IN_OUT.fits"
outin  = direc+"2019-04-08T025515_WR48a_IR-LM_OUT_IN.fits"
outout = direc+"2019-04-08T025515_WR48a_IR-LM_OUT_OUT.fits"
'''

'''


'''
inin   = direc+"2019-04-08T022901_tet_Cen_IR-LM_IN_IN.fits"
inout  = direc+"2019-04-08T022901_tet_Cen_IR-LM_IN_OUT.fits"
outin  = direc+"2019-04-08T022901_tet_Cen_IR-LM_OUT_IN.fits"
outout = direc+"2019-04-08T022901_tet_Cen_IR-LM_OUT_OUT.fits"

direc = "/mnt/c/DATA/2019-03-24_OIFITS/"


'''
inin   = direc+"2019-03-24T100706_V921_Sco_IR-LM_IN_IN.fits"
inout  = direc+"2019-03-24T100706_V921_Sco_IR-LM_IN_OUT.fits"
outin  = direc+"2019-03-24T100706_V921_Sco_IR-LM_OUT_IN.fits"
outout = direc+"2019-03-24T100706_V921_Sco_IR-LM_OUT_OUT.fits"
'''


inin   = direc+"2019-03-24T033720_HD_87643_IR-LM_IN_IN.fits"
inout  = direc+"2019-03-24T033720_HD_87643_IR-LM_IN_OUT.fits"
outin  = direc+"2019-03-24T033720_HD_87643_IR-LM_OUT_IN.fits"
outout = direc+"2019-03-24T033720_HD_87643_IR-LM_OUT_OUT.fits"

def calib_BCD(iifile, iofile, oifile, oofile, outputfile=os.getenv("HOME")+"/toto.fits",lim=180,plot=1):
    
    copyfile(iofile, outputfile)
    outhdu  = fits.open(outputfile, mode='update')

    dinin   = fits.open(iifile)
    dinout  = fits.open(iofile)
    doutin  = fits.open(oifile)
    doutout = fits.open(oofile)

    iitwl =   dinin['OI_WAVELENGTH'].data['EFF_WAVE'];
    
    iit3p =   dinin['OI_T3'].data['T3PHI'];
    iot3p =  dinout['OI_T3'].data['T3PHI'];
    oit3p =  doutin['OI_T3'].data['T3PHI'];
    oot3p = doutout['OI_T3'].data['T3PHI'];

    iidp =   dinin['OI_VIS'].data['VISPHI'];
    iodp =  dinout['OI_VIS'].data['VISPHI'];
    oidp =  doutin['OI_VIS'].data['VISPHI'];
    oodp = doutout['OI_VIS'].data['VISPHI'];

    iiv2 =   dinin['OI_VIS2'].data['VIS2DATA'];
    iov2 =  dinout['OI_VIS2'].data['VIS2DATA'];
    oiv2 =  doutin['OI_VIS2'].data['VIS2DATA'];
    oov2 = doutout['OI_VIS2'].data['VIS2DATA'];


                
    # First plot all original observables
    if plot:
        plt.figure (50)
        for i in np.arange(6):
            plt.subplot(3,2,i+1);    
            plt.plot(iitwl*1e6, iiv2[i,:])
        plt.figure (51)
        for i in np.arange(4):
            plt.subplot(2,2,i+1);    
            plt.plot(iitwl*1e6, iit3p[i,:])
        plt.figure (52)
        for i in np.arange(6):
            plt.subplot(3,2,i+1);    
            plt.plot(iitwl*1e6, iidp[i,:])  

    nwlen = np.shape(iit3p)[1];
    # addup the different exposures
    nrepeatii = int(np.shape(iit3p)[0]/4);
    print(nrepeatii)
    nrepeatii = int(np.shape(iiv2)[0]/6);
    nrepeatoi = int(np.shape(oiv2)[0]/6);
    nrepeatio = int(np.shape(iov2)[0]/6);
    nrepeatoo = int(np.shape(oov2)[0]/6);
    print(nrepeatii)
    print(nrepeatio)
    print(nrepeatoi)
    print(nrepeatoo)

    # Store multiple exposures data into the first 6 rows
    if nrepeatii > 1:
        for i in np.arange(nrepeatii-1):
            for j in np.arange(6):
                iiv2[j,:] += iiv2[(i+1)*6+j,:]
                iidp[j,:] += iidp[(i+1)*6+j,:]
            for j in np.arange(4):
                iit3p[j,:] += iit3p[(i+1)*4+j,:]
        
    if nrepeatio > 1:
        for i in np.arange(nrepeatio-1):
            for j in np.arange(6):
                iov2[j,:] += iov2[(i+1)*6+j,:]
                iodp[j,:] += iodp[(i+1)*6+j,:]
            for j in np.arange(4):
                iot3p[j,:] += iot3p[(i+1)*4+j,:]
        
    if nrepeatoi > 1:
        for i in np.arange(nrepeatoi-1):
            for j in np.arange(6):
                oiv2[j,:] += oiv2[(i+1)*6+j,:]
                oidp[j,:] += oidp[(i+1)*6+j,:]
            for j in np.arange(4):
                oit3p[j,:] += oit3p[(i+1)*4+j,:]
        
    if nrepeatoo > 1:
        for i in np.arange(nrepeatoo-1):
            for j in np.arange(6):
                oov2[j,:] += oov2[(i+1)*6+j,:]
                oodp[j,:] += oodp[(i+1)*6+j,:]
            for j in np.arange(4):
                oot3p[j,:] += oot3p[(i+1)*4+j,:]

#Treat closure phases
    idx  = np.array([[0,0, 3, 3],[1,2,1,2],[2,1,2,1],[3,3,0,0]])
    sign = np.array([[1,-1,1,-1],[1,1,-1,-1],[1,1,-1,-1],[1,-1,1,-1]])

#Initialize closure phase with same shape as input
    closfinal = np.zeros((4,nwlen));
    plt.figure (51)
    
    for i in np.arange(4):
        
        closfinal[i,:] = (sign[i,0] * iit3p[idx[i,0],:] + sign[i,1] * oit3p[idx[i,1],:] +\
                          sign[i,2] * iot3p[idx[i,2],:] + sign[i,3] * oot3p[idx[i,3],:])/\
                          (nrepeatii+nrepeatoi+nrepeatio+nrepeatoo)
        plt.subplot(2,2,i+1);
        plt.plot(iitwl*1e6, closfinal[i,:])
##        for j in np.arange(nrepeatii):
##            plt.plot(iitwl*1e6, sign[i,0] * iit3p[idx[i,0]+j*4,:])
##        for j in np.arange(nrepeatoi):
##            plt.plot(iitwl*1e6, sign[i,1] * oit3p[idx[i,1]+j*4,:])
##        for j in np.arange(nrepeatio):
##            plt.plot(iitwl*1e6, sign[i,2] * iot3p[idx[i,2]+j*4,:])
##        for j in np.arange(nrepeatoo):
##            plt.plot(iitwl*1e6, sign[i,3] * oot3p[idx[i,3]+j*4,:])
        plt.ylim(-lim,lim)

#Treat differential phases
    idx  = np.array([[0,0, 0, 0],[1,1,1,1],[2,3,4,5],[3,2,5,4],[4,5,2,3],[5,4,3,2]])
    sign = np.array([[1,-1,1,-1],[1,1,-1,-1],[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])

    dpfinal = np.zeros((6,nwlen));
    plt.figure (52)
    print(np.shape(idx))
    for i in np.arange(6):
        dpfinal[i,:] = (sign[i,0] * iidp[idx[i,0],:] + sign[i,1] * oidp[idx[i,1],:] +\
                        sign[i,2] * iodp[idx[i,2],:] + sign[i,3] * oodp[idx[i,3],:])/\
                        (nrepeatii+nrepeatoi+nrepeatio+nrepeatoo)
        plt.subplot(3,2,i+1); 
        plt.plot(iitwl*1e6, dpfinal[i,:])
        plt.ylim(-lim,lim)

#Treat visibilities
    v2final = np.zeros((6,nwlen));
    
    plt.figure (50)
    print(np.shape(idx))

    
        
    for i in np.arange(6):
        plt.subplot(3,2,i+1);  
        v2final[i,:] = (iiv2[i,:] + oiv2[idx[i,1],:] + iov2[idx[i,2],:] + oov2[idx[i,3],:])\
                       /(nrepeatii+nrepeatoi+nrepeatio+nrepeatoo)
        plt.plot(iitwl*1e6, v2final[i,:])
##        for j in np.arange(nrepeatii):
##            plt.plot(iitwl*1e6, sign[i,0] * iiv2[idx[i,0]+j*6,:])
##        for j in np.arange(nrepeatoi):
##            plt.plot(iitwl*1e6, sign[i,1] * oiv2[idx[i,1]+j*6,:])
##        for j in np.arange(nrepeatio):
##            plt.plot(iitwl*1e6, sign[i,2] * iov2[idx[i,2]+j*6,:])
##        for j in np.arange(nrepeatoo):
##            plt.plot(iitwl*1e6, sign[i,3] * oov2[idx[i,3]+j*6,:])
        plt.ylim(0,1)
        
    outhdu['OI_T3'].data['T3PHI']      = closfinal;
    outhdu['OI_VIS'].data['VISPHI']    = dpfinal;
    outhdu['OI_VIS2'].data['VIS2DATA'] = v2final;
    
    outhdu.flush()  # changes are written back to original.fits
    outhdu.close()

calib_BCD(inin, inout, outin, outout, outputfile=os.getenv("HOME")+"/toto.fits")

    
print(os.getenv("HOME")+"/toto.fits")


    
plt.show()
