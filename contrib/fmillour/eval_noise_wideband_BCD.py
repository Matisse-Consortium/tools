# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:44:16 2017

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

# LM band
wlen_lim  = [2.9,5.2]
wlen_use  = [3,5]
mjd_range = [57875,57876]
pattern = "*LAMP_L*.fits*"
filtype = "BCDIN"
demodtype = "false";
dirname = "/home/fmillour/MATISSE/TESTS_NICE/transfunc5_el"
V_spec = 0.075;
V_goal = 0.025;
V_expected = 0.006;
CP_spec = 40e-3*180/np.pi;
CP_goal = 1e-3*180/np.pi;
CP_expected = 8e-3*180/np.pi;
DP_spec = 30e-3*180/np.pi;
DP_goal = 1e-3*180/np.pi;
DP_expected = 11e-3*180/np.pi;
VD_spec     = 0.015
VD_goal     = 0.005
VD_expected = 0.002;


# M band
#wlen_lim  = [2.9,5.2]
#wlen_use  = [4.3,5]
#mjd_range = [57875,57876]
#pattern = "*LAMP_L*.fits*"
#filtype = "BCDIN"
#demodtype = "false";
#dirname = "/home/fmillour/MATISSE/TESTS_NICE/transfunc5_el"
#V_spec = 0.075;
#V_goal = 0.025;
#V_expected = 0.006;
#CP_spec = 40e-3*180/np.pi;
#CP_goal = 1e-3*180/np.pi;
#CP_expected = 8e-3*180/np.pi;
#DP_spec = 30e-3*180/np.pi;
#DP_goal = 1e-3*180/np.pi;
#DP_expected = 11e-3*180/np.pi;
#VD_spec     = 0.015
#VD_goal     = 0.005
#VD_expected = 0.002;

#### N band
#wlen_lim  = [7,14]
#wlen_use  = [8,13]
#mjd_range = [5785,578551.6]
#pattern = "*LAMP_N*.fits*"
#filtype = "BCDOUT"
#demodtype = "true";
##dirname = "/home/fmillour/MATISSE/TESTS_NICE/transfunc5_el"
#dirname = "/home/fmillour/MATISSE/TESTS_NICE/transfunc6"
#V_spec = 0.075;
#V_goal = 0.025;
#V_expected = 0.007;
#CP_spec = 40e-3*180/np.pi;
#CP_goal = 1e-3*180/np.pi;
#CP_expected = 8e-3*180/np.pi;
#DP_spec = 30e-3*180/np.pi;
#DP_goal = 1e-3*180/np.pi;
#DP_expected = 11e-3*180/np.pi;
#VD_spec     = 0.05
#VD_goal     = 0.02
#VD_expected = 0.002;



time_lim  = [0.5,mjd_range[1]-mjd_range[0]+0.5]

cols = ["YellowGreen", "Tomato", "SteelBlue", "Sienna", "Magenta", "Navy", "SeaGreen"]


files = os.listdir(dirname)

V2  = []
VD  = []
DP  = []
CP  = []
FL  = [];
MJD = [];
idx = wl = [];


count = 0;

for i,fil in enumerate(files):
    if fnmatch.fnmatch(fil,pattern):
        if(fnmatch.fnmatch(fil,"*"+filtype+"*")):
            col = 'b'
            
            print(fil)
            header = fits.getheader( dirname+"/"+fil)
            mjd    = header['MJD-OBS']
            date   = header['DATE-OBS']
            sensor = header['HIERARCH ESO INS SENS151 VAL']
            demod = header['HIERARCH ESO PRO REC4 PARAM7 VALUE']
            print(demod)
            print(mjd)
            
            if mjd > mjd_range[0] and mjd < mjd_range[1] and demod == demodtype:
                
                MJD.append(mjd)
                print(date)
                fh = fits.open(dirname+"/"+fil)
                
                #wl = fh['OI_WAVELENGTH'].data
                WL = (fh[3].data)['EFF_WAVE'] # assuming the first extension is a table
                wlc = fh[3].columns
                    
                IDX = np.where(ar(WL>wlen_use[0]*1e-6) & ar(WL<wlen_use[1]*1e-6))
                print(len(IDX[0]))
               
                if len(IDX[0]) > 0:
                    
                    wl  = WL
                    idx = IDX;
                    
                    v2 = (fh[4].data)['VIS2DATA'] # assuming the first extension is a table
                    #v2 = np.reshape(v2,[6,v2.shape[0]/6,v2.shape[1]])
                    print(v2.shape)               
                    
#                    print(len(v2))
 #                   print(len(v2[1]))
                    V2 = np.append(V2,v2)
                    print(V2.shape)               
                    vd = (fh[6].data)['VISAMP'] # assuming the first extension is a table
                    for k in range(0,vd.shape[0]):
                        vd[k] = vd[k] / np.average(vd[k][idx])
                    VD = np.append(VD,vd)
                    
                    dp = (fh[6].data)['VISPHI'] # assuming the first extension is a table
                    
                    for k in range(0,dp.shape[0]):
                        # Correct residual OPD
                        coef = np.polyfit(1./wl[idx],dp[k,:][idx],1)
                        p    = np.poly1d(coef)
                        dp[k] = dp[k] - p(1./wl) 
                        
                    DP = np.append(DP,dp)
                    
                    
                    cp = (fh[5].data)['T3PHI'] # assuming the first extension is a table
                    CP = np.append(CP,cp)
                    fl = (fh[7].data)['FLUXDATA'] # assuming the first extension is a table
                    FL = np.append(FL,fl)
                    
                    v2c = fh[4].columns
                    #print(v2c)
                    
                    fh.close()

################################################################

V2_spec = 2*V_spec;
V2_goal = 2*V_goal;
V2_expected = 2*V_expected;
print("V Spec",V_spec,"V goal",V_goal, "V expected",V_expected)
print("V2 Spec",V2_spec,"V2 goal",V2_goal, "V2 expected",V2_expected)
print("CP Spec",CP_spec,"CP goal",CP_goal, "CP expected",CP_expected)
print("DP Spec",DP_spec,"DP goal",DP_goal, "DP expected",DP_expected)
print("VD Spec",VD_spec,"VD goal",VD_goal, "VD expected",VD_expected)

V2        = np.reshape(V2,[V2.shape[0]/v2.shape[1]/6,6,v2.shape[1]])
V         = np.sqrt(V2)
V2rms     = np.std(V2,0)
Vrms      = np.std(V,0)
V2avg     = np.average(V2,0)
Vavg      = np.average(V,0)
V2band    = np.median(V2avg[:,idx][:,0,:],1)
Vband     = np.median(Vavg[:,idx][:,0,:],1)
V2rmsband = np.average(V2rms[:,idx][:,0,:],1)
Vrmsband  = np.average(Vrms[:,idx][:,0,:],1)
print("sigma V2")
print(V2rmsband)
print("sigma V2 / V2")
print(V2rmsband / V2band)
V2ptp = np.ptp(V2,0)
V2ptpband = np.max(V2ptp[:,idx][:,0,:],1)
print("V2 PTP")
print(V2ptpband)

print("sigma V")
print(Vrmsband)

print("sigma V / V")
print(Vrmsband / Vband)

Vptp = np.ptp(V,0)
Vptpband = np.max(Vptp[:,idx][:,0,:],1)
print("V PTP")
print(Vptpband)

fig1 = plt.figure(0) # Here's the part I need
plt.clf()
ax1 = fig1.add_subplot(111)
plt.axhline(y=V2_spec, color='r', linestyle='-')
plt.axhline(y=V2_goal, color='b', linestyle='--')
plt.axhline(y=V2_expected, color='k', linestyle='-.')
    
for k,j in enumerate(V2rms[:,1]):
    plt.plot(wl*1e6,V2rms[k,:], color=cols[k])
plt.ylim((0,0.16))
plt.xlim(wlen_lim)
plt.ylabel('$\sigma_{V^2}$')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_V2_Stdev_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')



fig1a = plt.figure(99) # Here's the part I need
plt.clf()
ax1a = fig1a.add_subplot(111)
plt.axhline(y=V_spec, color='r', linestyle='-')
plt.axhline(y=V_goal, color='b', linestyle='--')
plt.axhline(y=V_expected, color='k', linestyle='-.')
    
for k,j in enumerate(Vrms[:,1]):
    plt.plot(wl*1e6,Vrms[k,:]/Vavg[k,:], color=cols[k])
plt.ylim((0,0.16))
plt.xlim(wlen_lim)
plt.ylabel('$\sigma_V / V$')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_V_RelErr_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')    
    

fig1ap = plt.figure(999) # Here's the part I need
plt.clf()
ax1ap = fig1ap.add_subplot(111)
plt.axhline(y=V_spec, color='r', linestyle='-')
plt.axhline(y=V_goal, color='b', linestyle='--')
plt.axhline(y=V_expected, color='k', linestyle='-.')
    
for k,j in enumerate(Vrms[:,1]):
    plt.plot(wl*1e6,Vrms[k,:], color=cols[k])
plt.ylim((0,0.16))
plt.xlim(wlen_lim)
plt.ylabel('$\sigma_V$')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_V_Stdev_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')    
    
fig1b = plt.figure(10) # Here's the part I need
plt.clf()
ax1b = fig1b.add_subplot(111)
for k,j in enumerate(V2[:,1]):
    for l,i in enumerate(V2[1,:,1]):
        plt.plot(wl*1e6,V[k,l,:], color=cols[l])
plt.ylim((-0.1,1.1))
plt.xlim(wlen_lim)
plt.ylabel('$V$')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_V2_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    
        
################################################################

CP = np.reshape(CP,[CP.shape[0]/cp.shape[1]/4,4,cp.shape[1]])  
CPrms = np.std(CP,0)
CPrmsband = np.median(CPrms[:,idx][:,0,:],1)
print("CP RMS")
print(CPrmsband)
CPptp = np.ptp(CP,0)
CPptpband = np.max(CPptp[:,idx][:,0,:],1)
print("CP PTP")
print(CPptpband)

plt.figure(3) # Here's the part I need
plt.clf()
plt.axhline(y=CP_spec, color='r', linestyle='-')
plt.axhline(y=CP_goal, color='b', linestyle='--')
plt.axhline(y=CP_expected, color='k', linestyle='-.')
for k,j in enumerate(CPrms[:,1]):
    plt.plot(wl*1e6,CPrms[k,:], color=cols[k])
plt.ylim((0,3))
plt.xlim(wlen_lim)
plt.ylabel('Closure phase st. dev. ($^o$)')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_CP_Stdev_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    
fig4b = plt.figure(13) # Here's the part I need
plt.clf()
ax4b = fig4b.add_subplot(111)
for k,j in enumerate(CP[:,1]):
    for l,i in enumerate(CP[1,:,1]):
        plt.plot(wl*1e6,CP[k,l,:], color=cols[l])
plt.ylim((-10,10))
plt.xlim(wlen_lim)
plt.ylabel('Closure phase ($^o$)')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_CP_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')

################################################################
        
DP = np.reshape(DP,[DP.shape[0]/dp.shape[1]/6,6,dp.shape[1]])  
DPrms = np.std(DP,0)
DPrmsband = np.median(DPrms[:,idx][:,0,:],1)
print("DP RMS")
print(DPrmsband)
DPptp = np.ptp(DP,0)
DPptpband = np.max(DPptp[:,idx][:,0,:],1)
print("DP PTP")
print(DPptpband)


fig2 = plt.figure(1) # Here's the part I need
plt.clf()
ax2 = fig2.add_subplot(111)
plt.axhline(y=DP_spec, color='r', linestyle='-')
plt.axhline(y=DP_goal, color='b', linestyle='--')
plt.axhline(y=DP_expected, color='k', linestyle='-.')
for p in [
patches.Rectangle(
        (wlen_use[1], 0), 10, 5,
        hatch='\\',
        fill=False
    ),
patches.Rectangle(
        (wlen_use[0], 0), -10, 5,
        hatch='\\',
        fill=False
    ),
]:
    ax2.add_patch(p    )
for k,j in enumerate(DPrms[:,1]):
    plt.plot(wl*1e6,DPrms[k,:], color=cols[k])
plt.ylim((0,3))
plt.xlim(wlen_lim)
plt.ylabel('Differential phase st. dev. ($^o$)')
plt.xlabel('Wavelength ($\mu$m)')
    
plt.savefig(os.getenv("HOME")+'/MATISSE_DP_Stdev_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    
fig2b = plt.figure(11) # Here's the part I need
plt.clf()
ax2b = fig2b.add_subplot(111)
for p in [
patches.Rectangle(
        (wlen_use[1], -180), 10, 360,
        hatch='\\',
        fill=False
    ),
patches.Rectangle(
        (wlen_use[0], -180), -10, 360,
        hatch='\\',
        fill=False
    ),
]:
    ax2b.add_patch(p    )
for k,j in enumerate(DP[:,1]):
    for l,i in enumerate(DP[1,:,1]):
        plt.plot(wl*1e6,DP[k,l,:], color=cols[l])
plt.ylim((-10,10))
plt.xlim(wlen_lim)
plt.ylabel('Differential phase ($^o$)')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_DP_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    
################################################################
        
VD = np.reshape(VD,[VD.shape[0]/vd.shape[1]/6,6,vd.shape[1]])  
VDrms = np.std(VD,0)
VDrmsband = np.median(VDrms[:,idx][:,0,:],1)
print("VD RMS")
print(VDrmsband)
VDptp = np.ptp(VD,0)
VDptpband = np.max(VDptp[:,idx][:,0,:],1)
print("VD PTP")
print(VDptpband)
fig3 = plt.figure(2) # Here's the part I need
plt.clf()
ax3 = fig3.add_subplot(111)
plt.axhline(y=VD_spec, color='r', linestyle='-')
plt.axhline(y=VD_goal, color='b', linestyle='--')
plt.axhline(y=VD_expected, color='k', linestyle='-.')
for p in [
patches.Rectangle(
        (wlen_use[1], 0), 10, 5,
        hatch='\\',
        fill=False
    ),
patches.Rectangle(
        (wlen_use[0], 0), -10, 5,
        hatch='\\',
        fill=False
    ),
]:
    ax3.add_patch(p    )
for k,j in enumerate(VDrms[:,1]):
    plt.plot(wl*1e6,VDrms[k,:], color=cols[k])
plt.ylim((0,0.16))
plt.xlim(wlen_lim)
plt.ylabel('Differential visibility st. dev.')
plt.xlabel('Wavelength ($\mu$m)')
        
plt.savefig(os.getenv("HOME")+'/MATISSE_VD_Stdev_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    
    
    
fig3b = plt.figure(12) # Here's the part I need
plt.clf()
ax3b = fig3b.add_subplot(111)
for p in [
patches.Rectangle(
        (wlen_use[1], -180), 10, 360,
        hatch='\\',
        fill=False
    ),
patches.Rectangle(
        (wlen_use[0], -180), -10, 360,
        hatch='\\',
        fill=False
    ),
]:
    ax3b.add_patch(p    )
for k,j in enumerate(VD[:,1]):
    for l,i in enumerate(VD[1,:,1]):
        plt.plot(wl*1e6,VD[k,l,:], color=cols[l])
plt.ylim((0,2))
plt.xlim(wlen_lim)
plt.ylabel('Differential visibility')
plt.xlabel('Wavelength ($\mu$m)')

plt.savefig(os.getenv("HOME")+'/MATISSE_VD_demod-'+demodtype+"_"+filtype+"_"+'-'.join([str(wle) for wle in wlen_use])+'microns.png')
    