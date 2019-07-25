# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:06:27 2016

@author: fmillour
"""
#final size of the interpolated + zero-padded images-
band  = 2
dimfx = 1024
dimfy = 1024
#P0    = 3.72e-2*1.15 # separation between 2 peaks in pixels
P0    = 4.29e-2 # separation between 2 peaks in pixels

# Multithread processing
import poppy
poppy.conf.use_multiprocessing = True
poppy.conf.use_fftw            = False

from matissecubedata import *

#dir="/home/fmillour/MATISSE/TESTS_NICE/OPD/data/BSL1/Calib_OPD_IP1_20_30000/"
#file=dir+'BCD3_IP15_20_29960b.fits'

dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/"
file  = dir+'N_4beams_dark.fits'
file2 = dir+'N_4beams.fits'
#file = dir+'N_2beams_dark.fits'
#file2 = dir+'N_2beams.fits'

#dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/CPN1/"
#file = dir+'CPN1_CPN2_IP57_dark.fits'
#file2 = dir+'Calib_OPD_IP3_100_29000BIS/CPN3_IP35_100_28800.fits'

dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/CPN2/"
file = dir+'CPN1_CPN2_IP57_dark.fits'
file2 = dir+'Calib_OPD_IP5_250_28000/CPN2_IP35_250_29750b.fits'
#file2 = dir+'Calib_OPD_IP5_250_28000/CPN2_IP35_250_30000b.fits'

dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/CPN2/"
file = dir+'CPN2_CPN3_IP35_dark.fits'
file2 = dir+'Calib_OPD_IP5_250_56000/CPN2_IP35_20_57500b.fits'

dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/CPN4/"
file = dir+'CPN4_CPN3_IP13_dark.fits'
file2 = dir+'Calib_OPD_IP1_100_32500/CPN4_IP13_100_32100.fits'


# Load the data
data  = MatisseCubeData(file, band)
data2 = MatisseCubeData(file2, band)

#Wavelength calibration 
#y_in   = [230, 897]
#lam_in = [8.68,12.7]
y_in   = [230, 424, 471, 578, 695, 897, 950]
lam_in = [8.68, 9.86, 10.13, 10.8, 11.53, 12.7, 13.09]
y0     = 53
lam    = MatisseLamFit(y_in,lam_in,y0,data.dimy)

#data.prepare(lam,dimfx,dimfy,P0)


# Plot the image
plt.figure();
plt.imshow(data.data[1,:,:])
plt.figure();
#plt.imshow(data2.data[1,:,:]-data.data[1,:,:],vmin=0,vmax=100);

# Subtract the dark frame
data.interf = data2.interf - np.median(data.interf, axis=0);
data.data   = data2.data   - np.median(data.data, axis=0);

plt.figure();
plt.imshow(data.interf[1,:,:],vmin=0,vmax=100);

# prepare the data
data.prepare(lam,dimfx,dimfy,P0)

# Compute the FFT
data.computeAll()

# Plot the FFT2D
plt.figure();
plt.imshow(np.log(np.abs(data.fft2D[1,:,:])))
plt.plot([data.P,data.P], [0, data.dimfx], 'r-', lw=1) # Red straight line
plt.plot([data.P[1],data.P[1]], [0, data.dimfx], 'g-', lw=1) # Red straight line
plt.plot([data.P[2],data.P[2]], [0, data.dimfx], 'b-', lw=1) # Red straight line
          


# Plot the resampled fringes
plt.figure();
plt.imshow(data.interfpFinal[1,:,:],vmin=0,vmax=100);

# Plot 
plt.figure();
plt.imshow(data.data[1,:500,:],vmin=0,vmax=100);
filepng = "/home/fmillour/MATISSE_N_Fringes.png"
plt.savefig(filepng)