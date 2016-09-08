# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:06:27 2016

@author: fmillour
"""
#final size of the interpolated + zero-padded images
band  = 1
dimfx = 1024
dimfy = 1024
P0    = 3.72e-2 # separation between 2 peaks in pixels

# Multithread processing
import poppy
poppy.conf.use_multiprocessing = True
poppy.conf.use_fftw            = False

from matissecubedata import *

#dir="/home/fmillour/MATISSE/TESTS_NICE/OPD/data/BSL1/Calib_OPD_IP1_20_30000/"
#file=dir+'BCD3_IP15_20_29960b.fits'

dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/"
file2 = dir+'L_fringesAnthony.fits'
file  = dir+'L_hotdarkAnthony.fits'

#dir   = "/home/fmillour/MATISSE/TESTS_NICE/OPD/data/BCOP1/Calib_OPD_IP1_100_30000/"
#file2  = dir+'BCD3_IP15_100_29800.fits'

# Load the data
#data  = MatisseCubeData(file,band)
data = MatisseCubeData(file2, band)

#Wavelength calibration 
y_in   = [1754,1498,1303]
lam_in = [2.87,3.23,3.5]
y0     = 1201
dimy   = 650
lam    = MatisseLamFit(y_in,lam_in,y0,data.dimy)

# Subtract the dark frame
#data.interf = data2.interf - np.median(data.interf, axis=0);

plt.figure();
plt.imshow(data.interf[1,:,:],vmin=0);

xline = 1/P0*dimfx/data.dimx+dimfx/2;

# prepare the data
data.prepare(lam,dimfx,dimfy,P0)

# Compute the FFT
data.computeAll()

# Plot the FFT2D
plt.figure();
plt.imshow(np.log(np.abs(data.fft2D[1,:,:])))
xline = 1/P0*dimfx/data.dimx;
d = +dimfx/2;

plt.plot([data.P,data.P], [0, dimfx], 'r-', lw=1) # Red straight line
plt.plot([data.P[1],data.P[1]], [0, dimfx], 'g-', lw=1) # Red straight line
plt.plot([data.P[2],data.P[2]], [0, dimfx], 'b-', lw=1) # Red straight line


# Plot 
plt.figure();
plt.imshow(data.data[1,:,:],vmin=0);
filepng = "/home/fmillour/MATISSE_L_Fringes.png"
plt.savefig(filepng)