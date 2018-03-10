#
"""
Structure function calculation to check time variability of the background for MATISSE
@author: Krisztina Gabanyi
gabanyi@konkoly.hu
"""

#Currently only work if times between the observations are the same (within the given PREC)
#PREC cannot be changed yet
#Too slow (several hours) for the largest detector region (DATA11) (acceptable for the others)
#Therefore, it does not cycle automatically through all the regions, but asks the user which one the work with
#Should output filename include part of the input filename?

#Note: When calculating the difference between the images the "+ 0.0" was needed for Python 2.7.13 :: Anaconda, Inc. to use values from the fits file as float and not integers.

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import os
import numpy as np
from astropy.io import fits
import fnmatch
from numpy import asarray as ar

#Basic parameters
dirname = "/Users/gke/Downloads/MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001-10"
outdirname = "/Volumes/KINGSTON/MATISSE/output"

#Which fits file to read in
#fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001.fits"
fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0020.fits"

#Give a time-lag value (in sec) fow which the resulting 2 dimensional SF image is plotted
TL = 4.0

#Required precision for time-lag measurements, digit for seconds - currently cannot be changed!
PREC = 2

# open input fits file
hdul=fits.open(dirname + "/" + fil)

#How many fields are there in the FITS table of the data
NO_FIELDS = hdul['IMAGING_DATA'].header.get('TFIELDS')
#How many regions are there in the FITS table of the data
NO_REGIONS = hdul['IMAGING_DATA'].header.get('NREGION')
#How many timeslots/rows in the FITS table of the data
NO_ROWS = hdul['IMAGING_DATA'].header.get('NAXIS2')

#Which detector region you want to get the SF
print "There are " + str(NO_REGIONS) + " detector regions in this FITS file."
r = int(raw_input("Which region are you interested in? "))

#Create empty arrays to store the time-lags and the SF of the images
time_lags = np.zeros((NO_ROWS - 1))
images = np.zeros((NO_ROWS - 1, hdul['IMAGING_DATA'].data['DATA' + str(r)].shape[1],hdul['IMAGING_DATA'].data['DATA' + str(r)].shape[2]))

#Calculate the time-lags and the SF of the images
for i in range(NO_ROWS - 1):
	time_lags[i] = np.unique(np.round((-hdul['IMAGING_DATA'].data['TIME'][0:(NO_ROWS - 1 - i)] + hdul['IMAGING_DATA'].data['TIME'][(i + 1):NO_ROWS])*3600.0*24.0,PREC))
	images[i] = np.average(np.power(((0.0 + hdul['IMAGING_DATA'].data['DATA' + str(r)][0:(NO_ROWS - 1 - i)]) - (0.0 + hdul['IMAGING_DATA'].data['DATA' + str(r)][(i + 1):NO_ROWS])),2.0),axis = 0)

#Save the time_lags and the SF of the image averages and SF of the central pixels to a text file
f = open(outdirname + "/SF_data_region" + str(r) + '.dat', "w")
f.write("#time lag, SF value averaged over the whole image, SF value for the central pixel\n")
for i in range(NO_ROWS - 1):
    f.write(str(time_lags[i]) + " " + str(np.average(images, axis=(1,2))[i]) + " " + str(images[i,int(images.shape[1] * 0.5),int(images.shape[2] * 0.5)]) + "\n") 
f.close()

#Create a plot for the SF of the averaged image and the central pixel and save to a png file
fig = plt.figure()
fig.suptitle('Structure function of detector region ' + str(r))
ax = fig.add_subplot(111)
ax.set_xlabel('log$_{10}$(Time lag in sec)')
ax.set_ylabel('log$_{10}$(SF)')
ax.plot(np.log10(time_lags), np.log10(np.average(images, axis=(1,2))), label = 'Averaged image')
ax.plot(np.log10(time_lags), np.log10(images[:,int(images.shape[1] * 0.5),int(images.shape[2] * 0.5)]), label = 'Central pixel')
plt.legend(loc = 'best')
plt.savefig(outdirname + '/SF_region' + str(r) + '.png')

#Save the SF image for a particular time-lag (TL)
fig2 = plt.figure()
fig2.suptitle('Structure function of detector region ' + str(r) + ' for timelag ' + str(TL) + 's')
plt.imshow(images[np.where(np.round(time_lags - TL, 1) == 0)[0][1]], clim = (0.0,1000))
plt.colorbar()
plt.savefig(outdirname + '/SF_timelag' + str(TL) + '_region' + str(r) + '.png')

hdul.close()



