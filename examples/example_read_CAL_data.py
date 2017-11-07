# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 17:01:09 2016

This is a small example file of how to use the fits module of astropy in
application to matisse calibrated (cosmetic applied) data


@author: ame
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Choose a MATISSE file and do the test 
dir="/home/fmillour/MATISSE/postPAE/AtelierCom/allDataTypes/"
filename=dir+"CAL_MATISSE_GEN_LAMP_L108_0001.fits"


hdulist      = fits.open(filename)           # Opening the file and putting it into a header data unit list
hdunames     = [hdu.name for hdu in hdulist] # Return the name of all the "extensions" in the fits file
primary      = hdulist['PRIMARY']            # The first table is normally the primary table which contains the main header but no table in the case of ESO data                  
header       = primary.header                # Return the header of the primary table as a kind of python dictionnary
chip         = header['HIERARCH ESO DET CHIP NAME'] # Return the value of the keyword HIERARCH ESO DET CHIP NAME
keys         = header.keys()                 # Return the list of all keywords in the primary header
#---------------------------------------------------------------------------------
imdet        = hdulist['IMAGING_DETECTOR']   # The second table is usually the imaging_detector table listing the detector regions recroded in the file. 
himdet       = imdet.header                  # Header of the imaging_detector table. Use is the same way as the primary table header
cols         = imdet.columns                 # Name and format of the columns in the table : 'REGION' 'DETECTOR' 'PORTS' 'CORRELATION' 'REGNAME' 'CORNER' 'GAIN' 'NAXIS' 'CRVAL' 'CRPIX' 'CTYPE' 'CD' 'DMP' 'DMC'
print cols[0]                                # Return the name and format of the first column of the data => 'REGION' and 'I' for int
print imdet.data                                   # the data in the table          
nrow         = len(imdet.data)               # Number of row in the data. In the case of the imaging_detector it correspond to the number of detector regions in the file.

print nrow
colnames     = imdet.data['REGNAME']

interfsize   = imdet.data[2].field('NAXIS') # return the [x,y] size of the 11th region of the detector (usually the interferometric one)
interfcorner = imdet.data[2].field('CORNER') # return the upper left corner position in absolute coordinate of all 11th detector region.
regnames     = imdet.data.field('REGNAME')   # return the names of the all the regions of the detectors 'CAL1' ... 'PHOT1_1' 'INTERF_1' ....
#---------------------------------------------------------------------------------                                              
imdata       = hdulist['IMAGING_DATA']       # The third table is usually the imaging_data table. You can check before in the hdunames list.
himdet       = imdet.header                  # Header of the imaging_data table. Use is the same way as the primary table header
cols         = imdet.columns                 # Name and format of the columns in the table : 'TIME' 'DATA1' ....  'DATA21' 'EXPTIME' 'OPT_TRAIN' 'OPD' 'LOCALOPD' 'INS_TRAIN' 'OFFSET' 'ROTATION' 'STEPPING_PHASE' 'TARTYP' 'REFERENCE' 'TARGET' 
len(imdata.data)                             # The number of row in the data. In the case of the imaging_data table it correspond to the number of frames recorded.
datainterf   = imdata.data.field('DATA3')   # Return the data cube of corresponding to the 11th region of the detector (usually the interferometric one)
shape        = np.shape(datainterf)          # Return the tuple of the size of the three dimensions of datainterf 


print shape;

plt.imshow(datainterf[5,:,:])