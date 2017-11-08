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
dir      = "/home/fmillour/MATISSE/postPAE/AtelierCom/allDataTypes/"
filename = dir+"PHOTBEAMS_CAL_MATISSE_GEN_FLUX_N108_0001.fits.gz"


hdulist      = fits.open(filename)           # Opening the file and putting it into a header data unit list
hdunames     = [hdu.name for hdu in hdulist] # Return the name of all the "extensions" in the fits file
primary      = hdulist['PRIMARY']            # The first table is normally the primary table which contains the main header but no table in the case of ESO data                  
header       = primary.header                # Return the header of the primary table as a kind of python dictionnary

protype      = header['HIERARCH ESO PRO TYPE'] # Return the value of the keyword HIERARCH ESO DET CHIP NAME
recipe       = header['HIERARCH ESO PRO REC1 ID'] # Return the value of the keyword HIERARCH ESO DET CHIP NAME
print "This data was", protype,"by",recipe

chip         = header['HIERARCH ESO DET CHIP NAME'] # Return the value of the keyword HIERARCH ESO DET CHIP NAME
print "Detector:", chip

keys         = header.keys()                 # Return the list of all keywords in the primary header
#---------------------------------------------------------------------------------
imdet        = hdulist['PHOT_BEAMS']   # The second table is usually the imaging_detector table listing the detector regions recroded in the file. 
himdet       = imdet.header                  # Header of the imaging_detector table. Use is the same way as the primary table header
cols         = imdet.columns                 # Name and format of the columns in the table : 'REGION' 'DETECTOR' 'PORTS' 'CORRELATION' 'REGNAME' 'CORNER' 'GAIN' 'NAXIS' 'CRVAL' 'CRPIX' 'CTYPE' 'CD' 'DMP' 'DMC'
print cols[0]                                # Return the name and format of the first column of the data => 'REGION' and 'I' for int
nrow         = len(imdet.data)               # Number of row in the data. In the case of the imaging_detector it correspond to the number of detector regions in the file.

print nrow
phot1     = imdet.data['DATANOSKY1']
phot2     = imdet.data['DATANOSKY2']
phot3     = imdet.data['DATANOSKY3']
phot4     = imdet.data['DATANOSKY4']


plt.imshow(phot1[0,:,:])