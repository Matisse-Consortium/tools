# -*- coding: utf-8 -*-
"""
Routine converting MATISSE raw spectra (in ADU) in number of photons / s / m2 / um, and then in Jy
Created on Thu Mar 22 2018
@author: ama
"""

from   matplotlib import gridspec
import matplotlib.pyplot as plt
from   astropy.io import fits
import numpy as np
from   scipy.optimize import curve_fit
import scipy.ndimage.interpolation as ip
import glob,os
from scipy.constants import h,k,c


def convert_to_Jy(dir,nbeam,eta,phot_frac,tel_diam):
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Routine converting MATISSE raw spectra ('RAW_SPECTRUM_****.fits'), expressed in number of e-, in spectra expressed in Jy. 
#The Jy-converted raw spectra are saved in output fits files named 'RAW_SPECTRUM_****_Jy.fits'. 
#
#Info: An example of execution of this routine is given at the end of the file. Uncomment and edit the example if you want to run this routine directly from this file.
#
#Important : In order to be in unit of e-, the Matisse raw spectra must have been produced using the 'el' option of the mat_raw_estimates recipe. 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#INPUT:
#dir: directory where your MATISSE raw spectrum files ('RAW_SPECTRUM_****.fits') are put
#nbeam : number of beams (usually 4)
#eta : quantum efficiency of the detector (~0.9 for Hawaii-2RG, and ~0.5 for Aquarius)
#phot_frac : Fraction of flux sent to the photometric channels (~0.3 in Si_Phot mode, and =1 in High_Sens mode)  
#tel_diam : Telescope diameter (in m) 

#OUTPUT :
#Fits files with the following name:  'RAW_SPECTRUM_****_Jy.fits'
#

    raw_spectrum_Jy_files=glob.glob(dir+"RAW_SPECTRUM_*_Jy.fits")
    if (np.size(raw_spectrum_Jy_files) != -1):
        os.system("rm -f "+dir+"RAW_SPECTRUM_*_Jy.fits")
    raw_spectrum_files=glob.glob(dir+"RAW_SPECTRUM_*.fits")
    nfile_raw_spectrum = np.size(raw_spectrum_files)
    
    #Calculation of the telescope collecting area
    S_tel = np.pi*(tel_diam/2.)*(tel_diam/2.)
   
    
    #'Raw Spectrum' fits file reading
    for j in range(nfile_raw_spectrum):
        hdulist = fits.open(raw_spectrum_files[j])
        lam=hdulist['OI_WAVELENGTH'].data['EFF_WAVE'] # in m
        delta_lam=hdulist['OI_WAVELENGTH'].data['EFF_BAND'] # in m

        # wavelength array dimension
        dim_lam=np.size(lam)

        #Recalculation of the spectral channels width
        delta_lam[0] = abs(lam[1]-lam[0])
        for i in range(1,dim_lam):
            delta_lam[i] = abs(lam[i] - lam[i-1] )


        # Spectra extraction
        spec=np.zeros([nbeam,dim_lam],dtype=np.float)
        err_spec=np.zeros([nbeam,dim_lam],dtype=np.float)
        for i in range(nbeam):
            spec[i,:] = hdulist['OI_FLUX'].data['FLUXDATA'][i]
            err_spec[i,:] = hdulist['OI_FLUX'].data['FLUXERR'][i]

        # Detector integration time
        DIT=hdulist[0].header['HIERARCH ESO DET SEQ1 DIT']    #Detector Integration Time
        print('DIT = {0}').format(DIT)

        # Beam etendue
        #beam_etendue=np.square(1.5*np.pi/4)*np.square(3.5e-6)  #L band
        #Beam_etendue=np.square(2.0*np.pi/4.)*np.square(10.5e-6) # N band

        # Conversion to photons / s / m2 / m 
        conv_factor = DIT * delta_lam * eta * S_tel * phot_frac
        nphot= spec/conv_factor #number of photons / s / m2 / m
        err_nphot = err_spec/conv_factor

        # Conversion to W / m2 / m
        flux= nphot * h*c/lam  
        err_flux = err_nphot * h*c/lam  

        #Conversion to Jy
        flux_Jy = flux * (lam*lam/c) * 1e+26
        err_flux_Jy = err_flux * (lam*lam/c) * 1e+26
        
        #replacement of the raw spectra in the oiflux structure by the Jy-converted spectra 
        for i in range(nbeam):
            hdulist['OI_FLUX'].data['FLUXDATA'][i]=flux_Jy[i,:]  
            hdulist['OI_FLUX'].data['FLUXERR'][i]=err_flux_Jy[i,:]
            
        #Witing of the output files    
        filename = raw_spectrum_files[j].split("/")[-1]
        filename_str = filename.split(".")[0]
        filename_Jy=filename_str+'_Jy'+'.fits'
        hdulist.writeto(dir+filename_Jy)
        
        

    #Final Plots#    
    #plt.figure(0)
    #for i in range(nbeam-3):
    #    plt.plot(lam*1e+6,flux_Jy[i,:])    
    #    #plt.errorbar(lam,nphot[i,:],yerr=err_nphot)
    #    plt.xlabel('Wavelength [um]')
    #    plt.ylabel('Flux [Jy]')
    #plt.show()
    
    
#-------------------------------------
# Main (example of routine running)    
#-------------------------------------

#Data directory
#dir = '/home/amatter/Documents/Simu_Python/ThroughPut/2018-03-18/alf_ori_L/' 

#Number of beams
#nbeam=4

#quantum efficiency of the detector (=0.9 for Hawaii-2RG, and =0.5 for Aquarius)
#eta=0.9

#Fraction of flux sent to the photometric channels (=0.3 in Si_Phot, and =1 in High_Sens) 
#phot_frac=0.3

#Telescope diameter (in m)
#S_telescope = 1.8

#Call to the routine
#convert_to_Jy(dir,nbeam,eta,phot_frac,S_telescope)





    



