# -*- coding: utf-8 -*-
"""
Routines of general use for MATISSE data reduction
Created on Thu Jan 12 2017
Last version from Mon Oct 23 2017
@author: ama
"""

from   matplotlib import gridspec
import matplotlib.pyplot as plt
from   astropy.io import fits
import numpy as np
from   scipy.optimize import curve_fit
import scipy.ndimage.interpolation as ip
import glob,os


def MatisseLamFit(y_in,lam_in,y0,dimy):
#--------------------------------------------------------------------------------------------------------------
#Routine to create a dispersion law (quadratic) from the measurement of the position of a set of spectral lines
#--------------------------------------------------------------------------------------------------------------
 
#y_in: position in pixel of the spectral lines
#lam_in: corresponding wavelengths of the spectral lines
#y0: number of the first pixel of the spectral window
#dimy: number of pixels of the frame in the spectral direction
 
    a   = np.polyfit(y_in,lam_in,2);
    y   = np.arange(dimy)+y0;
    lam = a[2]+a[1]*y+a[0]*y**2;
    print("a = {0}".format(a[2]))    
    print("b = {0}".format(a[1]))
    print("c = {0}".format(a[0]))          
    return lam

  
def load_data_raw(filename):
#---------------------------------------------------------------------------------------------------
#Routine to load raw MATISSE fringe data from fits files produced by the transfer function template
#---------------------------------------------------------------------------------------------------

#filename: MATISSE raw data file (template format) before detector cosmetics

    hdulist = fits.open(filename)
    img = hdulist['IMAGING_DATA'].data['DATA11']
    nframe=np.size(img,axis=0)
    dimy=np.size(img,axis=1)
    dimx=np.size(img,axis=2)
    data=np.zeros([nframe,dimy,dimx],dtype=np.float)
    for i in range(10):
        data[i,:,:] = hdulist['IMAGING_DATA'].data['DATA11'][i]
    hdulist.close() 
    #Debug plot of the data frame
    #plt.figure()
    #plt.imshow(data[0,:,:],vmin=0)
    return data 
    

def load_data_cal_Siphot(filename):
#---------------------------------------------------------------------------------------------------------------------------
#Routine to load MATISSE fringe data (Si_phot mode) from fits files produced by the mat_cal_image (after detector cosmetics)
#---------------------------------------------------------------------------------------------------------------------------

#filename: MATISSE raw data file (template format) after detector cosmetics

    hdulist = fits.open(filename)
    img = hdulist['IMAGING_DATA'].data['DATA3']
    nframe=np.size(img,axis=0)
    nf=nframe/2
    dimy=np.size(img,axis=1)
    dimx=np.size(img,axis=2)
    data=np.zeros([nf,dimy,dimx],dtype=np.float)
    data_bg=np.zeros([nf,dimy,dimx],dtype=np.float)
    data_sub=np.zeros([nf,dimy,dimx],dtype=np.float)
    
    #Extraction of the raw frames
    for i in range(nf):
        data[i,:,:] = hdulist['IMAGING_DATA'].data['DATA3'][i]
    
    #Extraction of the background frames
    for i in range(nf):
        data_bg[i,:,:] = hdulist['IMAGING_DATA'].data['DATA3'][i+nf]
    hdulist.close()  
    
    #Subtraction of the background frames
    dark_mean=np.mean(data_bg,axis=0)
    for i in range(nf):
        data_sub[i,:,:]=data[i,:,:]-dark_mean
    
    #Debug plot of the data frame
    #plt.figure()
    #plt.imshow(data[0,:,:],vmin=0)
    return data_sub,nf,dimy,dimx
    
    

def load_data_cal_Highsens(filename):
#---------------------------------------------------------------------------------------------------------------------------
#Routine to load MATISSE fringe data (High_Sens mode) from fits files produced by the mat_cal_image (after detector cosmetics)
#---------------------------------------------------------------------------------------------------------------------------

#filename: MATISSE raw data file (template format) after detector cosmetics

    hdulist = fits.open(filename)
    img = hdulist['IMAGING_DATA'].data['DATA1']
    nframe=np.size(img,axis=0)
    nf=nframe/2
    dimy=np.size(img,axis=1)
    dimx=np.size(img,axis=2)
    data=np.zeros([nf,dimy,dimx],dtype=np.float)
    data_bg=np.zeros([nf,dimy,dimx],dtype=np.float)
    data_sub=np.zeros([nf,dimy,dimx],dtype=np.float)
    
    #Extraction of the raw frames
    for i in range(nf):
        data[i,:,:] = hdulist['IMAGING_DATA'].data['DATA1'][i]
    
    #Extraction of the background frames
    for i in range(nf):
        data_bg[i,:,:] = hdulist['IMAGING_DATA'].data['DATA1'][i+nf]
    hdulist.close()  
    
    #Subtraction of the background frames
    dark_mean=np.mean(data_dark,axis=0)
    for i in range(nf):
        data_sub[i,:,:]=data[i,:,:]-dark_mean
    
    #Debug plot of the data frame
    #plt.figure()
    #plt.imshow(data[0,:,:],vmin=0)
    return data_sub,nf,dimy,dimx 


def load_fft1D(filename):
#----------------------------------------------------------------------------------
#Routine to load the imaginary and real parts of the 1D FFT of MATISSE fringe data. 
#The 1D FFT is computed along the spatial direction for every spectral channel.
#----------------------------------------------------------------------------------

#filename: output file of the mat_est_corr recipe ("CORRFLUX_CAL_MATISSE_GEN_LAMP_*_*.fits")

    hdulist = fits.open(filename)
    real_part_fft1D = hdulist['OBJ_CORR_FLUX'].data['CORRFLUXREAL1']
    im_part_fft1D = hdulist['OBJ_CORR_FLUX'].data['CORRFLUXIMAG1']

    nframe=np.size(real_part_fft1D,axis=0) #Number of frames
    dimy=np.size(real_part_fft1D,axis=1) #Spectral direction
    dimf=np.size(real_part_fft1D,axis=2) #Spatial frequency direction

    return real_part_fft1D,im_part_fft1D 


def load_opd(filename,num_baseline):
#---------------------------------------------------------------------------------------------
#Routine to load the opd values (in um) computed for each fringe peak, without OPD modulation. 
#---------------------------------------------------------------------------------------------

#filename: output file of the mat_est_opd recipe ("OI_OPDWVPO_CORRFLUX_CAL_MATISSE_GEN_LAMP_*_*.fits")
#num_baseline: number of baselines

    hdulist = fits.open(filename)
    piston = hdulist['OI_OPD'].data['OPD']
    piston_err = hdulist['OI_OPD'].data['OPDERR']
    npiston = np.size(piston)
    nframe = npiston/num_baseline
    pistons = np.zeros([num_baseline,nframe],dtype=np.float)
    pistons_err = np.zeros([num_baseline,nframe],dtype=np.float)
    for i in range(num_baseline):
        pistons[i,:]=piston[i:npiston:num_baseline]
        pistons_err[i,:]=piston_err[i:npiston:num_baseline]
    
    return pistons,pistons_err 


def peak_fft1D(dimx,dimy,lam,lam_ref,peakNum):
#-------------------------------------------------------------------------------------------------------------------------
#Routine to compute the position (in pix) and width (in pix) of a given fringe peak at a given wavelength in the 1D FFT. 
#The 1D FFT is computed along the spatial direction for every spectral column.
#-------------------------------------------------------------------------------------------------------------------------

#dimx: dimension of the 1D FFT in the spatial (frequency) direction
#lam_tab: current wavelength [in um] 
#lam_ini: reference wavelength [in um] from which the position of the fringe peaks in the 1D FFT is computed 
#peakNum: peak number

    #Computation of position and width of the k^th fringe peak
    peakPosition=np.round((3*(peakNum)/(72.*(lam/lam_ini)))*dimx)+np.round(dimx/2)
    peakWidth=np.round((1./(72.*(lam/lam_ini)))*dimx)

    return peakPosition,peakWidth


       
def computeFft1D_spectral(interf,shift_im_min,shift_im_max,dimx,dimy,lam_tab,lam_ini):
#---------------------------------------------------------------------------------------------------------------------------------
#Routine to compute the 1D FFT (in the spatial direction) of a given MATISSE frame.
#Notably, this routine corrects the phase term introduced by a possible non-centering of the interferometric channel in the frame.
#---------------------------------------------------------------------------------------------------------------------------------  

#interf: MATISSE frame (2D array) corresponding to the interferometric channel
#shift_im_min: lower limit of the array of numerical shift values applied to the frame
#shift_im_max: upper limit of the array of numerical shift values applied to the frame
#dimx: dimension of the interferometric frame in the spatial direction
#dimy: dimension of the interferometric frame in the spectral direction
#lam_tab: wavelength array
#lam_ini: reference wavelength from which the position of the fringe peaks in the 1D FFT is computed 
    
    #Declaration of variables
    peakWidth=np.zeros(dimy)
    lower=np.zeros(dimy)
    upper=np.zeros(dimy)
    sigma_real=np.zeros(dimy)
    sigma_imag=np.zeros(dimy)
    npix=np.zeros(dimy,dtype=np.float)
    fft1D=np.zeros([dimy,dimx],dtype=np.complex)
    
    #Declaration of variables for chi2 calculation
    dimy_chi2=dimy#380
    shift_ima=np.zeros(dimy_chi2)
    
    #Declaration of variables for the numerical shift, denoted as "a", applied to recenter the interferometric channel     
    a=np.arange(shift_im_min,shift_im_max,0.1)
    dim_a=np.size(a)
    chi2_shift=np.zeros([dim_a,dimy_chi2],dtype=np.float)
    
    #Computation of the numerical shift (chi2 minimization on the imaginary part of the low-frequency peak in the 1D FFT) to be applied to the interferometric frame to recenter it
    for k in range(dimy_chi2):
        interf_dub=np.fft.fftshift(interf[k,:])
        fft1D_dub = np.fft.fftshift(np.fft.fft(interf_dub))
        peakWidth[k]=np.int((1./(72.*(lam_tab[k]/lam_ini)))*dimx)
        lower[k]=np.int(dimx/2)-peakWidth[k]
        upper[k]=np.int(dimx/2)+peakWidth[k]
        lower_short[k]=np.int(dimx/2)-peakWidth[k]/2
        upper_short[k]=np.int(dimx/2)+peakWidth[k]/2
        npix[k]=(np.int(upper[k])-np.int(lower[k]))
        sigma_real[k]=np.std(np.real(fft1D_dub[int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D_dub[int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma2_var=np.square(sigma_real[k])+np.square(sigma_imag[k])
        chi2_shift_min=1e10
        for l in range(dim_a):
            chi2_shift[l,k]=np.sum(np.square(np.imag(fft1D_dub[lower_short[k]:upper_short[k]+1]*np.exp(-1j*2.*np.pi*a[l]*np.arange(lower_short[k]-np.int(dimx/2),upper_short[k]+1-np.int(dimx/2))/dimx)))/sigma2_var)    
            #chi2_shift[l,k]=np.sum(np.square(np.imag(fft1D_dub[lower[k]:upper[k]]*np.exp(-1j*2.*np.pi*a[l]*np.arange(lower[k]-np.int(dimx/2),upper[k]-np.int(dimx/2))/dimx))))
            if chi2_shift[l,k] < chi2_shift_min:
                chi2_shift_min=chi2_shift[l,k]
                shift_ima[k]=a[l]
    shift_ima_mean=np.mean(shift_ima)

    
    #Centering of the interferometric frame by applying the correct numerical shift 'shift_ima_mean'
    for k in range(dimy):
        interf_final=np.fft.fftshift(ip.shift(interf[k,:],shift_ima_mean))
        fft1D[k,:] = np.fft.fftshift(np.fft.fft(interf_final))
    return fft1D
      
            
def computeFft1D(interf,dimx,dimy):
#----------------------------------------------------------------------------------
#Routine to compute the 1D FFT (in the spatial direction) of a given MATISSE frame.
#----------------------------------------------------------------------------------

#interf: MATISSE frame (2D array) corresponding to the interferometric channel
#dimx: dimension of the interferometric frame in the spatial direction
#dimy: dimension of the interferometric frame in the spectral direction
 
    fft1D=np.zeros([dimy,dimx],dtype=np.complex)
    for k in range(dimy):
        interf_dub=np.fft.fftshift(interf[k,:])
        fft1D[k,:] = np.fft.fftshift(np.fft.fft(interf_dub))
    return fft1D
  


def computeCorrFlux(fft1D,dimx,dimy,lam_tab,lam_ini,peakNum):
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#Routine to compute the unbiased energy (squared modulus) of the fringe peak of interest, and the energy of the low-frequency peak, at each wavelength. 
#It returns the correponding squared correlated flux (corr_flux), squared visibility (corr_flux_norm), and squared total flux (uncorr_flux).         
#----------------------------------------------------------------------------------------------------------------------------------------------------------

#fft1D: 1D FFT of the MATISSE interferometric frame for every spectral channel
#dimx: dimension of the 1D FFT in the spatial frequency direction
#dimy: dimension of the 1D FFT in the spectral direction
#lam_tab: wavelength array
#lam_ini: reference wavelength from which the position of the fringe peaks in the 1D FFT is computed 
#peakNum: number of the peak (=1, 2, 3, 4, 5, 6) for which the phasor is computed 

    #Declaration of variables
    peakPosition=np.zeros(dimy)
    peakWidth=np.zeros(dimy)
    lower=np.zeros(dimy)
    upper=np.zeros(dimy)
    lower_cent=np.zeros(dimy)
    upper_cent=np.zeros(dimy)
    upper=np.zeros(dimy)
    ind=np.zeros(dimy)
    npix=np.zeros(dimy,dtype=np.float) 
    sigma_corr_flux=np.zeros(dimy)
    bias_corr_flux=np.zeros(dimy)
    corr_flux=np.zeros(dimy,dtype=np.float)
    uncorr_flux=np.zeros(dimy,dtype=np.float)
    corr_flux_norm=np.zeros(dimy,dtype=np.float)
    shift_mean=100           
    for k in range(dimy):
        peakPosition[k]=np.round((3*(peakNum)/(72.*(lam_tab[k]/lam_ini)))*dimx)+np.round(dimx/2)
        peakWidth[k]=np.round((1./(72.*(lam_tab[k]/lam_ini)))*dimx)
        lower_cent[k]=np.round(dimx/2)-np.round(peakWidth[k])
        upper_cent[k]=np.round(dimx/2)+np.round(peakWidth[k]) 
        lower[k]=np.round(peakPosition[k]-peakWidth[k])
        upper[k]=np.round(peakPosition[k]+peakWidth[k])                
        ind[k]=np.argmax(abs(fft1D[np.roundint(lower[k]):np.round(upper[k])]))
        npix[k]=(np.round(upper[k])-np.round(lower[k]))
        int_mod=np.sum(np.square(np.abs(fft1D[k,np.round(lower[k]):np.round(upper[k])])))
        int_mod_cent=np.sum(np.square(np.abs(fft1D[k,np.round(lower_cent[k]):np.round(upper_cent[k])])))
        sigma_corr_flux[k]=np.std(np.square(np.abs(fft1D[k,np.round(lower[k]+2*peakWidth[k]+shift_mean):np.round(upper[k]+2*peakWidth[k]+shift_mean)])))*np.sqrt(npix[k])
        bias_corr_flux[k]=np.mean(np.square(np.abs(fft1D[k,np.round(lower[k]+2*peakWidth[k]+shift_mean):np.round(upper[k]+2*peakWidth[k]+shift_mean)])))
        corr_flux[k]=int_mod-bias_corr_flux[k]
        corr_flux_norm[k]=corr_flux[k]/int_mod_cent
        uncorr_flux[k]=int_mod_cent-bias_corr_flux[k]

    return (corr_flux,sigma_corr_flux,bias_corr_flux,corr_flux_norm,uncorr_flux)


        

def computePhasors_1D_spectral(fft1D,dimx,dimy,lam_tab,lam_ini,peakNum):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Routine to compute the complex phasor associated with the fringe peak of interest, at each wavelength. The fringe peak phasor is computed from the integration 
#of the real and imaginary parts of the 1D FFT over the fringe peak support. In order to remove the unknow phase origin 'affecting' the complex phasor,
# a spectral mathod is used. This consists in correlating the fringe peak phasor, at a given wavelength, with a reference phasor taken at a reference wavelength (band center).        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
#fft1D: 1D FFT of the MATISSE interferometric frame for every spectral channel
#dimx: dimension of the 1D FFT in the spatial frequency direction
#dimy: dimension of the 1D FFT in the spectral direction
#lam_tab: wavelength array
#lam_ini: reference wavelength from which the position of the fringe peaks in the 1D FFT is computed 
#peakNum: number of the peak (=1, 2, 3, 4, 5, 6) for which the phasor is computed 
    
    #Declaration of variables
    peakPosition=np.zeros(dimy)
    peakWidth=np.zeros(dimy)
    lower=np.zeros(dimy)
    upper=np.zeros(dimy)
    ind=np.zeros(dimy)
    npix=np.zeros(dimy,dtype=np.float) 
    sigma_real=np.zeros(dimy)
    sigma_imag=np.zeros(dimy)
    phasor1D=np.zeros(dimy,dtype=np.complex)
            
    for k in range(dimy):
        peakPosition[k]=np.int((3*(peakNum)/(72.*(lam_tab[k]/lam_ini)))*dimx)+np.int(dimx/2)
        peakWidth[k]=np.int((1./(72.*(lam_tab[k]/lam_ini)))*dimx)
        lower[k]=np.int(peakPosition[k]-peakWidth[k])
        upper[k]=np.int(peakPosition[k]+peakWidth[k])                
        ind[k]=np.argmax(abs(fft1D[int(lower[k]):int(upper[k])]))
        npix[k]=(np.int(upper[k])-np.int(lower[k]))

        #Computation of the phasor by integrating over the support of the fringe peak
        int_real=np.sum(np.real(fft1D[k,int(lower[k]):int(upper[k])]))
        int_imag=np.sum(np.imag(fft1D[k,int(lower[k]):int(upper[k])]))
        sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])        
        
        #Computation of the phasor by taking the maximum value of the fringe peak
        #int_real=np.real(fft1D[k,int(peakPosition[k])]) 
        #int_imag=np.imag(fft1D[k,int(peakPosition[k])]) 
        #sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))
        #sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)])) 
        phasor1D[k]=1j*int_imag+int_real
        
    phasor1D_ref=phasor1D[np.int(np.round(dimy/2))]  #Reference phasor
    phasor1D_diff=phasor1D*np.conjugate(phasor1D_ref)/np.abs(phasor1D_ref) #Intercorrelation                            

    return (phasor1D_diff,sigma_real,sigma_imag)
    
 
def computePhasors_1D(fft1D,dimx,dimy,lam_tab,lam_ini,peakNum):
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Routine to compute the complex phasor associated with the fringe peak of interest, at each wavelength. The fringe peak phasor is computed from the integration 
#of the real and imaginary parts of the 1D FFT over the fringe peak support). No correlation with a reference phasor at a reference wavelength is applied. However, 
#in order to remove the unknow phase origin 'affecting' the complex phasor, a differential mathod can be applied afterwards to the output of this routine.
#This consists in correlating the fringe peak phasor, at a given wavelength, with a reference phasor computed at the same wavelength on a reference frame 
#corresponding to a chosen voltage value, motor position.        
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
#fft1D: 1D FFT of the MATISSE interferometric frame for every spectral channel
#dimx: dimension of the 1D FFT in the spatial frequency direction
#dimy: dimension of the 1D FFT in the spectral direction
#lam_tab: wavelength array
#lam_ini: reference wavelength from which the position of the fringe peaks in the 1D FFT is computed 
#peakNum: number of the peak (=1, 2, 3, 4, 5, 6) for which the phasor is computed 

    #Declaration of variables
    peakPosition=np.zeros(dimy)
    peakWidth=np.zeros(dimy)
    lower=np.zeros(dimy)
    upper=np.zeros(dimy)
    ind=np.zeros(dimy)
    npix=np.zeros(dimy,dtype=np.float) 
    sigma_real=np.zeros(dimy)
    sigma_imag=np.zeros(dimy)
    phasor1D=np.zeros(dimy,dtype=np.complex)
            
    for k in range(dimy):
        peakPosition[k]=np.int((3*(peakNum)/(72.*(lam_tab[k]/lam_ini)))*dimx)+np.int(dimx/2)
        peakWidth[k]=np.int((1./(72.*(lam_tab[k]/lam_ini)))*dimx)
        lower[k]=np.int(peakPosition[k]-peakWidth[k])
        upper[k]=np.int(peakPosition[k]+peakWidth[k])                
        ind[k]=np.argmax(abs(fft1D[int(lower[k]):int(upper[k])]))
        npix[k]=(np.int(upper[k])-np.int(lower[k]))
        
        #Computation of the phasor by integrating over the support of the fringe peak
        int_real=np.sum(np.real(fft1D[k,int(lower[k]):int(upper[k])])) #Integration over the fringe peak support
        int_imag=np.sum(np.imag(fft1D[k,int(lower[k]):int(upper[k])])) #Integration over the fringe peak support
        sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        
        #Computation of the phasor by taking the maximum value of the fringe peak
        #int_real=np.real(fft1D[k,int(peakPosition[k])]) #value taken at the peak maximum
        #int_imag=np.imag(fft1D[k,int(peakPosition[k])]) #value taken at the peak maximum
        #sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))
        #sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))        
        
        phasor1D[k]=1j*int_imag+int_real

    return (phasor1D,sigma_real,sigma_imag)

        
    
def computePistons_1D_spectral(pistons,phasor1D,sigma_real,sigma_imag,lam_tab,dimy):
#----------------------------------------------------------------------------------------------------------------------
#Routine to compute the OPD (piston) associated with the differential complex phasor computed with the spectral method. 
#This piston is determined through a fit (over all the spectral channels) with a linear phasor model. The chi2 minimization 
#is performed by scanning a grid of piston values and finding the piston value that gives the lowest chi2 value.  
#----------------------------------------------------------------------------------------------------------------------        
    
#pistons: grid of piston values (in um) for the chi2 minimization (parameter space of the fit)
#phasor1D: differential phasor computed with the spectral method
#sigma_real: uncertainty on the real part of the differential complex phasor
#sigma_imag: uncertainty on the imaginary part of the differential complex phasor
#lam_tab: wavelength array
#dimy: dimension of the 1D FFT in the spectral direction
    
    # Wide scan
    npistons=np.size(pistons)
    chi2_phasor1D=np.zeros(npistons)
    for i in range(npistons):
        phasor1D_mod=np.exp(-1j*2.0*np.pi*pistons[i]*(1./lam_tab-1./lam_tab[np.int(np.round(dimy/2))]))
        phasor_tot=phasor1D*phasor1D_mod
        sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
        chi2_phasor1D[i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)            
    chi2_min=np.amin(chi2_phasor1D)
    ind=np.argmin(chi2_phasor1D)
    piston1D=pistons[ind]
    
    # Narrow scan around the global minimum
    pistons = np.arange(piston1D-5.,piston1D+5.,0.001,dtype=np.float)    
    npistons=np.size(pistons)
    chi2_phasor1D=np.zeros(npistons)
    for i in range(npistons):
        phasor1D_mod=np.exp(-1j*2.0*np.pi*pistons[i]*(1./lam_tab-1./lam_tab[np.int(np.round(dimy/2))]))
        phasor_tot=phasor1D*phasor1D_mod
        sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
        chi2_phasor1D[i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)
    
    #Chi2 minimization        
    chi2_min=np.amin(chi2_phasor1D)
    ind=np.argmin(chi2_phasor1D)
    piston1D=pistons[ind] #Piston value minimizing the chi2        
    phasor1D_mod_min=(np.exp(1j*2.0*np.pi*piston1D*(1./lam_tab-1./lam_tab[np.int(np.round(dimy/2))]))) #Corresponding differential phasor 
    #sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
    sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
    phasor_tot_min=(phasor1D*phasor1D_mod_min)
    sigma2_var_min=np.square(sigma_real*np.imag(phasor1D_mod_min))+np.square(sigma_imag*np.real(phasor1D_mod_min))  

    return(piston1D,sigma2_piston1D,phasor1D_mod_min,phasor_tot_min)
    

      
def computePistons_1D(pistons,phasor1D,sigma_real,sigma_imag,lam_tab,dimy):
#----------------------------------------------------------------------------------------------------------------------------
#Routine to compute the OPD (piston) associated with the absolute complex phasor computed with the computePhasors_1D routine. 
#This piston is determined through a fit (over all the spectral channels) with a linear phasor model. The chi2 minimization 
#is performed by scanning a grid of piston values and finding the piston value giving the lowest chi2 value.  
#----------------------------------------------------------------------------------------------------------------------------    
    
#pistons: grid of piston values (in um) for the chi2 minimization (parameter space for the fit)
#phasor1D: Absolute phasor provided by the computePhasors_1D routine
#sigma_real: uncertainty on the real part of the differential complex phasor
#sigma_imag: uncertainty on the imaginary part of the differential complex phasor
#lam_tab: wavelength array
#dimy: dimension of the 1D FFT in the spectral direction 
   
    # Wide scan
    npistons=np.size(pistons)
    chi2_phasor1D=np.zeros(npistons)
    for i in range(npistons):
        phasor1D_mod=np.exp(-1j*2.0*np.pi*pistons[i]*(1./lam_tab))
        phasor_tot=phasor1D*phasor1D_mod
        sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
        chi2_phasor1D[i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)            
    chi2_min=np.amin(chi2_phasor1D)
    ind=np.argmin(chi2_phasor1D)
    piston1D=pistons[ind]
    
    # Narrow scan around the global minimum
    pistons = np.arange(piston1D-5.,piston1D+5.,0.001,dtype=np.float)    
    npistons=np.size(pistons)
    chi2_phasor1D=np.zeros(npistons)
    for i in range(npistons):
        phasor1D_mod=np.exp(-1j*2.0*np.pi*pistons[i]*(1./lam_tab))
        phasor_tot=phasor1D*phasor1D_mod
        sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
        chi2_phasor1D[i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)

    #Chi2 minimization        
    chi2_min=np.amin(chi2_phasor1D)
    ind=np.argmin(chi2_phasor1D)
    piston1D=pistons[ind] #Piston value minimizing the chi2 
    phasor1D_mod_min=np.exp(1j*2.0*np.pi*piston1D*(1./lam_tab)) #Corresponding phasor
    #sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
    sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
    phasor_tot_min=(phasor1D/np.abs(phasor1D))*phasor1D_mod_min
    sigma2_var_min=np.square(sigma_real*np.imag(phasor1D_mod_min))+np.square(sigma_imag*np.real(phasor1D_mod_min))  
            
    return(piston1D,sigma2_piston1D,phasor1D_mod_min,phasor_tot_min)        
    
    

