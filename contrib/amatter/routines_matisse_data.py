# -*- coding: utf-8 -*-
"""
Routines of general use for MATISSE data reduction
Created on Thu Jan 12 2017
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
    

def load_data_cal(filename):
#------------------------------------------------------------------------------------------------------------
#Routine to load MATISSE fringe data from fits files produced by the mat_cal_image (after detector cosmetics)
#------------------------------------------------------------------------------------------------------------

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
    hdulist.close()
    
    #Extraction of the background frames
    for i in range(nf):
        data_bg[i,:,:] = hdulist['IMAGING_DATA'].data['DATA3'][i+nf]
    hdulist.close()  
    
    #Subtraction of the background frames
    dark_mean=np.mean(data_dark,axis=0)
    for i in range(nf):
        data_sub[i,:,:]=data[i,:,:]-dark_mean
    
    #Debug plot of the data frame
    #plt.figure()
    #plt.imshow(data[0,:,:],vmin=0)
    return data_sub 

    
       
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
    print("shift_ima_mean = {0}").format(shift_ima_mean)
    #print("shift_ima = {0}").format(shift_ima)
    
    #Centering of the interferometric frame by applying the correct numerical shift 'shift_ima_mean'
    for k in range(dimy):
        interf_final=np.fft.fftshift(ip.shift(interf[k,:],shift_ima_mean))
        fft1D[k,:] = np.fft.fftshift(np.fft.fft(interf_final))
        #print("shift_ima[k] = {0}").format(shift_ima[k])
        #plt.figure()
        #plt.plot(a,chi2_shift[:,200])
        #tab_x=range(dimx)
        #plt.figure()
        #plt.imshow(np.sqrt(np.abs(fft1D)),vmin=0)
        #plt.plot(tab_x[int(lower[100]):int(upper[100])],np.abs(fft1D[100,int(lower[100]):int(upper[100])]))
        #plt.figure()
        #plt.plot(np.sqrt(np.abs(fft1D_dub)))
        #plt.plot(np.real(fft1D[10,:]))
        #plt.figure()
        #plt.plot(np.abs(fft1D[200,:]))
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
#----------------------------------------------------------------------------------------------------------------
#Routine to compute the energy (correlated flux) associated with the fringe peak of interest, at each wavelength.         
#----------------------------------------------------------------------------------------------------------------

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
    corr_flux=np.zeros(dimy,dtype=np.float)
            
    for k in range(dimy):
        peakPosition[k]=np.int((3*(peakNum)/(72.*(lam_tab[k]/lam_ini)))*dimx)+np.int(dimx/2)
        peakWidth[k]=np.int((1./(72.*(lam_tab[k]/lam_ini)))*dimx)
        lower[k]=np.int(peakPosition[k]-peakWidth[k])
        upper[k]=np.int(peakPosition[k]+peakWidth[k])                
        ind[k]=np.argmax(abs(fft1D[int(lower[k]):int(upper[k])]))
        npix[k]=(np.int(upper[k])-np.int(lower[k]))
        #int_real=np.sum(np.real(fft1D[k,int(lower[k]):int(upper[k])]))
        #int_imag=np.sum(np.imag(fft1D[k,int(lower[k]):int(upper[k])]))
        int_mod=np.sum(np.square(np.abs(fft1D[k,int(lower[k]):int(upper[k])])))
        sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        corr_flux[k]=int_mod

    #Various debug plots                            
    #tab_x=range(dimx)
    #plt.figure()
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(np.sqrt(np.square(sigma_imag)+np.square(sigma_real)))
    #print("sigma_imag[420] = {0}").format(sigma_imag[420])
    #print("sigma_imag[420] = {0}").format(sigma_real[420])
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(tab_x[int(lower[100]):int(upper[100])],np.abs(fft1D[100,int(lower[100]):int(upper[100])]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,int(lower[420]+2*peakWidth[420]+5):int(upper[420]+2*peakWidth[420]+5)]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,:]))
    return (corr_flux,sigma_real,sigma_imag)

        

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
        int_real=np.sum(np.real(fft1D[k,int(lower[k]):int(upper[k])]))
        int_imag=np.sum(np.imag(fft1D[k,int(lower[k]):int(upper[k])]))
        sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])        
        
        #int_real=np.real(fft1D[k,int(peakPosition[k])]) #Value taken at the maxium of the fringe peak
        #int_imag=np.imag(fft1D[k,int(peakPosition[k])]) #Value taken at the maxium of the fringe peak
        #sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))
        #sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)])) 
        phasor1D[k]=1j*int_imag+int_real
        
    phasor1D_ref=phasor1D[np.int(np.round(dimy/2))]  #Reference phasor
    phasor1D_diff=phasor1D*np.conjugate(phasor1D_ref)/np.abs(phasor1D_ref) #Intercorrelation                            

    #Various debug plots
    #tab_x=range(dimx)
    #plt.figure()
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(np.sqrt(np.square(sigma_imag)+np.square(sigma_real)))
    #print("sigma_imag[420] = {0}").format(sigma_imag[420])
    #print("sigma_imag[420] = {0}").format(sigma_real[420])
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(tab_x[int(lower[100]):int(upper[100])],np.abs(fft1D[100,int(lower[100]):int(upper[100])]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,int(lower[420]+2*peakWidth[420]+5):int(upper[420]+2*peakWidth[420]+5)]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,:]))
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
        
        int_real=np.sum(np.real(fft1D[k,int(lower[k]):int(upper[k])])) #Integration over the fringe peak support
        int_imag=np.sum(np.imag(fft1D[k,int(lower[k]):int(upper[k])])) #Integration over the fringe peak support
        sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))*np.sqrt(npix[k])
        
        #int_real=np.real(fft1D[k,int(peakPosition[k])]) #value taken at the peak maximum
        #int_imag=np.imag(fft1D[k,int(peakPosition[k])]) #value taken at the peak maximum
        #sigma_real[k]=np.std(np.real(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))
        #sigma_imag[k]=np.std(np.imag(fft1D[k,int(lower[k]+2*peakWidth[k]+5):int(upper[k]+2*peakWidth[k]+5)]))        
        
        phasor1D[k]=1j*int_imag+int_real

    #Various debug plots                            
    #tab_x=range(dimx)
    #plt.figure()
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(np.sqrt(np.square(sigma_imag)+np.square(sigma_real)))
    #print("sigma_imag[420] = {0}").format(sigma_imag[420])
    #print("sigma_imag[420] = {0}").format(sigma_real[420])
    #plt.plot(np.imag(phasor1D_diff))
    #plt.plot(tab_x[int(lower[100]):int(upper[100])],np.abs(fft1D[100,int(lower[100]):int(upper[100])]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,int(lower[420]+2*peakWidth[420]+5):int(upper[420]+2*peakWidth[420]+5)]))
    #plt.figure()
    #plt.plot(np.imag(fft1D[420,:]))
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
        #phasor1D_mod=np.exp(-1j*2.0*np.pi*pistons[i]*(1./lam_tab))
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
    
    #Various debug plots        
    #plt.figure()        
    #plt.plot(lam_tab,np.imag(phasor1D)/np.abs(phasor1D))    
    #plt.plot(lam_tab,np.imag(phasor1D_mod_min))    
    #plt.figure()        
    #plt.plot(lam_tab,np.real(phasor1D)/np.abs(phasor1D))    
    #plt.plot(lam_tab,np.real(phasor1D_mod_min))
    #plt.figure()    
    #plt.plot(1./lam_tab,np.angle(phasor1D_mod_min))
    #plt.plot(1./lam_tab,np.angle(phasor1D))
    #plt.figure()
    #plt.plot(pistons,chi2_phasor1D)    
    #plt.plot(lam_tab,np.square(np.imag(phasor_tot_min)))
    #plt.plot(lam_tab,sigma2_var_min)
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
            
    #plt.figure()        
    #plt.plot(lam_tab,np.imag(phasor1D)/np.abs(phasor1D))    
    #plt.plot(lam_tab,np.imag(phasor1D_mod_min))    
    #plt.figure()        
    #plt.plot(lam_tab,np.real(phasor1D)/np.abs(phasor1D))    
    #plt.plot(lam_tab,np.real(phasor1D_mod_min))           
    #plt.figure()    
    #plt.plot(1./lam_tab,np.angle(phasor1D_mod_min))
    #plt.plot(1./lam_tab,np.angle(phasor1D))
    #plt.figure()
    #plt.plot(pistons,chi2_phasor1D)    
    #plt.plot(lam_tab,np.square(np.imag(phasor_tot_min)))
    #plt.plot(lam_tab,sigma2_var_min)
    return(piston1D,sigma2_piston1D,phasor1D_mod_min,phasor_tot_min)        
    
    

#def computePistons_1D_spectral_wv(piston,phasor1D,sigma_real,sigma_imag,lam_tab,dimy):
##-----------------------------------------------------------------------------------------------------------------------------
##Routine to compute the OPD (piston) associated with the differential complex phasor computed with the spectral method. 
##This piston is determined through a fit (over all the spectral channels) with a quadratic phasor model. The chi2 minimization 
##is performed by scanning a grid of piston values and finding the piston value giving the lowest chi2 value.  
##-----------------------------------------------------------------------------------------------------------------------------    
#    
##pistons: grid of piston values (in um) for the chi2 minimization (parameter space for the fit)
##phasor1D: differential phasor computed with the spectral method
##sigma_real: uncertainty on the real part of the differential complex phasor
##sigma_imag: uncertainty on the imaginary part of the differential complex phasor
##lam_tab: wavelength array
##dimy: dimension of the 1D FFT in the spectral direction
#
#    pistons = np.arange(piston-15.0,piston+15.0,0.01,dtype=np.float)
#    npistons=np.size(pistons)
#    a = np.arange(-70.0,-49.0,0.5,dtype=np.float)
#    dima=np.size(a)
#    chi2_phasor1D=np.zeros([npistons,dima],dtype=float)   
#    chi2_min=1e12
#    for i in range(dima):
#        for j in range(npistons):
#            phasor1D_mod=np.exp(-1j*2.0*np.pi*(pistons[j]*(1./lam_tab-1./lam_tab[np.int(dimy/2)])+a[i]*(np.square(1./lam_tab)-np.square(1./lam_tab[np.int(dimy/2)]))))
#            phasor_tot=phasor1D*phasor1D_mod
#            sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
#            chi2_phasor1D[j,i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)
#            if (chi2_phasor1D[j,i] < chi2_min):
#                chi2_min=chi2_phasor1D[j,i]
#                j_min=j
#                i_min=i  
#    piston1D=pistons[j_min]
#    a1D=a[i_min]
#    phasor1D_mod_min=np.exp(1j*2.0*np.pi*(piston1D*(1./lam_tab-1./lam_tab[np.int(dimy/2)])+a1D*(np.square(1./lam_tab)-np.square(1./lam_tab[np.int(dimy/2)]))))
#    sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
#    phasor_tot_min=(phasor1D*phasor1D_mod_min)
#    sigma2_var_min=np.square(sigma_real*np.imag(phasor1D_mod_min))+np.square(sigma_imag*np.real(phasor1D_mod_min))  
#            
#    #plt.figure()        
#    #plt.plot(lam_tab,np.imag(phasor1D)/np.abs(phasor1D))    
#    #plt.plot(lam_tab,np.imag(phasor1D_mod_min))    
#    #plt.figure()        
#    #plt.plot(lam_tab,np.real(phasor1D)/np.abs(phasor1D))    
#    #plt.plot(lam_tab,np.real(phasor1D_mod_min))           
#    #plt.figure()    
#    #plt.plot(1./lam_tab,np.angle(phasor1D_mod_min))
#    #plt.plot(1./lam_tab,np.angle(phasor1D))
#    #plt.figure()
#    #plt.plot(a,chi2_phasor1D[j_min,:])    
#    #plt.plot(lam_tab,np.square(np.imag(phasor_tot_min)))
#    #plt.plot(lam_tab,sigma2_var_min)
#    return(piston1D,sigma2_piston1D,phasor1D_mod_min,phasor_tot_min)

#
##-----------------------------------------------------------------------------------------------------------------------------
##Routine to compute the OPD (piston) associated with the absolute complex phasor computed with the computePhasors_1D routine. 
##This piston is determined through a fit (over all the spectral channels) with a quadratic phasor model. The chi2 minimization 
##is performed by scanning a grid of piston values and finding the piston value giving the lowest chi2 value.  
##-----------------------------------------------------------------------------------------------------------------------------      
#def computePistons_1D_wv(piston,phasor1D,sigma_real,sigma_imag,lam_tab,dimy):
##pistons: grid of piston values (in um) for the chi2 minimization (parameter space for the fit)
##phasor1D: Absolute phasor provided by the computePhasors_1D routine
##sigma_real: uncertainty on the real part of the differential complex phasor
##sigma_imag: uncertainty on the imaginary part of the differential complex phasor
##lam_tab: wavelength array
##dimy: dimension of the 1D FFT in the spectral direction
# 
#    pistons = np.arange(piston-15.0,piston+15.0,0.01,dtype=np.float)
#    npistons=np.size(pistons)
#    a = np.arange(-70.0,-49.0,0.5,dtype=np.float)
#    dima=np.size(a)
#    chi2_phasor1D=np.zeros([npistons,dima],dtype=float)   
#    chi2_min=1e12
#    for i in range(dima):
#        for j in range(npistons):
#            phasor1D_mod=np.exp(-1j*2.0*np.pi*(pistons[j]*(1./lam_tab)+a[i]*(np.square(1./lam_tab))))
#            phasor_tot=phasor1D*phasor1D_mod
#            sigma2_var=np.square(sigma_real)+np.square(sigma_imag)                
#            chi2_phasor1D[j,i]=np.sum(np.square(np.imag(phasor_tot))/sigma2_var)
#            if (chi2_phasor1D[j,i] < chi2_min):
#                chi2_min=chi2_phasor1D[j,i]
#                j_min=j
#                i_min=i  
#    piston1D=pistons[j_min]
#    a1D=a[i_min]
#    phasor1D_mod_min=np.exp(1j*2.0*np.pi*(piston1D*(1./lam_tab)+a1D*(np.square(1./lam_tab))))
#    sigma2_piston1D=1./(np.sum((1./sigma2_var)*np.square(np.real(phasor1D)*np.real(phasor1D_mod_min)-np.imag(phasor1D)*np.imag(phasor1D_mod_min))))
#    phasor_tot_min=(phasor1D*phasor1D_mod_min)
#    sigma2_var_min=np.square(sigma_real*np.imag(phasor1D_mod_min))+np.square(sigma_imag*np.real(phasor1D_mod_min)) 
#            
#    #plt.figure()        
#    #plt.plot(lam_tab,np.imag(phasor1D)/np.abs(phasor1D))    
#    #plt.plot(lam_tab,np.imag(phasor1D_mod_min))    
#    #plt.figure()        
#    #plt.plot(lam_tab,np.real(phasor1D)/np.abs(phasor1D))    
#    #plt.plot(lam_tab,np.real(phasor1D_mod_min))
#    #plt.figure()    
#    #plt.plot(1./lam_tab,np.angle(phasor1D_mod_min))
#    #plt.plot(1./lam_tab,np.angle(phasor1D))
#    #plt.figure()
#    #plt.plot(pistons,chi2_phasor1D)    
#    #plt.plot(lam_tab,np.square(np.imag(phasor_tot_min)))
#    #plt.plot(lam_tab,sigma2_var_min)
#    return(piston1D,sigma2_piston1D,phasor1D_mod_min,phasor_tot_min) 
#
#
#    
#    
#def fit_plot_piston(dirdat,dirplot,fileForth,fileBack,loc_legend):
##    
#    xForth,pistonForth,pistonForthErr = np.loadtxt(fileForth)
#    xBack,pistonBack,pistonBackErr = np.loadtxt(fileBack)
#    
#    div    = fileForth.split("/")
##    #dstep  = div[-3]
#    #device = div[-2]
#    device = div[-2].split("_")[0]
#    filename_plot0=div[-1].split(".")[0].split("_")[0]
#    filename_plot1=div[-1].split(".")[0].split("_")[1]
#    #Fit de la loi voltage/piston
#    aForth,covForth  = np.polyfit(xForth,pistonForth,1,w=pistonForthErr,cov=True)
#    yForth           = aForth[0]*xForth+aForth[1]
#    
#    aBack,covBack  = np.polyfit(xBack,pistonBack,1,w=pistonBackErr,cov=True)
#    yBack           = aBack[0]*xBack+aBack[1]
#    
#    xGlobal=np.append(xForth,xBack)
#    pistonGlobal=np.append(pistonForth,pistonBack)
#    pistonGlobalErr=np.append(pistonForthErr,pistonBackErr)
#    aGlobal,covGlobal  = np.polyfit(xGlobal,pistonGlobal,1,w=pistonGlobalErr,cov=True)
#    yGlobal           = aGlobal[0]*xGlobal+aGlobal[1]    
#    
#    #Print of the results
#    print("Forth: a = {0:2.4f}$\pm${1:2.4f}".format(aForth[0],np.sqrt(covForth[0,0])))
#    print("Forth; b = {0:2.4f}$\pm${1:2.4f}".format(aForth[1],np.sqrt(covForth[1,1])))
#    print("Back: a = {0:2.4f}$\pm${1:2.4f}".format(aBack[0],np.sqrt(covBack[0,0])))
#    print("Back; b = {0:2.4f}$\pm${1:2.4f}".format(aBack[1],np.sqrt(covBack[1,1])))
#    print("Global: a = {0:2.4f}$\pm${1:2.4f}".format(aGlobal[0],np.sqrt(covGlobal[0,0])))
#    print("Global; b = {0:2.4f}$\pm${1:2.4f}".format(aGlobal[1],np.sqrt(covGlobal[1,1])))
#
##    
#    #Plot de la loi voltage/piston pour direction Forth et Back
#    fig  = plt.figure()
#    gs   = gridspec.GridSpec(3, 1)
#    plt1 = fig.add_subplot(gs[0:2,:])
#    plt.title("Calibration of Piezo-actuator {0}".format(device)) 
#    plt.ylabel(r"Measured piston ($\mu$m)")
#    dot_Forth,    = plt.plot(xForth,pistonForth, 'ro')
#    line_Forth,   = plt.plot(xForth,yForth,'r')
#    dot_back,  = plt.plot(xBack,pistonBack, 'go')
#    line_back, = plt.plot(xBack,yBack,'g')
#    plt.legend([(line_Forth,dot_Forth), (line_back,dot_back)], ['Forth $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aForth[0],np.sqrt(covForth[0,0])), 'Back $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aBack[0],np.sqrt(covBack[0,0]))],loc=loc_legend)
#    plt.legend([(line_Forth,dot_Forth), (line_back,dot_back)], ['Forth $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aForth[0],np.sqrt(covForth[0,0])), 'Back $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aBack[0],np.sqrt(covBack[0,0]))],loc=loc_legend)
##    
#    plt.tick_params(axis='x',which='both', bottom='off',top='off',labelbottom='off')
##    
#    #plt1.set_xticklabels([]);
#    fig.add_subplot(gs[2,:],sharex=plt1)
#    plt.plot(xForth,yForth-pistonForth,'ro')
#    plt.plot(xBack,yBack-pistonBack,'go')
#    plt.plot([xBack[0],xBack[-1]],[0,0],'k--')
#    plt.xlabel(r"Piezo Voltage (V)")               
#    plt.ylabel(r"Residual ($\mu$m)")    
#    
#    #Sauvegarde des plots
#    #fileeps = dirplot+"MATISSE_PiezoCalib_{0}.eps".format(device)
#    #plt.savefig(fileeps)
#    #filepng = dirplot+"MATISSE_PiezoCalib_{0}.png".format(device)
#    #plt.savefig(filepng)
#    fileeps = dirplot+filename_plot0+"_"+filename_plot1+"_{0}.eps".format(device)
#    plt.savefig(fileeps)
#    filepng = dirplot+filename_plot0+"_"+filename_plot1+"_{0}.png".format(device)
#    plt.savefig(filepng)
#    
#    #Plot de la loi voltage/piston global
#    fig  = plt.figure()
#    gs   = gridspec.GridSpec(3, 1)
#    plt1 = fig.add_subplot(gs[0:2,:])
#    plt.title("Calibration of Piezo-actuator {0}".format(device)) 
#    plt.ylabel(r"Measured piston ($\mu$m)")
#    dot_Global,    = plt.plot(xGlobal,pistonGlobal, 'ro')
#    line_Global,   = plt.plot(xGlobal,yGlobal,'r')
#    plt.legend([(line_Global,dot_Global)], ['Global $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aGlobal[0],np.sqrt(covGlobal[0,0]))],loc=loc_legend)
#    plt.legend([(line_Global,dot_Global)], ['Global $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/Volt'.format(aGlobal[0],np.sqrt(covGlobal[0,0]))],loc=loc_legend)
##    
#    plt.tick_params(axis='x',which='both', bottom='off',top='off',labelbottom='off')
##    
#    #plt1.set_xticklabels([]);
#    fig.add_subplot(gs[2,:],sharex=plt1)
#    plt.plot(xGlobal,yGlobal-pistonGlobal,'ro')
#    #plt.plot(xBack,yBack-pistonBack,'go')
#    plt.plot([xGlobal[0],xGlobal[-1]],[0,0],'k--')
#    plt.xlabel(r"Piezo Voltage (V)")
#    plt.ylabel(r"Residual ($\mu$m)")    
#    
#    #Sauvegarde des plots
#    #fileeps_Global = dirplot+"MATISSE_PiezoCalib_Global_{0}.eps".format(device)
#    #plt.savefig(fileeps_Global)
#    #filepng_Global = dirplot+"MATISSE_PiezoCalib_Global_{0}.png".format(device)
#    #plt.savefig(filepng_Global)
#    
#    fileeps_Global = dirplot+filename_plot0+"_"+filename_plot1+"_Global_{0}.eps".format(device)
#    plt.savefig(fileeps_Global)
#    filepng_Global = dirplot+filename_plot0+"_"+filename_plot1+"_Global_{0}.png".format(device)
#    plt.savefig(filepng_Global)
#       
#    
#    #Sauvegarde des résidus du fit
#    res_Forth=yForth-pistonForth
#    tabForth_residual   = np.array([xForth,res_Forth])
#    div_forth=fileForth.split("/")
#    fileforth_residual=dirdat+"residual_"+div_forth[-1]
#    np.savetxt(fileforth_residual,tabForth_residual)
#    res_Back=yBack-pistonBack
#    tabBack_residual   = np.array([xBack,res_Back])
#    div_back=fileBack.split("/")
#    fileback_residual=dirdat+"residual_"+div_back[-1]
#    np.savetxt(fileback_residual,tabBack_residual)
#    
#    
#    #Sauvegarde des plots des histogrammes des résidus du fit
#    plt.figure()
#    plt.hist(res_Forth,bins=10)    
#    fileeps = dirplot+filename_plot0+"_"+filename_plot1+"_res_Forth_{0}.eps".format(device)
#    plt.savefig(fileeps)
#    filepng = dirplot+filename_plot0+"_"+filename_plot1+"_res_Forth_{0}.png".format(device)
#    plt.savefig(filepng)
#    
#    plt.figure()
#    plt.hist(res_Back,bins=10)    
#    fileeps = dirplot+filename_plot0+"_"+filename_plot1+"_res_Back_{0}.eps".format(device)
#    plt.savefig(fileeps)
#    filepng = dirplot+filename_plot0+"_"+filename_plot1+"_res_Back _{0}.png".format(device)
#    plt.savefig(filepng)
#    
#    print("Forth: std deviation = {0}".format(np.std(yForth-pistonForth)))
#    print("Back: std deviation = {0}".format(np.std(yBack-pistonBack)))
#    print("Global: std deviation = {0}".format(np.std(yGlobal-pistonGlobal)))    
#    
#    
#def save_piston(filename,step,piston,piston_err):
#        #Ecriture des pistons dans un fichier
#        tab = np.array([step,piston,piston_err])
#        #print("stepBack = {0}".format(self.stepBack))
#        np.savetxt(filename,tab)
#        print("Saving piston data to {0}".format(filename))
#
#def save_piston_res(filename_res,step,piston_res):
#        #Ecriture des pistons dans un fichier
#        tab = np.array([step,piston_res])
#        #print("stepBack = {0}".format(self.stepBack))
#        np.savetxt(filename_res,tab)
#        print("Saving residual piston data to {0}".format(filename_res))
#
