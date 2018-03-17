#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
  $Id: SNR_acquis.py 2018-03-14 author: jvarga

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2018- Observatoire de la Côte d'Azur

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

  $Author: jvarga $
  $Date: 2018-02-15 $
  $Revision: 3 $
  $Name:  $
  first line in case of Windows: #!/cygdrive/d/jvarga/Programok/Anaconda/python

  Signal to noise ratio calculationg tool for MATISSE.
  e-mail: varga.jozsef@csfk.mta.hu
"""

# SNR_calc
# Input parameters (all are optional)
# - input_path: path to the input fits file
# - output_path: path to the folder where the results are stored. if not given, then a folder is created with the name of the input file
# - subwindows: list containing the ID numbers of subwindows (starting from 1) in which we are interested (default: [9,10,12,13], only good for SiPhot)
# - min_frame: calculation runs from frame No. min_frame (default = first frame)
# - max_frame: calculation runs till frame No. max_frame (default = last frame)
# - noise_calc_mode: mode of the calculation of noise (default = 'std')
#   'std': noise is calculated by taking the std of the pixel values
#   'hist': noise is calculated by fitting Gaussian to the pixel histogram (slower)
# - signal_calc_mode: how to calculate the signal (default = 'int')
#   'max': signal is the max. value of the pixels
#   'gauss': signal is the max. value of a fitted 2D Gaussian
#   'int': signal is measured by integrating over an aperture
# - make_plots: indicate whether to make plots from the results (True/False, default: True)
# The script generates some plots.
#
# Output parameters:
# - SNR_frames: SNR of the individual frames, size: [n_subwindows,n_frames]
# - SNR_stacked: SNR of the stacked image (all frames summed), size: [n_subwindows]
# - std_noise_arr: standard deviation of the individual frames, size: [n_subwindows,n_frames]
# - median_background_arr, median level of the individual frames, size: [n_subwindows,n_frames]
# - std_noise_stacked_cal, standard deviation of the stacked image, size: [n_subwindows]
# - median_background_stacked_arr, median level of the stacked image, size: [n_subwindows]

# Chanegelog:
# - 2018-03-15: detect outlier pixels, and replace them using convolution
# - 2018-03-16: new signal measurement mode: 'int': the signal is integrated over in an elliptical aperture

#TODO: test and debug
#TODO: test if SNR scaling factor for "int" mode is suitable:
#TODO:    currently the integrated signal is divided by (2.0*np.pi *sigma_x*sigma_y) to estimate the
#TODO:    height of the signal. sigma_x and sigma_y are constants, and if the real PSF sigmas are different,
#TODO:    it will introduce an error in the SNR calculation
#TODO: select only central detector area to fit the Gaussian (to speed up the fitting)

import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits
import scipy.misc
from scipy.optimize import curve_fit
import sys
import wx
from astropy.stats import sigma_clip
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

# basic parameters
output_plot1_name = "_SNR_ACQ_per_frame"
output_plot2_name = "_SNR_ACQ_avg"
output_plot3_name = "_Noise_histogram"
output_plot4_name = "_stacked_frame"
#sigma values of the acquisiton image
sigma_x = 12*0.42 #0.42*lambda/D is the sigma of the Gaussian function which approximates the Airy-profile
sigma_y = 3*0.42
apert = 1.58 #optimum aperture radius for a Gaussian profile: 1.58 sigma
apert_factor = 2.0
bg_mask_factor = 4.0

#source: https://stackoverflow.com/questions/44865023/circular-masking-an-image-in-python-using-numpy-arrays
def createEllipticMask(h, w, center=None, ax=None,ay=None):
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if ax is None: # use the smallest distance between the center and image walls
        ax = min(center[0], center[1], w-center[0], h-center[1])
        ay = ax
    if ay is None:
        ay = ax

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2/ax**2 + (Y-center[1])**2/ay**2)

    mask = dist_from_center <= 1.0
    return mask

# Define model function to be used to fit to the data:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# Define model function to be used to fit to the data:
def gauss2d((x,y), *p):
    A, mu_x, sigma_x, mu_y, sigma_y,offset = p
    g = A*np.exp(-(x-mu_x)**2/(2.*sigma_x**2)-(y-mu_y)**2/(2.*sigma_y**2))+offset
    return g.ravel()

#open a file browse dialog to select which file to open
def get_path():
    app = wx.App(None)
    style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    dialog = wx.FileDialog(None, 'Open', style=style) #wildcard=wildcard
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path

#source: https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def SNR_acquis(input_path="", output_path="", \
               subwindows=[9,10,12,13], min_frame=0, max_frame=-1, \
               noise_calc_mode='std', signal_calc_mode='int',make_plots=True):
    if input_path == "":
          input_path = get_path()
    abs_path = os.path.abspath(input_path)
    input_path_dir = os.path.dirname(abs_path)
    input_filename = os.path.basename(abs_path)
    if make_plots == True:
        if output_path == "":
            #if output dir is not given, create a folder in the input directory to store the results
            try:
                output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SNR'
                os.makedirs(output_dir)
            except OSError:
                #creating output_directory_failed
                if os.path.isdir(output_dir):
                    output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SNR'
                else:
                    output_dir = input_path_dir
        else:
            output_dir = output_path
    # open input fits file
    hdu_list = fits.open(abs_path, mode='update')
    hdr = hdu_list[0].header
    n_frames = hdu_list['IMAGING_DATA'].header['NAXIS2']
    if max_frame == -1:
        max_frame = n_frames
    n_subwindows = hdu_list['IMAGING_DETECTOR'].header['NREGION']
    subwindows = np.array(subwindows)-1
    n_signal_subwindows = 0
    catg = subwindows*0 #0: CAL, 1: PHOT, 2: INT
    for i in range(len(subwindows)):
        if subwindows[i] > n_subwindows-1:
            subwindows[i] = n_subwindows-1
        if ('CAL' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[i]]):
            catg[i] = 0
        if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[i]]):
            catg[i] = 1
            n_signal_subwindows = n_signal_subwindows + 1
        if ('INT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[i]]):
            catg[i] = 2
            n_signal_subwindows = n_signal_subwindows + 1
    n_subwindows = len(subwindows)
    n_wavelengths = hdu_list['IMAGING_DATA'].data['DATA9'][0].shape[0]
    x_wl = np.arange(n_wavelengths)

    # Create empty arrays
    std_noise_arr = np.zeros([n_subwindows, len(range(min_frame, max_frame))]) * np.nan
    median_background_arr = np.zeros([n_subwindows, len(range(min_frame, max_frame))]) * np.nan
    signal_arr = np.zeros([n_subwindows, len(range(min_frame, max_frame))]) * np.nan
    stacked_signal_arr = np.zeros([n_subwindows]) * np.nan
    std_noise_stacked_arr = np.zeros([n_subwindows])*np.nan
    median_background_stacked_arr = np.zeros([n_subwindows])*np.nan

    # It is a 9x9 array
    kernel = Gaussian2DKernel(1.0)

    # loop through detector subwindows
    for j in range(n_subwindows):
        print 'Subwindow ', subwindows[j]+1,' ',hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]
        data_id_str = '%d' % (subwindows[j] + 1)
        # get imaging data
        img_data = 1.0*hdu_list['IMAGING_DATA'].data['DATA' + data_id_str]

        #determine and exclude outliers with sigma clipping
        print 'determine outliers with sigma clipping'
        for i in range(min_frame, max_frame):
            img_frame = img_data[i]
            clipped_frame = sigma_clip(img_frame, sigma=6)
            #check if there are outliers
            if clipped_frame.mask.any():
                print "Frame "+str(i)+': Outliers detected.'
                #set outliers to NaN
                img_frame[clipped_frame.mask] = np.NaN
                #convolve image with astropy's convolution tool
                # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
                img_frame_convolved = convolve(img_frame, kernel)
                #replace NaNs with the convolved values
                img_frame[clipped_frame.mask] = img_frame_convolved[clipped_frame.mask]

        x = np.arange(img_frame.shape[1])
        y = np.arange(img_frame.shape[0])
        x, y = np.meshgrid(x, y)
        len_x = img_frame.shape[0]
        len_y = img_frame.shape[1]
        x = x - len_y / 2.0
        y = y - len_x / 2.0
        if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]):
            imagemask = createEllipticMask(len_x, len_y,  ax=apert_factor*apert*sigma_x,ay=apert_factor*apert*sigma_y)
            background_mask = createEllipticMask(len_x, len_y, ax=bg_mask_factor*apert*sigma_x, ay=bg_mask_factor*apert*sigma_y)
            background_mask = ~background_mask
            n_pixels_mask = np.nansum(imagemask)

            #print (len_x * len_y)
            #print len(img_frame[background_mask])
            #print n_pixels_mask
            #plt.imshow(imagemask*img_frame)
            #plt.show()
        else:
            background_mask = np.ones((len_x, len_y), dtype=bool)

        # create empty array to store the stacked frames
        img_frame_summed = 1.0 * np.zeros_like(img_frame)
        # loop through all frames
        for i in range(min_frame, max_frame):
            # print "Frame ",i
            img_frame = img_data[i]
            # test image:
            #img_frame = np.random.normal(1000.0, 30.0, img_frame.shape)
            # add a Gaussian
            #gauss_img = 80000.0 * np.exp(-((x - 0.0) ** 2 / (2 * sigma_x ** 2) + (y - 0.0) ** 2 / (2 * sigma_y ** 2)))
            #img_frame = img_frame + gauss_img
            #plt.imshow(img_frame*imagemask)
            #plt.show()
            # plt.imshow(img_frame*background_mask)
            # plt.show()
            img_frame_summed = img_frame_summed + img_frame

            # if signal_calc_mode == 'gauss':
            #     #extract the center part of the image to speed up the fit
            #     print x.shape
            #     img_frame_cut = np.reshape( sigma_x*10.0,np.abs(y[0,:])<sigma_y*10.0)
            #     print  img_frame_cut.shape
            #     plt.imshow(img_frame_cut)
            #     plt.show()

            #measure noise
            if noise_calc_mode == 'gauss':
                # calculate the noise by fitting the histogram of all frames in a given subwindow
                #  histogram of the image
                hist,bin_edges = np.histogram(img_frame,bins='auto')
                bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

                # location of the maximum of the histogram
                idx_max = np.argmax(hist)

                # filter histogram in order to select only the noise part: select the part up to +1 HWHM
                F_halfmax = hist[idx_max] / 2.0            # Flux at the half maximum
                idx_minus_1HWHM = np.where(hist > F_halfmax)[0][0]
                diff_idx_1sigma = int(np.round(np.abs(idx_max-idx_minus_1HWHM)))

                #filtered histogram (only to +1 HWHM)
                hist_filt = hist[0:idx_max+diff_idx_1sigma]
                bin_centres_filt = bin_centres[0:idx_max+diff_idx_1sigma]

                #fit filtered histogram with a Gaussian
                # p0 is the initial guess for the fitting coefficients (A, mu and sigma)
                p0 = [np.nanmax(hist_filt),bin_centres_filt[hist_filt == np.nanmax(hist_filt)][0],np.nanstd(hist_filt)]
                try:
                    coeff, var_matrix = curve_fit(gauss, bin_centres_filt, hist_filt, p0=p0)
                    gauss_fit = gauss(bin_centres, *coeff)
                    #the width (std) of the fitted Gaussian is the noise
                    std_noise_arr[j,i] = coeff[2]
                    median_background_arr[j,i] = coeff[1]
                except RuntimeError:
                    print 'Noise calc: Curve fit failed (to histogram)'
                    std_noise_arr[j,i] = np.nanstd(img_frame)
                    median_background_arr[j,i] = np.nanmedian(img_frame)

                print 'Fitted mean = ', coeff[1]
                print 'Fitted standard deviation = ', coeff[2]
            elif noise_calc_mode == 'std':
                std_noise_arr[j, i] = np.nanstd(img_frame[background_mask])
                median_background_arr[j, i] = np.nanmedian(img_frame[background_mask])

            #measure signal
            if signal_calc_mode == "gauss":
                if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]):
                    #fit a 2D Gaussian to the image to measure the peak value of the signal
                    p0 = [np.nanmax(img_frame), 0.0, sigma_x, 0.0, sigma_y, np.nanmedian(img_frame)]
                    try:
                        coeff_signal, var_matrix_signal = curve_fit(gauss2d, (x,y), img_frame.flatten(),p0=p0)
                        signal_arr[j,i] = coeff_signal[0]+coeff_signal[5]
                        #print coeff_signal
                    except RuntimeError:
                        print 'Fit signal with Gaussian: fit failed'
                        signal_arr[j,i] = np.nanmax(img_frame)
                else:
                    signal_arr[j, i] = np.nanmax(img_frame)
            elif signal_calc_mode == "int":
                if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]):
                    signal_arr[j, i] = np.nansum(imagemask*(img_frame-median_background_arr[j, i]))
                else:
                    signal_arr[j, i] = np.nansum(img_frame-median_background_arr[j, i])
            elif signal_calc_mode == "max":
                signal_arr[j,i] = np.nanmax(img_frame)
            else:
                signal_arr[j, i] = np.nanmax(img_frame)

        #now analyze the stacked signal
        # calculate noise
        if noise_calc_mode == 'gauss':
            # calculate the noise of the stacked frame by fitting the histogram
            hist_s, bin_edges_s = np.histogram(img_frame_summed, bins='auto')
            bin_centres_s = (bin_edges_s[:-1] + bin_edges_s[1:]) / 2
            p0 = [np.nanmax(hist_s), bin_centres_s[hist_s == np.nanmax(hist_s)][0], np.nanstd(hist_s)]
            try:
                coeff_s, var_matrix_s = curve_fit(gauss, bin_centres_s, hist_s, p0=p0)
                gauss_fit_s = gauss(bin_centres_s, *coeff_s)
                std_noise_stacked_arr[j] = coeff_s[2]
                median_background_stacked_arr[j] = coeff_s[1]
            except RuntimeError:
                print 'Curve fit failed'
                std_noise_stacked_arr[j] = np.nanstd(img_frame_summed)
                median_background_stacked_arr[j] = np.nanmedian(img_frame_summed)
        elif noise_calc_mode == 'std':
            std_noise_stacked_arr[j] = np.nanstd(img_frame_summed[background_mask])
            median_background_stacked_arr[j] = np.nanmedian(img_frame_summed[background_mask])

        #calculate signal
        if signal_calc_mode == "gauss":
            if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]):
                # calculate the signal
                # fit a 2D Gaussian to the image to measure the peak value of the signal
                p0 = [np.nanmax(img_frame_summed), 0.0, sigma_x, 0.0,
                      sigma_y, np.nanmedian(img_frame_summed)]
                try:
                    coeff_signal, var_matrix_signal = curve_fit(gauss2d, (x, y), img_frame_summed.flatten(), p0=p0)
                    stacked_signal_arr[j] = coeff_signal[0]+coeff_signal[5]
                except RuntimeError:
                    print 'Fit stacked signal with Gaussian: fit failed'
                    stacked_signal_arr[j] = np.nanmax(img_frame_summed)
            else:
                stacked_signal_arr[j] = np.nanmax(img_frame_summed)
        elif signal_calc_mode == "int":
            if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]):
                stacked_signal_arr[j] = np.nansum(imagemask * (img_frame_summed-median_background_stacked_arr[j]))
            else:
                stacked_signal_arr[j] = np.nansum(img_frame_summed-median_background_stacked_arr[j])
        elif signal_calc_mode == "max":
            stacked_signal_arr[j] = np.nanmax(img_frame_summed)
        else:
            stacked_signal_arr[j] = np.nanmax(img_frame_summed)

        if make_plots == True:
            #figure: stacked image
            fig_sh = plt.figure()
            ax = fig_sh.add_subplot(1, 1, 1)
            plt.imshow(img_frame_summed,interpolation=None)
            plt.gca().invert_yaxis()
            plt.tight_layout()
            # save the stacked image
            outputfig = output_dir + '/' + input_filename.split('.')[0] + output_plot4_name + '_' + \
                        hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]
            fig_sh.savefig(outputfig + '.png', dpi=300)
            fig_sh.savefig(outputfig + '.eps', format='eps', dpi=300)
            plt.close(fig_sh)

    #calculate signal-to-noise ratio (SNR)
    #the final noise estimation will be just the average of the noise values over each subwindow
    std_noise_cal_arr = std_noise_arr #[catg==0,:]
    std_noise_cal = np.nanmedian(std_noise_cal_arr,axis=0)
    median_noise_cal_arr = median_background_arr #[catg==0,:]
    median_background_cal = np.nanmedian(median_noise_cal_arr,axis=0)
    std_noise_stacked_cal = np.nanmedian(std_noise_stacked_arr) #[catg==0]
    median_background_stacked_cal = np.nanmedian(median_background_stacked_arr) #[catg==0]

    print n_pixels_mask
    print 'std_noise = ',std_noise_cal
    print 'median_background = ', median_background_cal
    print 'std_noise_stacked = ', std_noise_stacked_cal
    print 'median_background_stacked = ', median_background_stacked_cal
    SNR_frames = np.zeros([n_subwindows, len(range(min_frame, max_frame))]) * np.nan
    for j in range(n_subwindows):
        print 'Subwindow ', subwindows[j] + 1, ' ', hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]
        print 'signal_frames = ', signal_arr[j,:]
        if signal_calc_mode == 'int':
            #SNR_frames[j, :] = (signal_arr[j, :]/(n_pixels_mask*1.0)) / std_noise_cal
            SNR_frames[j, :] = (signal_arr[j, :] / (2.0*np.pi *sigma_x*sigma_y)) / std_noise_cal #problematic if the real width of the brightness distribution is different
        else:
            SNR_frames[j,:] = (signal_arr[j,:]-median_background_cal)/std_noise_cal
        print 'SNR_frames = ', SNR_frames[j,:]
    if signal_calc_mode == 'int':
        #SNR_stacked = (stacked_signal_arr/(n_pixels_mask*1.0)) / std_noise_stacked_cal
        SNR_stacked = (stacked_signal_arr / (2.0*np.pi *sigma_x*sigma_y)) / std_noise_stacked_cal #problematic if the real width of the brightness distribution is different
    else:
        SNR_stacked = (stacked_signal_arr - median_background_stacked_cal) / std_noise_stacked_cal
    print 'SNR_stacked = ', SNR_stacked


    #make plots
    if make_plots == True:
        #SNR plot for each subwindow
        for j in range(n_subwindows):
            fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(10, 8))
            ax1.plot(np.arange(min_frame, max_frame), SNR_frames[j, :], '-', lw=1.5, alpha=0.9, \
                     label=hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]])
            #ax1.set_xlabel(r"$\lambda\mathrm{\ (}\mu\mathrm{m)}$")
            ax1.set_xlabel(r'$Frame$')
            ax1.set_ylabel(r"${\mathrm{SNR}}$")
            ax1.set_title(r"$\mathrm{SNR\ (individual\ frames)}$")
            leg = ax1.legend(loc='best') #loc='upper left')
            try:
                leg.get_frame().set_alpha(0.5)
            except AttributeError:
                pass
            #ax1.set_ylim([0, 3.0*np.nanmedian(SNR[SNR>0])])
            #ax1.set_xlim([np.nanmin(Q)-0.1,1.2])
            #ax1.invert_yaxis()
            plt.tight_layout()

            outputfig = output_dir + '/' + input_filename.split('.')[0] + output_plot1_name + '_' + hdu_list['IMAGING_DETECTOR'].data['REGNAME'][subwindows[j]]
            fig1.savefig(outputfig + '.png', dpi=300)
            fig1.savefig(outputfig + '.eps', format='eps', dpi=300)
            plt.close(fig1)

        #plot a histogram of a selected subwindow to demonstrate the fitting and noise estimation
        if noise_calc_mode == 'hist':
            fig3 = plt.figure()
            ax = fig3.add_subplot(1, 1, 1)
            plt.plot(bin_centres,1.0*hist/np.nanmax(hist),label='Histogram')
            plt.plot(bin_centres_filt,1.0*hist_filt/np.nanmax(hist_filt),label='Clipped histogram')
            # Now find the cdf
            cdf = np.cumsum(hist)
            # And finally plot the cdf
            plt.plot(bin_centres, 1.0*cdf/np.nanmax(cdf),label='CDF')
            #plt.plot(bin_centres[0:idx_max], 1.0*hist[0:idx_max]/np.nanmax(hist[0:idx_max]))
            #plt.plot((bin_centres[idx_max],bin_centres[idx_max]),(0,1))
            #plt.plot((bin_centres[int(idx_minus_1HWHM)],bin_centres[int(idx_minus_1HWHM)]),(0,1))
            plt.plot(bin_centres,gauss_fit/np.nanmax(gauss_fit),label='Gaussian fit')
            leg = ax.legend(loc='best') #loc='upper left')
            leg.get_frame().set_alpha(0.5)
            ax.set_yscale('log')

            outputfig = output_dir + '/' + input_filename.split('.')[0] + output_plot3_name
            fig3.savefig(outputfig + '.png', dpi=300)
            fig3.savefig(outputfig + '.eps', format='eps', dpi=300)
            plt.close(fig3)

    hdu_list.close()
    print 'READY'
    return SNR_frames,SNR_stacked,std_noise_arr,median_background_arr,std_noise_stacked_cal,median_background_stacked_arr

if __name__ == '__main__':
    list_arg = sys.argv
    #print list_arg
    nargs = len(list_arg)
    if nargs == 1:
        SNR_acquis()
    elif nargs == 2:
        SNR_acquis(list_arg[1])
    elif nargs == 3:
        SNR_acquis(list_arg[1],list_arg[2])
    elif nargs == 4:
        SNR_acquis(list_arg[1],list_arg[2],list_arg[3])
    elif nargs == 5:
        SNR_acquis(list_arg[1], list_arg[2],list_arg[3], list_arg[4])
    elif nargs == 6:
        SNR_acquis(list_arg[1], list_arg[2],list_arg[3], list_arg[4],list_arg[5])
    elif nargs == 7:
        SNR_acquis(list_arg[1], list_arg[2], list_arg[3], list_arg[4], list_arg[5],list_arg[6])
    elif nargs == 8:
        SNR_acquis(list_arg[1], list_arg[2], list_arg[3], list_arg[4], list_arg[5], list_arg[6],list_arg[7])
    elif nargs == 9:
        SNR_acquis(list_arg[1], list_arg[2], list_arg[3], list_arg[4], list_arg[5], list_arg[6],list_arg[7],list_arg[8])
    else:
        print "Wrong number of arguments."
