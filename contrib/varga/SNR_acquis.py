# 
"""
Signal to noise ratio calculationg tool for MATISSE.
@author: Jozsef Varga
varga.jozsef@csfk.mta.hu
"""

#TODO: output data files
#TODO: measure signal proerly in case of low flux when a Gaussian cannot be fitted
#TODO: calculate SNR of interferometric signal

import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits
import scipy.misc
from scipy.optimize import curve_fit

# basic parameters
basedir = 'D:/jvarga/Dokumentumok/MATISSE/data/MATISSE_ACQ/'
inputdir = 'D:/jvarga/Dokumentumok/MATISSE/data/MATISSE_ACQ/'
input_filename = 'MATISSE_ACQ_STD033_0002_SIM.fits'
output_plot1_name = "SNR_ACQ_per_frame"
output_plot2_name = "SNR_ACQ_avg"
output_plot3_name = "Noise_histogram"
output_plot4_name = "Signal_histogram"
colors = ['red','orange','green','blue','purple']

# open input fits file (acquisition images)
hdu_list = fits.open(inputdir + input_filename, mode='update')
hdr = hdu_list[0].header
n_frames = hdu_list['IMAGING_DATA'].header['NAXIS2']
n_subwindows = hdu_list['IMAGING_DETECTOR'].header['NREGION']
n_wavelengths = hdu_list['IMAGING_DATA'].data['DATA9'][0].shape[0]
n_signal_subwindows = 5
x_wl = np.arange(n_wavelengths)

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# Create empty arrays
std_noise_arr = np.zeros([n_subwindows])
mean_noise_arr = np.zeros([n_subwindows])
signal_arr = np.zeros([n_signal_subwindows,n_frames,n_wavelengths])
stacked_signal_arr = np.zeros([n_signal_subwindows,n_wavelengths])
std_noise_stacked_arr = np.zeros([n_signal_subwindows])
mean_noise_stacked_arr = np.zeros([n_signal_subwindows])

# open figure: Histogram of the stacked histogram
fig_sh = plt.figure()
ax = fig_sh.add_subplot(1, 1, 1)

# loop through all detector subwindows
for j in range(n_subwindows):
    print 'Subwindow ', j
    data_id_str = '%d' % (j + 1)
    # get imaging data
    img_data = hdu_list['IMAGING_DATA'].data['DATA' + data_id_str]

    # calculate the noise by fitting the histogram of all frames in a given subwindow
    #  histogram of the image
    hist,bin_edges = np.histogram(img_data,bins='auto')
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
        std_noise_arr[j] = coeff[2]
    except RuntimeError:
        print 'Curve fit failed'
        std_noise_arr[j] = np.nanstd(img_data)

    print 'Fitted mean = ', coeff[1]
    print 'Fitted standard deviation = ', coeff[2]

    # check if the the subwindow contains the photometric signal
    if ('PHOT' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j]):
        print 'Photometry subwindow'
        if (hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j] == 'PHOT1_1'):
            idx = 0
        elif (hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j] == 'PHOT2_1'):
            idx = 1
        elif (hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j] == 'PHOT3_1'):
            idx = 2
        elif (hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j] == 'PHOT4_1'):
            idx = 3
        else:
            idx = 0
        # measure the signal
        img_frame = hdu_list['IMAGING_DATA'].data['DATA' + data_id_str][0]
        # create empty array to store the stacked frames
        img_frame_summed = 1.0*np.zeros_like(img_frame)
        # loop through all frames
        for i in range(n_frames):
            #print "Frame ",i
            img_frame = hdu_list['IMAGING_DATA'].data['DATA' + data_id_str][i]
            img_frame_summed = img_frame_summed + img_frame
            #loop through each wavelength
            for k in range(n_wavelengths):
                #fit a 1D Gaussian to the image profile at the given wavelength to measure the peak value of the signal
                p0 = [np.amax(img_frame[k, :]),len(img_frame[k, :])/2.0,len(img_frame[k, :])/2.0]
                try:
                    coeff_signal, var_matrix_signal = curve_fit(gauss, np.arange(len(img_frame[k, :])), img_frame[k, :],p0=p0)
                    signal_arr[idx,i,k] = coeff_signal[0]
                except RuntimeError:
                    print 'Curve fit failed'
                    signal_arr[idx,i,k] = np.amax(img_frame[k, :])

        #now analyze the stacked signal
        # loop through each wavelength
        for k in range(n_wavelengths):
            # calculate the signal
            # fit a 1D Gaussian to the image profile at the given wavelength to measure the peak value of the signal
            p0 = [np.amax(img_frame_summed[k, :]), len(img_frame_summed[k, :]) / 2.0, len(img_frame_summed[k, :]) / 2.0]
            try:
                coeff_signal, var_matrix_signal = curve_fit(gauss, np.arange(len(img_frame_summed[k, :])), img_frame_summed[k, :],
                                                            p0=p0)
                stacked_signal_arr[idx, k] = coeff_signal[0]
            except RuntimeError:
                print 'Curve fit failed'
                stacked_signal_arr[idx, k] = np.amax(img_frame_summed[k, :])

        # calculate the noise of the stacked frame by fitting the histogram
        hist_s, bin_edges_s = np.histogram(img_frame_summed, bins='auto')
        bin_centres_s = (bin_edges_s[:-1] + bin_edges_s[1:]) / 2
        p0 = [np.nanmax(hist_s), bin_centres_s[hist_s == np.nanmax(hist_s)][0], np.nanstd(hist_s)]
        try:
            coeff_s, var_matrix_s = curve_fit(gauss, bin_centres_s, hist_s, p0=p0)
            gauss_fit_s = gauss(bin_centres_s, *coeff_s)
            std_noise_stacked_arr[idx] = coeff_s[2]
        except RuntimeError:
            print 'Curve fit failed'
            std_noise_stacked_arr[idx] = np.nanstd(img_frame_summed)

        #plot the histogram of the stacked frame along with the Gaussian fit
        plt.plot(bin_centres_s, 1.0 * hist_s / np.amax(hist_s),'-',color=colors[idx], label=('PHOT%d_1' % (idx+1)))
        plt.plot(bin_centres_s, gauss_fit_s / np.amax(gauss_fit_s),'--',color=colors[idx], label=('Gaussian fit PHOT%d_1' % (idx+1)))

    if ('INTERF' in hdu_list['IMAGING_DETECTOR'].data['REGNAME'][j]):
        print 'Interferometry subwindow'
        idx = 4
        # ...
        pass

hdu_list.close()

#save the stacked image histogram plot
ax.set_ylim([0.01,1.1])
leg = ax.legend(loc='best')  # loc='upper left')
leg.get_frame().set_alpha(0.5)
ax.set_yscale('log')
outputfig = basedir + output_plot4_name
fig_sh.savefig(outputfig + '.png', dpi=300)
fig_sh.savefig(outputfig + '.eps', format='eps', dpi=300)
plt.close(fig_sh)

#calculate signal-to-nois ratio (SNR)
#the final noise estimation will be just the average of the noise values over each subwindow
std_noise = np.median(std_noise_arr)
mean_noise = np.median(mean_noise_arr)
#print 'std_noise = ',std_noise
#SNR of the individual frames
SNR = (signal_arr-mean_noise)/std_noise
#SNR of the stacked frame
SNR_stacked = stacked_signal_arr
for j in range(n_signal_subwindows):
    SNR_stacked[j,:] = (stacked_signal_arr[j,:]-mean_noise_stacked_arr[j])/std_noise_stacked_arr[j]

#make plots
#SNR plot for each frame
fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(10, 8))
for j in range(n_signal_subwindows-1):
    ax1.plot(x_wl,SNR[j,0,:],'-',color=colors[j],lw=1,alpha=0.5,label=('PHOT%d_1' % (j+1)))
for i in range(1, n_frames):
    for j in range(n_signal_subwindows-1):
        ax1.plot(x_wl, SNR[j, i, :], '-', color=colors[j], lw=1, alpha=0.5, label='_')
#ax1.set_xlabel(r"$\lambda\mathrm{\ (}\mu\mathrm{m)}$")
ax1.set_xlabel(r'$Y$')
ax1.set_ylabel(r"${\mathrm{SNR}}$")
ax1.set_title(r"$\mathrm{SNR\ (individual\ frames)}$")
leg = ax1.legend(loc='best') #loc='upper left')
leg.get_frame().set_alpha(0.5)
ax1.set_ylim([0, 3.0*np.nanmedian(SNR[SNR>0])])
#ax1.set_xlim([np.nanmin(Q)-0.1,1.2])
#ax1.invert_yaxis()
plt.tight_layout()

outputfig = basedir + output_plot1_name
fig1.savefig(outputfig + '.png', dpi=300)
fig1.savefig(outputfig + '.eps', format='eps', dpi=300)
plt.close(fig1)

#SNR plot for the stacked frame
fig2, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(10, 8))
for j in range(n_signal_subwindows-1):
    ax1.plot(x_wl,SNR_stacked[j,:],'-',color=colors[j],lw=2,alpha=0.5,label=('PHOT%d_1' % (j+1)))
#ax1.set_xlabel(r"$\lambda\mathrm{\ (}\mu\mathrm{m)}$")
ax1.set_xlabel(r'$Y$')
ax1.set_ylabel(r"${\mathrm{SNR}}$")
ax1.set_title(r"$\mathrm{SNR\ (stacked)}$")
leg = ax1.legend(loc='best') #loc='upper left')
leg.get_frame().set_alpha(0.5)
ax1.set_ylim([0, 2.0*np.nanmedian(SNR_stacked)])
#ax1.set_xlim([np.nanmin(Q)-0.1,1.2])
#ax1.invert_yaxis()
plt.tight_layout()

outputfig = basedir + output_plot2_name
fig2.savefig(outputfig + '.png', dpi=300)
fig2.savefig(outputfig + '.eps', format='eps', dpi=300)
plt.close(fig2)

#plot a histogram of a selected subwindow to dfitting and noise estimationemonstrate the
fig3 = plt.figure()
ax = fig3.add_subplot(1, 1, 1)
plt.plot(bin_centres,1.0*hist/np.amax(hist),label='Histogram')
plt.plot(bin_centres_filt,1.0*hist_filt/np.amax(hist_filt),label='Clipped histogram')
# Now find the cdf
cdf = np.cumsum(hist)
# And finally plot the cdf
plt.plot(bin_centres, 1.0*cdf/np.amax(cdf),label='CDF')
#plt.plot(bin_centres[0:idx_max], 1.0*hist[0:idx_max]/np.amax(hist[0:idx_max]))
#plt.plot((bin_centres[idx_max],bin_centres[idx_max]),(0,1))
#plt.plot((bin_centres[int(idx_minus_1HWHM)],bin_centres[int(idx_minus_1HWHM)]),(0,1))
plt.plot(bin_centres,gauss_fit/np.amax(gauss_fit),label='Gaussian fit')
leg = ax.legend(loc='best') #loc='upper left')
leg.get_frame().set_alpha(0.5)
ax.set_yscale('log')

outputfig = basedir + output_plot3_name
fig3.savefig(outputfig + '.png', dpi=300)
fig3.savefig(outputfig + '.eps', format='eps', dpi=300)
plt.close(fig3)

print 'READY'