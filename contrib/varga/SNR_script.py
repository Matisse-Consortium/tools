# A script for the SNR analysis, using the SNR_acquis function.
# The result is the SNR of the acquisition images as function of the source flux.
# Jozsef Varga
# Konkoly Observatory
# e-mail: varga.jozsef@csfk.mta.hu

import numpy as np
import openpyxl as xl
import matplotlib.pyplot as plt
import SNR_acquis
import os
from astropy.io import fits

basedir = "D:/jvarga/Dokumentumok/MATISSE/"
datadir = basedir + "data/"
outputdir = basedir+"data/MATISSE_ACQ/"
input_xls_filename = "ObservationLogMatched1.xlsx"
sheetname = 'filelist_N'
outputfig1 = "MATISSE_ACQ_SNR_sec"
outputfig2 = "MATISSE_ACQ_SNR"
outputfile1_name = "MATISSE_ACQ_SNR_data.csv"
max_frame = -1
min_frame = 0

lw_global = 1.5
tfs = 10
tfs2 = 10
figx = 10 #(cm)
figy = 8 #(cm)
ms = 10 #markersize

min_row = 2
max_row = 19

#plot_tags = [r'$\mathrm{Total\ flux}$', r'$B_\mathrm{p} = 15\mathrm{\ m}$', r'$B_\mathrm{p} = 30\mathrm{\ m}$',
#             r'$B_\mathrm{p} = 60\mathrm{\ m}$', r'$B_\mathrm{p} = 100\mathrm{\ m}$', r'$B_\mathrm{p} = 130\mathrm{\ m}$',
#             r'$B_\mathrm{p} = 150\mathrm{\ m}$']

filename_col = 'D'
target_col = 'H'
L_col = 'I'
sigma_L_col = 'J'
N_col = 'K'
sigma_N_col = 'L'
date_folder_col = 'S'

#Second run: load excel sheet, and plot data
wb = xl.load_workbook(basedir+ 'comm/AIV_sky/' + input_xls_filename, data_only=True)
sheet = wb.get_sheet_by_name(sheetname)

filenames = np.chararray([max_row - min_row + 1], 32)
targets = np.chararray([max_row - min_row + 1], 32)
n_subwindows = 4
subwindow_labels=['PHOT1','PHOT2','PHOT3','PHOT4']
colors = ['red','orange','green','blue']

L = np.zeros([max_row - min_row + 1])
N = np.zeros([max_row - min_row + 1])
sigma_L = np.zeros([max_row - min_row + 1])
sigma_N = np.zeros([max_row - min_row + 1])
date_folders = np.chararray([max_row - min_row + 1], 32)
SNR_median_frames = np.zeros([max_row - min_row + 1,n_subwindows])
SNR_stacked_all = np.zeros([max_row - min_row + 1,n_subwindows])
SNRs_median_frames = np.zeros([max_row - min_row + 1,n_subwindows])
SNRs_stacked_all = np.zeros([max_row - min_row + 1,n_subwindows])
std_noise_stacked = np.zeros([max_row - min_row + 1,n_subwindows])
signal_stacked = np.zeros([max_row - min_row + 1,n_subwindows])
median_background_stacked = np.zeros([max_row - min_row + 1,n_subwindows])
DITs = np.zeros([max_row - min_row + 1])
n_files=len(filenames)
i=0
for row in range(min_row, max_row + 1):
    filenames[i] = sheet[filename_col + str(row)].value
    targets[i] = sheet[target_col + str(row)].value
    try:
        filenames[i] = sheet[filename_col + str(row)].value
        targets[i] = sheet[target_col + str(row)].value
        date_folders[i] = sheet[date_folder_col + str(row)].value
    except ValueError:
        filenames[i] = ""
        targets[i] = ""
        date_folders[i] = ""
    try:
        L[i] = sheet[L_col + str(row)].value
        sigma_L[i] = sheet[sigma_L_col + str(row)].value
    except ValueError:
        L[i] = np.NaN
        sigma_L[i] = np.NaN
    try:
        N[i] = sheet[N_col + str(row)].value
        sigma_N[i] = sheet[sigma_N_col + str(row)].value
    except ValueError:
        N[i] = np.NaN
        sigma_N[i] = np.NaN
    i = i + 1

for i in range(len(filenames)):
    filename = filenames[i]
    print i,filename
    input_path = datadir+date_folders[i]+'/'+filename
    if os.path.isfile(input_path):
        # open input fits file
        hdu_list = fits.open(input_path)
        DITs[i] =  hdu_list[0].header["HIERARCH ESO DET SEQ1 DIT"] #(s)
        print targets[i],DITs[i]
        hdu_list.close()

        [SNR_frames, SNR_stacked, std_noise_arr, median_background_arr, \
            std_noise_stacked_arr, median_background_stacked_arr,n_frames, \
            n_target_frames_stacked, n_sky_frames_stacked] = \
            SNR_acquis.SNR_acquis(input_path=input_path, output_path="", \
                   subwindows=[9,10,12,13], min_frame=min_frame, max_frame=max_frame, \
                   noise_calc_mode='hist', signal_calc_mode='gauss',stack_frames=True,make_plots=True)
        SNR_median_frames[i,:] = np.nanmedian(SNR_frames,axis=1)
        SNR_stacked_all[i,:] = SNR_stacked
        SNRs_median_frames[i, :] = np.nanmedian(SNR_frames, axis=1)*np.sqrt(1.0/DITs[i])
        SNRs_stacked_all[i, :] = SNR_stacked*np.sqrt(1.0/(n_target_frames_stacked*DITs[i]))
        std_noise_stacked[i, :] = std_noise_stacked_arr
        median_background_stacked[i, :] = median_background_stacked_arr

for j in range(n_subwindows):
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(figx, figy))
    #ax1.errorbar(L ,SNRs_median_frames[:,j] ,xerr=sigma_L, fmt='o', lw=lw_global, alpha=0.9, \
    #     label=subwindow_labels[j],color=colors[j],markersize=ms)
    ax1.errorbar(L, SNRs_stacked_all[:, j],xerr=sigma_L, fmt='d', lw=lw_global, alpha=0.9, \
                label=subwindow_labels[j]+' stacked', color=colors[j],markersize=ms)
    for k in range(n_files):
        if ((not np.isnan(L[k])) and (not np.isnan(SNRs_median_frames[k,j]))):
            ax1.text(L[k],SNRs_median_frames[k,j],targets[k].replace('_',' '),fontsize=tfs2)
        if ((not np.isnan(L[k])) and (not np.isnan(SNRs_stacked_all[k,j]))):
            ax1.text(L[k], SNRs_stacked_all[k, j], targets[k].replace('_',' '),fontsize=tfs2)
    #ax1.set_ylim((0, None))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$L\mathrm{\ (Jy)}$')
    ax1.set_ylabel(r"${\mathrm{SNR/s}}$")
    ax1.set_title(r"$\mathrm{SNR/s\ (vs.\ flux)}$")
    leg = ax1.legend(loc='best')  # loc='upper left')
    plt.tight_layout()

    outputfig = outputdir + '/' + outputfig1 + '_' + subwindow_labels[j]
    fig1.savefig(outputfig + '.png', dpi=300)
    fig1.savefig(outputfig + '.eps', format='eps', dpi=300)
    plt.close(fig1)


for j in range(n_subwindows):
    fig1, (ax1) = plt.subplots(1, 1, sharey=False, figsize=(figx, figy))
    #ax1.errorbar(L ,SNR_median_frames[:,j] ,xerr=sigma_L, fmt='o', lw=lw_global, alpha=0.9, \
    #     label=subwindow_labels[j],color=colors[j],markersize=ms)
    ax1.errorbar(L, SNR_stacked_all[:, j],xerr=sigma_L, fmt='d', lw=lw_global, alpha=0.9, \
                label=subwindow_labels[j]+' stacked', color=colors[j],markersize=ms)
    for k in range(n_files):
        if ((not np.isnan(L[k])) and (not np.isnan(SNR_median_frames[k,j]))):
            ax1.text(L[k],SNR_median_frames[k,j],targets[k].replace('_',' '),fontsize=tfs2)
        if ((not np.isnan(L[k])) and (not np.isnan(SNR_stacked_all[k,j]))):
            ax1.text(L[k], SNR_stacked_all[k, j], targets[k].replace('_',' '),fontsize=tfs2)
    #ax1.set_ylim((0,None))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$L\mathrm{\ (Jy)}$')
    ax1.set_ylabel(r"${\mathrm{SNR}}$")
    ax1.set_title(r"$\mathrm{SNR\ (vs.\ flux)}$")
    leg = ax1.legend(loc='best')  # loc='upper left')
    plt.tight_layout()

    outputfig = outputdir + '/' + outputfig2 + '_' + subwindow_labels[j]
    fig1.savefig(outputfig + '.png', dpi=300)
    fig1.savefig(outputfig + '.eps', format='eps', dpi=300)
    plt.close(fig1)

#open output data file
outputfile1 = outputdir + outputfile1_name
f_out1 = open(outputfile1, 'w+')
f_out1.write('#filename;date;target;L_flux(Jy);N_flux(Jy);sigma_L_flux(Jy);sigma_N_flux(Jy);'+\
        'subwindow;SNR_frame;SNR_stacked;SNRs_frame;SNRs_stacked;DIT;n_target_frames_stacked;std_noise_stacked;median_background_stacked\n')
for i in range(n_files):
    for j in range(n_subwindows):
        f_out1.write("%s;%s;%s;%f;%f;%f;%f;" % (filenames[i],date_folders[i],targets[i],\
            L[i],N[i],sigma_L[i],sigma_N[i]))
        f_out1.write("%s;%f;%f;%f;%f;%f;%d;%f;%f\n" % (subwindow_labels[j], \
            SNR_median_frames[i,j],SNR_stacked_all[i,j],\
            SNRs_median_frames[i,j],SNRs_stacked_all[i,j],DITs[i],n_target_frames_stacked, \
            std_noise_stacked[i, j],median_background_stacked[i,j] ) )
f_out1.close()

print 'READY.'