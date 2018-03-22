#!/usr/bin/python
"""
Structure function calculation to check time variability of the background for MATISSE
@author: Krisztina Gabanyi
gabanyi@konkoly.hu
"""

# Changelog:
# 2018-03-15: rewritten into function form (to be able to call it from system shell) (jvarga)
# 2018-03-15: images are not interpolated anymore (jvarga)
# 2018-03-15: separately plot SF for averaged image and for central pixel (jvarga)
# 2018-03-15: fix error when given timelag (timelag_img) is larger than the maximum timelag of the data (jvarga)
# 2018-03-21: new function arguments: window_area,n_max_frames, which can be used to speed up the analysis (jvarga)

# Caveats:
# - Currently only works if the observations are equally sequenced in time (within the given PREC)
# - When calculating the difference between the images the "+ 0.0" was needed for Python 2.7.13 :: Anaconda, Inc. to use values from the fits file as float and not integers.

# input parameters (all are optional)
# input_path: path to the input fits file
# output_path: path to the folder where the results are stored. if not given, then a folder is created with the name of the input file
# region_id: number of the detector region (starting from 1) (default = 1)
# timelag_img: Give a time-lag value (in sec) fow which the resulting 2 dimensional SF image is plotted (default = 0.5 s)
# bin_precision: Required precision for time-lag measurements, digit for seconds (default = 2, which means 0.01 s)
# window_area: edges of a selected window within the image [xmin,xmax,ymin,ymax] (default=None, which means the entire image)
# n_max_frames: number of frames to take into account (default=None, which means all frames)

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import os
import numpy as np
from astropy.io import fits
import fnmatch
from numpy import asarray as ar
import wx
import wx.xrc
import sys
import string

printable = set(string.printable)
# Basic parameters
# dirname = "/Users/gke/Downloads/MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001-10"
# outdirname = "/Volumes/KINGSTON/MATISSE/output"
# Which fits file to read in
# fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001.fits"
# fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0020.fits"

# open a file browse dialog to select which file to open
def get_path():
    app = wx.App(None)
    style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    dialog = wx.FileDialog(None, 'Open', style=style)  # wildcard=wildcard
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path

# window_area: [xmin,xmax,ymin,ymax]
def SF_calc(input_path="", output_path="", region_id=1, timelag_img=0.5, bin_precision=2, window_area=None, n_max_frames = None, \
            output_name_tag=""):
    if input_path == "":
        input_path = get_path()
    abs_path = os.path.abspath(input_path)
    input_path_dir = os.path.dirname(abs_path)
    input_filename = os.path.basename(abs_path)
    if output_path == "":
        # if output dir is not given, create a folder in the input directory to store the results
        try:
            output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SF'
            os.makedirs(output_dir)
        except OSError:
            # creating output_directory_failed
            if os.path.isdir(output_dir):
                output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SF'
            else:
                output_dir = input_path_dir
    else:
        output_dir = output_path
    # open input fits file
    hdul = fits.open(abs_path)  # dirname + "/" + fil)
    print 'Open file: ', abs_path
    imdata_header = hdul['IMAGING_DATA'].header
    # How many fields are there in the FITS table of the data
    NO_FIELDS = imdata_header.get('TFIELDS')
    # How many regions are there in the FITS table of the data
    NO_REGIONS = imdata_header.get('NREGION')
    # How many timeslots/rows in the FITS table of the data
    NO_ROWS = imdata_header.get('NAXIS2')
    regname = hdul['IMAGING_DETECTOR'].data['REGNAME']
    for i in range(len(regname)):
        regname[i] = filter(lambda x: x in printable, regname[i])
    if n_max_frames == None:
        n_max_frames = NO_ROWS
    print 'Number of frames:', NO_ROWS
    print 'Number of regions:', NO_REGIONS
    # Which detector region you want to get the SF
    # print "There are " + str(NO_REGIONS) + " detector regions in this FITS file."
    # r = int(raw_input("Which region are you interested in? "))
    if region_id > NO_REGIONS:
        print "Region ",region_id,". No such detector region. Using region ",NO_REGIONS," instead."
        region_id = NO_REGIONS
    print 'Analysing region ' + str(region_id) + ' - ' + regname[region_id - 1]

    if NO_ROWS > 5:
        # Create empty arrays to store the time-lags and the SF of the images
        time_lags = np.zeros((n_max_frames - 1))
        cdim = hdul['IMAGING_DATA'].data['DATA' + str(region_id)].shape
        xs = cdim[1]
        ys = cdim[2]
        print "Data cube dimensions (original): ", cdim
        if window_area == None:
            images = np.zeros((n_max_frames - 1, xs,ys))
        else:
            if window_area[3] > xs:
                window_area[3] = xs - 1
            if window_area[2] > xs:
                window_area[2] = xs - 1
            if window_area[2] >= window_area[3]:
                window_area[2] = window_area[3] - 1
            if window_area[1] > ys:
                window_area[1] = ys - 1
            if window_area[0] > ys:
                window_area[0] = ys - 1
            if window_area[0] >= window_area[1]:
                window_area[0] = window_area[1] - 1
            images = np.zeros((n_max_frames - 1, window_area[3]-window_area[2],window_area[1]-window_area[0]))


        # Calculate the time-lags and the SF of the images
        img_data = hdul['IMAGING_DATA'].data['DATA' + str(region_id)] * 1.0
        #print  img_data.shape
        if window_area == None:
            img_data = img_data[0:n_max_frames, :, :]
        else:
            img_data = img_data[0:n_max_frames, window_area[2]:window_area[3],window_area[0]:window_area[1]]
        print "Data cube dimensions (after windowing): ", images.shape
        #img_data =np.random.normal(0.0,1.0,img_data.shape) #for testing
        print  img_data.shape
        for i in range(n_max_frames - 1):
            print 'Frame: ', i
            # print 'hu'
            dt = np.unique(np.round((-hdul['IMAGING_DATA'].data['TIME'][0:(n_max_frames - 1 - i)] +
                                               hdul['IMAGING_DATA'].data['TIME'][(i + 1):n_max_frames]) * 3600.0 * 24.0, bin_precision))
            time_lags[i] =  dt[0]
            # print 'ha'
            images[i] = np.average(np.power((img_data[0:(n_max_frames - 1 - i)] - img_data[(i + 1):n_max_frames]), 2.0), axis=0)
        print "Timelag range: ", np.nanmin(time_lags), ' - ', np.nanmax(time_lags)
        if np.nanmax(time_lags) < timelag_img:
            timelag_img = np.nanmax(time_lags)
            # Save the time_lags and the SF of the image averages and SF of the central pixels to a text file
        f = open(output_dir + "/SF_data_region_" + regname[region_id - 1] + '_' + output_name_tag + '.dat', "w")
        f.write("#time lag, SF value averaged over the whole image, SF value for the central pixel\n")
        for i in range(n_max_frames - 1):
            f.write(str(time_lags[i]) + " " + str(np.average(images, axis=(1, 2))[i]) + " " + str(
                images[i, int(images.shape[1] * 0.5), int(images.shape[2] * 0.5)]) + "\n")
        f.close()

        # Create a plot for the SF of the averaged image and save to a png file
        fig = plt.figure()
        fig.suptitle('Structure function of detector region ' + regname[region_id - 1] + ' ' + output_name_tag)
        ax = fig.add_subplot(111)
        ax.set_xlabel('log$_{10}$(Time lag in sec)')
        ax.set_ylabel('log$_{10}$(SF)')
        ax.plot(np.log10(time_lags), np.log10(np.average(images, axis=(1, 2))), label='Averaged image')
        plt.legend(loc='best')
        plt.savefig(output_dir + '/SF_avg_region_' + regname[region_id - 1] + '_' + output_name_tag + '.png')

        # Create a plot for the SF of the central pixel and save to a png file
        fig = plt.figure()
        fig.suptitle('Structure function of detector region ' + regname[region_id - 1] + ' ' + output_name_tag)
        ax = fig.add_subplot(111)
        ax.set_xlabel('log$_{10}$(Time lag in sec)')
        ax.set_ylabel('log$_{10}$(SF)')
        ax.plot(np.log10(time_lags), np.log10(images[:, int(images.shape[1] * 0.5), int(images.shape[2] * 0.5)]),
                label='Central pixel')
        plt.legend(loc='best')
        plt.savefig(output_dir + '/SF_center_region_' + regname[region_id - 1] + '_' + output_name_tag + '.png')

        # Save the SF image for a particular time-lag (TL)
        fig2 = plt.figure()  # figsize=(10.0, ys/xs*10.0)
        fig2.suptitle('Structure function of detector region ' + regname[region_id - 1] + ' ' + output_name_tag + ', dt = ' + str(timelag_img) + 's (arcsinh scale)')
        # print images.shape
        # print time_lags.shape
        idx = (np.abs(time_lags - timelag_img)).argmin()
        # print idx
        plt.imshow(np.arcsinh(images[idx, :, :]), interpolation='None')  # clim = (0.0,1000)
        plt.colorbar()
        plt.savefig(
            output_dir + '/SF_timelag' + str(timelag_img) + '_region_' + regname[region_id - 1] + '_' + output_name_tag + '.png')

        hdul.close()

    else:
        print 'Too few frames (< 5) to calculate the structrure function.'
    print 'READY.'

##########################################################
# example usage:
#inpath=r"D:\jvarga\Dokumentumok\MATISSE\data\2018-03-16\MATISSE_GEN_MATISSE_SIPHOT_N_SKY_075_0001.fits"
#inpath=r"D:\jvarga\Dokumentumok\MATISSE\data\2018-03-16\MATISSE_OBS_HISENS_N_SKY_075_0001.fits"
# inpath=r"D:\jvarga\Dokumentumok\MATISSE\data\2018-03-16\MATISSE_OBS_SIPHOT_LM_SKY_075_0001.fits"
# output_name_tag = "corner" #"center","corner"
# window_center = [200, 280, 100, 170] #HISENS INTERF: [180, 290, 90, 160] #SIPHOT INTERF (reg 5):[200, 280, 100, 170] #xmin,xmax,ymin,ymax
# window_corner = [0, 100, 200, 250]
# window = window_corner
# SF_calc(input_path=inpath, output_path="", region_id=11, timelag_img=0.5, bin_precision=2, window_area=window,n_max_frames = None,\
#             output_name_tag=output_name_tag)

if __name__ == '__main__':
    list_arg = sys.argv  # print list_arg
    nargs = len(list_arg)
    if nargs == 1:
        SF_calc()
    elif nargs == 2:
        in_path = list_arg[1]
        SF_calc(in_path)
    elif nargs == 3:
        in_path = list_arg[1]
        out_path = list_arg[2]
        SF_calc(in_path, out_path)
    elif nargs == 4:
        in_path = list_arg[1]
        out_path = list_arg[2]
        r = list_arg[3]
        SF_calc(in_path, out_path, r)
    elif nargs == 5:
        in_path = list_arg[1]
        out_path = list_arg[2]
        r = list_arg[3]
        tl = list_arg[4]
        SF_calc(in_path, out_path, r, tl)
    elif nargs == 6:
        in_path = list_arg[1]
        out_path = list_arg[2]
        r = list_arg[3]
        tl = list_arg[4]
        prec = list_arg[5]
        SF_calc(in_path, out_path, r, tl, prec)
    elif nargs == 7:
        in_path = list_arg[1]
        out_path = list_arg[2]
        r = list_arg[3]
        tl = list_arg[4]
        prec = list_arg[5]
        win = list_arg[6]
        SF_calc(in_path, out_path, r, tl, prec,win)
    elif nargs == 8:
        in_path = list_arg[1]
        out_path = list_arg[2]
        r = list_arg[3]
        tl = list_arg[4]
        prec = list_arg[5]
        win = list_arg[6]
        n_max_frames = list_arg[7]
        SF_calc(in_path, out_path, r, tl, prec,win,n_max_frames)
    else:
        print "Wrong number of arguments."
