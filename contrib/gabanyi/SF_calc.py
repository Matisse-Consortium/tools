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

# Caveats:
# - Currently only works if the observations are equally sequenced in time (within the given PREC)
# - When calculating the difference between the images the "+ 0.0" was needed for Python 2.7.13 :: Anaconda, Inc. to use values from the fits file as float and not integers.

# input parameters (all are optional)
# input_path: path to the input fits file
# output_path: path to the folder where the results are stored. if not given, then a folder is created with the name of the input file
# r: number of the detector region (starting from 1) (default = 1)
# TL: Give a time-lag value (in sec) fow which the resulting 2 dimensional SF image is plotted (default = 4)
# PREC: Required precision for time-lag measurements, digit for seconds (default = 2, which means 0.01 s)

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import os
import numpy as np
from astropy.io import fits
import fnmatch
from numpy import asarray as ar
import wx
import sys

#Basic parameters
#dirname = "/Users/gke/Downloads/MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001-10"
#outdirname = "/Volumes/KINGSTON/MATISSE/output"
#Which fits file to read in
#fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0001.fits"
#fil = "MATISSE_GEN_MATISSE_SIPHOT_N_SKY_053_0020.fits"

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

def SF_calc(input_path="",output_path="",r=1,TL=4.0,PREC=2):
	if input_path == "":
		input_path = get_path()
	abs_path = os.path.abspath(input_path)
	input_path_dir = os.path.dirname(abs_path)
	input_filename = os.path.basename(abs_path)
	if output_path == "":
		# if output dir is not given, create a folder in the input directory to store the results
		try:
			output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SNR'
			os.makedirs(output_dir)
		except OSError:
			# creating output_directory_failed
			if os.path.isdir(output_dir):
				output_dir = input_path_dir + '/' + input_filename.split('.')[0] + '_SNR'
			else:
				output_dir = input_path_dir
	else:
		output_dir = output_path
	# open input fits file
	hdul=fits.open(abs_path) #dirname + "/" + fil)
	imdata_header = hdul['IMAGING_DATA'].header
	#How many fields are there in the FITS table of the data
	NO_FIELDS = imdata_header.get('TFIELDS')
	#How many regions are there in the FITS table of the data
	NO_REGIONS = imdata_header.get('NREGION')
	#How many timeslots/rows in the FITS table of the data
	NO_ROWS = imdata_header.get('NAXIS2')

	#Which detector region you want to get the SF
	#print "There are " + str(NO_REGIONS) + " detector regions in this FITS file."
	#r = int(raw_input("Which region are you interested in? "))
	print 'Analysing region '+ str(r) + ' - ' + hdul['IMAGING_DETECTOR'].data['REGNAME'][r-1]

	if NO_ROWS > 1:
		#Create empty arrays to store the time-lags and the SF of the images
		time_lags = np.zeros((NO_ROWS - 1))
		images = np.zeros((NO_ROWS - 1, hdul['IMAGING_DATA'].data['DATA' + str(r)].shape[1],hdul['IMAGING_DATA'].data['DATA' + str(r)].shape[2]))

		#Calculate the time-lags and the SF of the images
		for i in range(NO_ROWS - 1):
			time_lags[i] = np.unique(np.round((-hdul['IMAGING_DATA'].data['TIME'][0:(NO_ROWS - 1 - i)] + hdul['IMAGING_DATA'].data['TIME'][(i + 1):NO_ROWS])*3600.0*24.0,PREC))
			images[i] = np.average(np.power(((0.0 + hdul['IMAGING_DATA'].data['DATA' + str(r)][0:(NO_ROWS - 1 - i)]) - (0.0 + hdul['IMAGING_DATA'].data['DATA' + str(r)][(i + 1):NO_ROWS])),2.0),axis = 0)

		#Save the time_lags and the SF of the image averages and SF of the central pixels to a text file
		f = open(output_dir + "/SF_data_region" + str(r) + '.dat', "w")
		f.write("#time lag, SF value averaged over the whole image, SF value for the central pixel\n")
		for i in range(NO_ROWS - 1):
			f.write(str(time_lags[i]) + " " + str(np.average(images, axis=(1,2))[i]) + " " + str(images[i,int(images.shape[1] * 0.5),int(images.shape[2] * 0.5)]) + "\n")
		f.close()

		#Create a plot for the SF of the averaged image and save to a png file
		fig = plt.figure()
		fig.suptitle('Structure function of detector region ' + str(r))
		ax = fig.add_subplot(111)
		ax.set_xlabel('log$_{10}$(Time lag in sec)')
		ax.set_ylabel('log$_{10}$(SF)')
		ax.plot(np.log10(time_lags), np.log10(np.average(images, axis=(1,2))), label = 'Averaged image')
		plt.legend(loc = 'best')
		plt.savefig(output_dir + '/SF_avg_region' + str(r) + '.png')

		#Create a plot for the SF of the central pixel and save to a png file
		fig = plt.figure()
		fig.suptitle('Structure function of detector region ' + str(r))
		ax = fig.add_subplot(111)
		ax.set_xlabel('log$_{10}$(Time lag in sec)')
		ax.set_ylabel('log$_{10}$(SF)')
		ax.plot(np.log10(time_lags), np.log10(images[:,int(images.shape[1] * 0.5),int(images.shape[2] * 0.5)]), label = 'Central pixel')
		plt.legend(loc = 'best')
		plt.savefig(output_dir + '/SF_center_region_' + str(r) + '.png')

		#Save the SF image for a particular time-lag (TL)
		fig2 = plt.figure()
		fig2.suptitle('Structure function of detector region ' + str(r) + ' for timelag ' + str(TL) + 's')
		plt.imshow(images[np.where(np.round(time_lags - TL, 1) == 0)[0][1]], clim = (0.0,1000),interpolation='None')
		plt.colorbar()
		plt.savefig(output_dir + '/SF_timelag' + str(TL) + '_region' + str(r) + '.png')

		hdul.close()
	else:
		print 'Not enough data to calculate the structrure function.'

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
		SF_calc(in_path, out_path,r)
	elif nargs == 5:
		in_path = list_arg[1]
		out_path = list_arg[2]
		r = list_arg[3]
		tl = list_arg[4]
		SF_calc(in_path, out_path,r,tl)
	elif nargs == 6:
		in_path = list_arg[1]
		out_path = list_arg[2]
		r = list_arg[3]
		tl = list_arg[4]
		prec = list_arg[5]
		SF_calc(in_path, out_path,r,tl,prec)
	else:
		print "Wrong number of arguments."
