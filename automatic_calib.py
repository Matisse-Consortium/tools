#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------------------
# Script to run all of the mat_cal_xxx  on the files produced by the automaticPipeline
#---------- Version 0, 10 Nov 2018 -------------------------------------------------------
#-----------------------------------------------------------------------------------------
from subprocess import call
import argparse
import glob
from astropy.io import fits

#--------------------------------------------------------------------------------
#------------ Set up the argument parser for command line arguments -------------
#---- This really shouldn't take many options besides the config file -----------
#---- and the input and output directories for data -----------------------------
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Wrapper to run the calibration steps of the MATISSE DRL.')
parser.add_argument('-o', '--out_dir', dest='out_dir', metavar='Output Directory', type=str, default='.', \
	help='The path to the directory you want results to be stored in (defaults to current directory).')
#-------------------------------------------------------------------------
#--- The SOF should contain TARGET_RAW_INT and CALIB_RAW_INT -------------
#-------------------------------------------------------------------------
#parser.add_argument('sof', type=str, help='The filepath containing your sof file. REQUIRED. Should contain TARGET_RAW_INT and CALIB_RAW_INT files.')
parser.add_argument('--target_dir', dest='targ_dir', metavar='Target Directory', type=str, default='.', \
	help='The path to the directory containing your TARGET intermediate products. Absolute filepath.')

parser.add_argument('--calib_dir', dest='calib_dir', metavar='Calibrator Directory', type=str, default='.', \
	help='The path to the directory containing your CALIB intermediate products. Absolute filepath.')
args = parser.parse_args()


def make_sof(input_dir, output_dir, targ_type):
	files =  glob.glob(input_dir+'/*.fits')
	print input_dir#, files
	fname ='%s/%s.sof'%(output_dir, targ_type)
	outfile = open(fname, 'w')
	for f in files:
		try:
			hdu = fits.open(f)
			obstype = hdu[0].header['pipefile'].split('0')[0]
			fmat = obstype[:-1]
			if 'RAW_INT' not in fmat:
				f = '#' + f
			outfile.write('{} \t {} \n'.format(f,fmat))
				#print '{} \t {}'.format(f,fmat), ' added'
		except:
			continue
	outfile.close()
	return fname



#---------------------------------------------------------------------------
#---- Make the SOF files ---------------------------------------------------
targsof = make_sof(args.targ_dir, args.out_dir, 'TARGET')
calibsof = make_sof(args.calib_dir, args.out_dir, 'CALIB')
combinedsof = '%s/%s.sof'%(args.out_dir, 'COMBINED')
call('cat %s %s > %s'%(targsof, calibsof, combinedsof), shell=True)




#--------------------------------------------------------------------------
#----- Run the Recipes ----------------------------------------------------
sof = combinedsof
outdir = args.out_dir
#mat_cal_cphase = "esorex --output-dir=%s mat_cal_cphase %s"%(outdir,sof)
#mat_cal_dphase = "esorex --output-dir=%s mat_cal_dphase %s"%(outdir,sof)
#mat_cal_vis = "esorex --output-dir=%s mat_cal_vis %s"%(outdir,sof)
mat_cal_oifits = "esorex --output-dir=%s mat_cal_oifits %s"%(outdir,sof)


#print 'Running mat_cal_cphase on sof:%s'%(sof)
#call(mat_cal_cphase, shell=True)
#print 'Running mat_cal_dphase on sof:%s'%(sof)
#call(mat_cal_dphase, shell=True)
#print 'Running mat_cal_vis on sof:%s'%(sof)
#call(mat_cal_vis, shell=True)
print 'Running mat_cal_oifits on sof:%s'%(sof)
call(mat_cal_oifits, shell=True)




