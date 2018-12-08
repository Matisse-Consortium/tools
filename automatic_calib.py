#!/usr/bin/env python
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# Script to run all of the mat_cal_xxx  on the files produced by the
#                            automaticPipeline
#                         Version 0, 10 Nov 2018
#------------------------------------------------------------------------------

from subprocess import call 
import argparse
import glob
import os
import numpy as np
from astropy.io import fits

#------------------------------------------------------------------------------
#------------ Set up the argument parser for command line arguments -----------
#---- This really shouldn't take many options besides the config file ---------
#---- and the input and output directories for data ---------------------------
#------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Wrapper to run the calibration steps of the MATISSE DRL.')

parser.add_argument('-o', '--out_dir', dest='out_dir', metavar='Output Directory', type=str, default='./CALIBRATED', \
    help='The path to the directory you want results to be stored in (defaults to current directory).')
#-------------------------------------------------------------------------
#--- The SOF should contain TARGET_RAW_INT and CALIB_RAW_INT -------------
#-------------------------------------------------------------------------
#parser.add_argument('sof', type=str, help='The filepath containing your sof file. REQUIRED. Should contain TARGET_RAW_INT and CALIB_RAW_INT files.')
parser.add_argument('--in_dir', dest='in_dir', metavar='Working Directory', type=str, default='.', \
    help='The path to the directory containing your oifits data.')


parser.add_argument('--timespan', dest='timespan', metavar='Time span of calibrators', type=str, default='.', \
    help='The time search interval for selecting calibrators around the science star')

args = parser.parse_args()


def make_sof(input_dir, output_dir, timespan=1./24.):

    SOFFILE = [];

    files =  glob.glob(input_dir+'/*.fits')
    print input_dir#, files

    DIC = []
    # First read all files
    print("Reading all files keywords...")
    for f in files:
        hdr     = fits.open(f)[0].header
        dic = {'hdr': hdr}
        DIC.append(dic)
        #try:
    print("Scanning files...")
    for i,f in enumerate(files):
            hdri = DIC[i]['hdr']
            obstype = hdri['ESO PRO CATG']

            if obstype == 'TARGET_RAW_INT':
                #print("\nFound a TARGET file! Working on it...")
                print("Working on", f)
                mjd  = hdri['MJD-OBS']
                bcd1 = hdri['ESO INS BCD1 NAME']
                bcd2 = hdri['ESO INS BCD2 NAME']
                chip = hdri['ESO DET CHIP NAME']
                dit  = hdri['ESO DET SEQ1 DIT']
                #print(mjd)

                filename  = os.path.basename(f)
                name, ext = os.path.splitext(filename)

                #print(filename)
                #print(name)
                #print(output_dir)
                fname ='%s/%s_cal_oifits.sof'%(output_dir, name)

                soffile = open(fname, 'w')
                SOFFILE.append(fname)
                soffile.write('{} \t {} \n'.format(f,obstype))

                calcount = 0;
                for j,fcal in enumerate(files):
                    if fcal != f:
                        hdrj = DIC[j]['hdr']
                        obstypec = hdrj['ESO PRO CATG']
                        mjdc  = hdrj['MJD-OBS']
                        bcd1c = hdrj['ESO INS BCD1 NAME']
                        bcd2c = hdrj['ESO INS BCD2 NAME']
                        chipc = hdrj['ESO DET CHIP NAME']
                        ditc  = hdrj['ESO DET SEQ1 DIT']
                        dif = mjd - mjdc
                        absdif = np.abs(dif)

                        if obstypec == 'CALIB_RAW_INT' and bcd1 == bcd1c and bcd2 == bcd2c and chip == chipc and dit == ditc:
                            if absdif < timespan:
                                #print(fcal)
                                #print(mjdc)
                                #print(dif)
                                soffile.write('{} \t {} \n'.format(fcal,obstypec))
                                calcount+=1

                soffile.close()
                print("Found",calcount,"suitable calibrators")

            fmat = obstype[:-1]
            if 'RAW_INT' not in fmat:
                f = '#' + f
                #print '{} \t {}'.format(f,fmat), ' added'
        #except:
         #   continue
    return SOFFILE

outdir = args.out_dir

if not os.path.exists(outdir):
    os.makedirs(outdir)

#---------------------------------------------------------------------------
#---- Make the SOF files ---------------------------------------------------
targsof = make_sof(args.in_dir, args.out_dir)


#--------------------------------------------------------------------------
#----- Run the Recipes ----------------------------------------------------
for isof in targsof:
    print 'Running mat_cal_oifits on sof:%s'%(isof)
    call("esorex --output-dir=%s mat_cal_oifits %s"%(outdir,isof), shell=True)

    name, ext = os.path.splitext(isof)
    #print(name)

    # Rename files
    resultFiles = glob.glob(outdir+'/TARGET_CAL_INT_????.fits')
    #print(resultFiles)
    for idx,fi in enumerate(resultFiles):
        print("renaming",fi, name+"_"+str(idx+1)+'.fits')
        os.rename(fi, name+"_"+str(idx)+'.fits')
    #except:
     #   print("WARNING: No product!")