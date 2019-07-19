#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
  $Id$

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2018- Observatoire de la Côte d'Azur

  Created in 2016
  @author: J. Isbell, F. Millour

  Automatic MATISSE calibration !

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
"""

from   subprocess import call 
import argparse
import glob
import os
import sys
import numpy as np
from astropy.io import fits


#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

if __name__ == '__main__':
    print("Starting...")
    
    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Wrapper to run the calibration steps of the MATISSE DRS on a given directory containong raw OIFITS files.')

    #--------------------------------------------------------------------------
    parser.add_argument('in_dir', metavar='in_dir', type=str, \
    help='The path to the directory containing your oifits data.', default=None)

    #--------------------------------------------------------------------------
    parser.add_argument('-o', '--out_dir', dest='out_dir', metavar='Out Dir', type=str, default='./CALIBRATED', \
    help='The path to the directory you want results to be stored in (defaults to current directory).')

    #--------------------------------------------------------------------------
    #parser.add_argument('-i', '--in_dir', dest='in_dir', metavar='Working Directory', type=str, \
    #help='The path to the directory containing your oifits data.')

    #--------------------------------------------------------------------------
    parser.add_argument('--timespan', dest='timespan', metavar='Timespan of calibs', type=str, default='.', \
    help='The time search interval for selecting calibrators around the science star')
    #--------------------------------------------------------------------------

    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_autoCalib.py --help to be kind with you:\033[0m\n")
        parser.print_help()
	sys.exit(0)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
        
    #----------------------------------------------------------------------
    #---- Make the SOF files ----------------------------------------------
    targsof = make_sof(args.in_dir, args.out_dir)
    
    #--------------------------------------------------------------------------
    #----- Run the Recipes ----------------------------------------------------
    for isof in targsof:
        print 'Running mat_cal_oifits on sof:%s'%(isof)
        call("esorex --output-dir=%s mat_cal_oifits %s"%(args.out_dir,isof), shell=True)

        name, ext = os.path.splitext(isof)
        #print(name)
        
        # Rename files
        resultFiles = glob.glob(args.out_dir+'/TARGET_CAL_INT_????.fits')
        #print(resultFiles)
        for idx,fi in enumerate(resultFiles):
            print("renaming",fi, name+"_"+str(idx+1)+'.fits')
            os.rename(fi, name+"_"+str(idx)+'.fits')
