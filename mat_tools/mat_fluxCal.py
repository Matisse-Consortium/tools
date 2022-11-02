#!/usr/bin/env python3
# -*- coding: utf-8 -*-
########################################################################

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from shutil import copyfile
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.special import j0,j1
from scipy.interpolate import interp1d
import math
import os
from astroquery.simbad import Simbad
from numpy.polynomial.polynomial import polyval
from astropy.convolution import Gaussian1DKernel,Box1DKernel,convolve
import scipy.stats
from libFluxCal import *
import argparse
import sys
import imp
import shutil
from operator import itemgetter


###################### Main ################################################################

if __name__ == '__main__':
    print("Starting...")

    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Spectrophotometric calibration of MATISSE total and correlated spectra.')

    #--------------------------------------------------------------------------
    parser.add_argument('dir_oifits', default="",  \
    help='The path to the directory containing the MATISSE oifits files of the science target and of the spectrophotometric calibrator')

    #--------------------------------------------------------------------------
    #parser.add_argument('dir_caldatabases', default="",  \
    #help='The path to the directory containing the calibrator synthetic spectra databases')

    #--------------------------------------------------------------------------
    parser.add_argument('--sciname', default="",  \
    help='Name of the science target (as indicated in the oifits file name).')

    #--------------------------------------------------------------------------
    parser.add_argument('--calname', default="",  \
    help='Name of the spectrophotometric calibrator (as indicated in the oifits file name).')
        
    #--------------------------------------------------------------------------
    parser.add_argument('--mode', default="flux", \
                        help='Type of MATISSE spectra you want to calibrate. "flux": total spectra or "corrflux": correlated spectra.')

    #--------------------------------------------------------------------------
    parser.add_argument('--airmassCorr',  \
                        help='Do airmass correction between the science and the calibrator', action='store_true')
    #--------------------------------------------------------------------------

    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_fluxcal.py --help to be kind with you:\033[0m\n")
        #parser.print_help()
        print("\n This routine can produce two types of calibrated oifits files depending on the selected mode (either 'flux' or 'corrflux':\n")
        print("\n - ***_calflux.fits: only total flux is calibrated (incoherently processed oifits file expected) and stored in the OI_FLUX table (FLUXDATA column).\n")
        print("\n - ***_calcorrflux.fits: only correlated fluxes are calibrated (coherently processed oifits file expected) and stored in the OI_VIS table (VISAMP column).\n")
        print("\n Example of calibration of the total flux of a MATISSE science oifits files with airmass correction:\n") 
        print(" mat_fluxcal dir --sciname='sci' --calname='cal' --mode='flux' --airmassCorr\n")       
        sys.exit(0)

    #-----------------------------------------
    #Path to the calibrators spectra databases 
    #-----------------------------------------
    a=imp.find_module("libFluxCal")
    dir_caldatabases=os.path.dirname(a[1])+'/calib_spec_databases'
    #dir_caldatabases='/data/users/ama/dev_python/tools/mat_tools/calib_spec_databases'

    #---------------------
    # Oifits files sorting
    #---------------------
    
    args.dir_oifits = os.path.abspath(args.dir_oifits)+"/"
    scifiles=glob.glob(args.dir_oifits+'*'+args.sciname+'*Chop.fits')
    calfiles=glob.glob(args.dir_oifits+'*'+args.calname+'*Chop.fits')
    nfiles_sci=np.size(scifiles)
    nfiles_cal=np.size(calfiles)
    list_of_dicts_sci=[]
    list_of_dicts_cal=[]
    for i in range(nfiles_sci):
        hdul_sci=fits.open(scifiles[i])
        hdr_sci=hdul_sci[0].header
        dateobs_sci=hdr_sci['MJD-OBS']
        bcd_pos1=hdr_sci['HIERARCH ESO INS BCD1 ID']
        bcd_pos2=hdr_sci['HIERARCH ESO INS BCD2 ID']
        chop_status=hdr_sci['HIERARCH ESO DET CHOP ST']
        tpl_start=hdr_sci['HIERARCH ESO TPL START']
        filename=scifiles[i].split("/")[-1]
        dic_sci = {'DATE_OBS':dateobs_sci}
        dic_sci['TPL_START']=tpl_start
        dic_sci['FILENAME']=filename
        dic_sci['BCD']=bcd_pos1+'-'+bcd_pos2
        dic_sci['CHOP_ST']=chop_status
        list_of_dicts_sci.append(dic_sci)

    for i in range(nfiles_cal):
        hdul_cal=fits.open(calfiles[i])
        hdr_cal=hdul_cal[0].header
        dateobs_cal=hdr_cal['MJD-OBS']
        bcd_pos1=hdr_cal['HIERARCH ESO INS BCD1 ID']
        bcd_pos2=hdr_cal['HIERARCH ESO INS BCD2 ID']
        chop_status=hdr_cal['HIERARCH ESO DET CHOP ST']
        tpl_start=hdr_cal['HIERARCH ESO TPL START']
        filename=calfiles[i].split("/")[-1]
        dic_cal = {'DATE_OBS':dateobs_cal}
        dic_cal['TPL_START']=tpl_start
        dic_cal['FILENAME']=filename
        dic_cal['BCD']=bcd_pos1+'-'+bcd_pos2
        dic_cal['CHOP_ST']=chop_status
        list_of_dicts_cal.append(dic_cal)

    sorted_list_of_dicts_sci = sorted(list_of_dicts_sci,key=itemgetter('DATE_OBS'))
    sorted_scifiles=[]
    for dic in sorted_list_of_dicts_sci:
        files = dic['FILENAME']
        tpl_start_sci = dic['TPL_START']
        #print('files={0}').format(files)
        sorted_scifiles.append(files)
        #sorted_tpl_start_sci.append(files)
    sorted_list_of_dicts_cal = sorted(list_of_dicts_cal,key=itemgetter('DATE_OBS'))
    sorted_calfiles=[]
    for dic in sorted_list_of_dicts_cal:
        files = dic['FILENAME']
        tpl_start_cal = dic['TPL_START']
        #print('files={0}').format(files)
        sorted_calfiles.append(files)
        #sorted_tpl_start_cal.append(files)

        
    
    #-----------------
    # Flux calibration
    #-----------------

    nfiles=np.size(sorted_scifiles)
    for i in range(nfiles):    
        print('Sci = {0}'.format(sorted_scifiles[i]))
        print('Cal = {0}'.format(sorted_calfiles[i]))
        print('Cal database = {0}'.format(dir_caldatabases))
        if args.mode == 'flux':
            outputfile=sorted_scifiles[i].split(".")[0]+'_calflux.fits'
        elif args.mode == 'corrflux':
            outputfile=sorted_scifiles[i].split(".")[0]+'_calcorrflux.fits'
        fluxcal(args.dir_oifits+sorted_scifiles[i],args.dir_oifits+sorted_calfiles[i],args.dir_oifits+outputfile, dir_caldatabases, mode='flux',output_fig_dir='',match_radius=25.0,do_airmass_correction=args.airmassCorr)


