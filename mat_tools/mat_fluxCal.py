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
    parser.add_argument('--scifile', default="",  \
    help='Uncalibrated oifits file of the science target.')

    #--------------------------------------------------------------------------
    parser.add_argument('--calfile', default="",  \
    help='Oifits file of the spectrophotometric calibrator.')
        
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
        print(" mat_fluxcal dir --scifile=sci.fits --calfile=cal.fits --mode='flux' --airmassCorr\n")       
        sys.exit(0)

    #-----------------------------------------
    #Path to the calibrators spectra databases 
    #-----------------------------------------
    a=imp.find_module("libFluxCal")
    dir_caldatabases=os.path.dirname(a[1])+'/calib_spec_databases'
    #dir_caldatabases='/data/users/ama/dev_python/tools/mat_tools/calib_spec_databases'

    #-----------------
    # Flux calibration
    #-----------------

    args.dir_oifits = os.path.abspath(args.dir_oifits)+"/"
    if args.mode == 'flux':
        outputfile=args.scifile.split(".")[0]+'_calflux.fits'
    elif args.mode == 'corrflux':
        outputfile=args.scifile.split(".")[0]+'_calcorrflux.fits'
        
    print('Sci = {0}'.format(args.dir_oifits+args.scifile))
    print('Cal = {0}'.format(args.dir_oifits+args.calfile))
    print('Cal database = {0}'.format(dir_caldatabases))
        
    fluxcal(args.dir_oifits+args.scifile, args.dir_oifits+args.calfile, args.dir_oifits+outputfile, dir_caldatabases, mode='flux',output_fig_dir='',match_radius=25.0,do_airmass_correction=args.airmassCorr)


