# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 09:59:00 2016

@author: pbe
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calibrate correlated flux of calibrators using black body law.')
    parser.add_argument('dirIn', default="",help='Input directory')
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_MergeAllOiFits.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        print("\n     Example : python mat_calCalibratorCorFlux.py .")
        sys.exit(0)

    dicTeff={"* phi01 Lup": 4200., "* del Sgr": 4350.}
        
    oifitsfiles=[args.dirIn+"/"+fi for fi in os.listdir(args.dirIn) if ".fits" in fi]

    dirOut=args.dirIn+"_CORRECTED/"
    if (os.path.isdir(dirOut)):
        shutil.rmtree(dirOut)
    os.mkdir(dirOut)
    for elt in oifitsfiles:
        hdu=fits.open(elt)
        hdr=hdu[0].header
        #chip=hdr["ESO DET CHIP NAME"]
        amptype=hdu['OI_VIS'].header['AMPTYP']
        if ("ESO PRO JSDC NAME" in hdr and amptype == "correlated flux"):
            print amptype
            wave=hdu["OI_WAVELENGTH"].data["EFF_WAVE"]
            name = hdr["ESO PRO JSDC NAME"]
            temp = dicTeff[name]
            diam_mas = hdr["ESO PRO JSDC DIAMETER"]
            print elt,name,temp
            flux=(np.pi*(diam_mas*np.pi/180/3600/1000)**2/4) * ( 2*6.62E-34* 3E8/(wave)**3 / ( np.exp(6.62E-34*3E8/wave/1.38E-23/temp) - 1 ) * 1E26)
            hdu["OI_VIS"].data["VISAMP"]*=flux
            hdu["OI_VIS"].data["VISAMPERR"]*=flux
            hdu["OI_VIS2"].data["VIS2DATA"]*=flux**2
            hdu["OI_VIS2"].data["VIS2ERR"]*=flux**2
        hdu.writeto(dirOut+os.path.basename(elt))
        hdu.close()

    #print dicTeff["* del Sgr"]
