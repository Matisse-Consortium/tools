#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 09:28:53 2021

@author: ame
"""

import os
import argparse
import sys
from astropy.io import fits
import numpy as np
from libFluxCalib import mat_calibrateTotalFlux,mat_CreateTFinFlux
import shutil


def _findCalIndex(i0,tts,way):
    l=len(tts)
    i=i0
    stop=False
    while(not(stop)):
        i+=way
        if ((i<0) | (i>=l)):
            stop=True
            i=-1
        else:
            if (tts[i]=="CALIB_RAW_INT"):
                stop=True
    return i

def mat_autoCalibFlux(dirUncalibrated,dirCalibrated,dirout="_FLXCALIBRATED",avgTel=True,verbose=True):

    filenames=[os.path.join(dirUncalibrated,fi) for fi in os.listdir(dirUncalibrated) if ".fits" in fi]
    nfiles=len(filenames)

    headers=[fits.getheader(fi) for fi in filenames]
    BCD_list=[h["HIERARCH ESO DET BCD STATE"] for h in headers]
    chop_list=[h["HIERARCH ESO ISS CHOP ST"]=="T" for h in headers]
    #tplstart_list=[h["HIERARCH ESO TPL START"] for h in headers]
    grouped_filenames=[[] for i in range(8)]
    grouped_headers=[[] for i in range(8)]

    for ifiles in range(nfiles):
        igroup=(chop_list[ifiles])*4+BCD_list[ifiles]
        grouped_filenames[igroup].append(filenames[ifiles])
        grouped_headers[igroup].append(headers[ifiles])
     
    BCDs=["OUT-OUT","IN-IN  ","IN-OUT ","OUT-IN "]    
    
    
    sci_filenames=[os.path.join(dirCalibrated,fi) for fi in os.listdir(dirCalibrated) if ".fits" in fi]
    #sci_nfiles=len(sci_filenames)
    sci_headers=[fits.getheader(fi) for fi in sci_filenames]
    sci_BCD_list=[h["HIERARCH ESO DET BCD STATE"] for h in sci_headers]
    sci_chop_list=[h["HIERARCH ESO ISS CHOP ST"]=="T" for h in sci_headers]
    sci_tplstart_list=[h["HIERARCH ESO TPL START"] for h in sci_headers]
    
    
    for (igroup,grouped_filenamei) in enumerate(grouped_filenames):
        nfilei=len(grouped_filenamei)
        chopi="No Chopping" if igroup<4 else "   Chopping"        
        BCDi=BCDs[igroup % 4]
        if nfilei!=0:
            if verbose==True :
                print("\n\n**************************************************************************")  
                print("************************** {0} BCD {1} ************************".format(chopi,BCDi)) 
                print("**************************************************************************")              
                print("{0} files found".format(nfilei))
        
            names=[]
            for ifile in range(nfilei):
                hi=grouped_headers[igroup][ifile]
                try:
                    names.append(hi["HIERARCH ESO PRO JSDC NAME"])
                except:
                    names.append(hi["HIERARCH ESO OBS TARG NAME"])
            names=np.array(names)       
            MJDs=np.array([hi["MJD-OBS"] for hi in grouped_headers[igroup]]).astype(float)  
            tartypes=np.array([hi["HIERARCH ESO PRO CATG"] for hi in grouped_headers[igroup]])   
            filesi=np.array(grouped_filenamei)
            headersi=np.array(grouped_headers[igroup])
           
            
            #sorting infos by MJD
            idx=np.argsort(MJDs)       
            names=names[idx]
            MJDs=MJDs[idx]
            tartypes=tartypes[idx]
            filesi=filesi[idx]
            headersi=headersi[idx]
            
            if verbose==True :
                print(filesi)
                print(MJDs)
                print(names)
            
            if verbose==True :
                print("\n---------------------------------------------------------------------------")
                print("------------------------------ Creating TFFs ------------------------------")  
            #creating the TFF for the CALIB files (modif saved on TF_FLUX table)
            for ifile in range(nfilei):
                if tartypes[ifile]=="CALIB_RAW_INT":
                    if verbose:
                        print("------------------------------------------------------------------------\n")
                        print("Computing TF in Flux for {0}".format(os.path.basename(filesi[ifile])))
                    
                    mat_CreateTFinFlux(filesi[ifile],fluxModel="bb",fit=True,Teff="Gaia",
                            diam="JSDC",filename=None,filenameThFlux=None,query="cds",
                            verbose=verbose,save=True,overwrite=True)
            if verbose==True :
                print("---------------------------------------------------------------------------")
                print("-------------------------- Calibrating Targets ----------------------------")  
            try:
                os.mkdir(dirout)
            except:
                for fi in os.listdir(dirout):
                    os.remove(os.path.join(dirout,fi)) 
            #Calibrating TARGET files
            for ifile in range(nfilei):   
                if verbose==True :
                    print("------------------------------------------------------------------------")
                    print("Checking {0} idx={1}".format(filesi[ifile],ifile))
                if tartypes[ifile]=="TARGET_RAW_INT": 
                    if verbose==True :
                        print("Calibrating Flux for {0} idx={1}".format(filesi[ifile],ifile))
                    #determining the CALIB files to use and their weights
                    im=_findCalIndex(ifile,tartypes,-1)
                    ip=_findCalIndex(ifile,tartypes,1)                
                    if (ip!=-1) & (im!=-1):
                        weights=np.array([MJDs[ip]-MJDs[ifile],MJDs[ifile]-MJDs[im]])
                        weights=weights/np.sum(weights) 
                        ids=np.array([im,ip])
                    elif ip==-1:
                        weights=np.array([1])
                        ids=np.array([im])
                    else: #im==-1
                        weights=np.array([1])
                        ids=np.array([ip])
                     
                     
                    if verbose:
                        weigthTxt=["{0:.1f}".format(wi) for wi in weights]
                        if verbose==True :
                            print("Using calibrators at indices {0} with weights={1}".format(ids,weigthTxt))
                    #find the calibrated sci file corresponding to the uncalibrated one
                    hi=headersi[ifile]
                    tplstart=hi["HIERARCH ESO TPL START"]

                    w1=(np.array(sci_tplstart_list)==tplstart)
                    w2=(np.array(sci_BCD_list)==(igroup%4))
                    w3=(np.array(sci_chop_list)==(igroup//4))
                    """
                    print(sci_tplstart_list)
                    print(sci_BCD_list)
                    print(sci_chop_list)
                    print(w1)
                    print(w2)
                    print(w3)
                    print(tplstart)
                    """
                    
                    
                    
                    ww=np.where(w1 & w2 & w3)[0]

                    if len(ww)==0:
                        if verbose==True :
                            print("Skipping {0} idx={1}".format(scifname,ifile))
                    else:
                        ww=ww[0]
                        scifname=sci_filenames[ww]
                        calfnames=[filesi[idi] for idi in ids]
                        if verbose:
                            print("Calibrating Flux for {0} idx={1}".format(scifname,ifile))
                        mat_calibrateTotalFlux(scifname,calfnames,outdir=dirout,weight=weights,verbose=verbose,avgTel=avgTel)
                else:
                    if verbose:
                        print("=>Calibrator")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calibrate Targets total flux\n '
              'The theoretical flux of the calibrators is computed using a blackbody'
              'fitted on the simbad & vizier fluxes and magnitudes\n'
              'The Transfer-Function in flux (TFF) for each calibrator is stored in the TF_FLUX table\n'
              'The TFF is linearly interpolated during the night.'
              ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('dirUncalibrated', type=str, help='Directory of the uncalibrated data')
    parser.add_argument('dirCalibrated', type=str, help='Directory of the calibrated data')
   
    parser.add_argument('--dirOut', default="_FLXCALIBRATED",  \
    help='directory where the flux-calibrated fits files will be stored')
 
    parser.add_argument('--sepTels', type=str,default="FALSE", \
    help='Separate the fluxes on the 4 telescopes (TRUE/FALSE)')
        
    parser.add_argument('--verbose', default=0,  \
    help='print information on the flux calibration', action='store_true')

    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_autoCalibFlux.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        sys.exit(0)
        
    avgTel=True if args.sepTels=="FALSE" else False
    
    print(avgTel)
    mat_autoCalibFlux(args.dirUncalibrated,args.dirCalibrated,dirout=args.dirOut,
                      avgTel=avgTel,verbose=args.verbose)