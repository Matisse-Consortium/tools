#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.io import fits
import os
import glob
import sys


arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_change_target_catergory.py script to switch between SCIENCE and CALIB"
    print "Usage : dir tarname type"

else:   
    dir=arg[1]
    tarname=arg[2].replace('"','').strip()
    tartype=arg[3]



os.chdir(dir)
for file in glob.glob("*.fits"):
    d=fits.open(file,mode='update') 
    try :
        if d[0].header['ESO OBS TARG NAME'].strip()==tarname: 
            print("changing type of {0} to {1}".format(file,tartype))
            d[0].header['ESO DPR CATG'] = tartype
            if tartype=="SCIENCE" and d[0].header['ESO DO CATG']=="CALIB_RAW":
                d[0].header['ESO DO CATG'] = "TARGET_RAW"
            elif tartype=="CALIB" and d[0].header['ESO DO CATG']=="TARGET_RAW":
                 d[0].header['ESO DO CATG'] = "CALIB_RAW"
            d.flush()
    except:
        print("bad {0}".format(file))
    d.close()    
    
 
