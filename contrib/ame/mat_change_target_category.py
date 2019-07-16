#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.io import fits
import os
import glob
import sys


arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_change_target_catergory.py script to switch between SCI and CAL"
    print "Usage : dir tarname type"

else:   
    dir=arg[1]
    tarname=arg[2].replace('"','').strip()
    tartype=arg[3]
    
    

os.chdir(dir)
for file in glob.glob("*.fits"):    
    d=fits.open(file,mode='update')
    #print("{0}=>{1}".format(file, d[0].header['ESO OBS TARG NAME'].strip()))
    if d[0].header['ESO OBS TARG NAME'].strip()==tarname:    
        print("changing type of {0} to {1}".format(file,tartype))

	if tartype=="SCI":	
       		d[0].header['ESO PRO CATG']= "TARGET_RAW_INT"
	elif tartype=="CAL":
       		d[0].header['ESO PRO CATG']= "CALIB_RAW_INT"			        
	d.flush()
    d.close()
