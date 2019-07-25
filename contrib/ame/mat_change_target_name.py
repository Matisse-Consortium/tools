#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.io import fits
import os
import glob
import sys


arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_change_target_name.py script to change target name on MATISSE raw files"
    print "Usage : dir currentname newname"

else:   
    dir=arg[1]
    currentname=arg[2].replace('"','').strip()
    newname=arg[3].replace('"','').strip()

    

os.chdir(dir)
for file in glob.glob("*.fits"): 
    print(file)   
    d=fits.open(file,mode='update')
    try:
	print(d[0].header['ESO OBS TARG NAME'].strip())
    except:
	pass
    if d[0].header['ESO OBS TARG NAME'].strip()==currentname:    
        print("changing tarname of {0} to {1}".format(file,newname))
        d[0].header['ESO OBS TARG NAME']= newname
        d.flush()
    d.close()



