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
    if (d['OI_TARGET'].data['TARGET'][0]).strip()==tarname:    
        print("changing type of {0} to {1}".format(file,tartype))
        d['OI_TARGET'].data['CATEGORY'][0]= tartype
        d.flush()
    d.close()
