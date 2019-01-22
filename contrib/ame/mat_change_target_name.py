from astropy.io import fits
import os
import glob
import sys


arg=sys.argv

if (arg[1]=="--help" or arg[1]== "-h"):
    print "mat_change_target_catergory.py script to switch between SCI and CAL"
    print "Usage : dir currentname newname"

else:   
    dir=arg[1]
    currentname=arg[2].replace('"','').strip()
    newname=arg[3].replace('"','').strip()

    

os.chdir(dir)
for file in glob.glob("*.fits"):    
    d=fits.open(file,mode='update')
    if (d['OI_TARGET'].data['TARGET'][0]).strip()==currentname:    
        print("changing tarname of {0} to {1}".format(file,newname))
        d['OI_TARGET'].data['TARGET'][0]= newname
        d.flush()
    d.close()
