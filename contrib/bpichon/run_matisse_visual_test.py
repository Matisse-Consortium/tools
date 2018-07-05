"""
Created on 2017-10-05
@author: bpichon for matisse
"""

#==============================================================================
# MAIN MATISSE QUALITY REPORT PROGRAM    
#==============================================================================

import os
import sys
##import traceback

try:
   import pyfits
except:
   from astropy.io import fits as pyfits
import matisse_visual
from glob import glob
#from optparse import OptionParser

# Add options
usage = """
        usage:  %prog
"""
#parser = OptionParser(usage)
#parser.add_option("-o","--overwrite", action="store_true", dest="overwrite_flag", default=False,
#                  help="Overwrite existing PDFs")

#(argoptions, args) = parser.parse_args()

filelist=[]
## If the user specifies a file name or wild cards ("*_0001.fits")
if len(sys.argv) > 1 :
    longnames = [f for files in sys.argv[1:] for f in glob(files)]
    filelist = [os.path.splitext(f)[0] for f in longnames]
## Processing of the full current directory
else :
    for file in os.listdir("."):
        if file.endswith(".fits"):
            filelist.append(os.path.splitext(file)[0])
#
#BP filelist.sort() # process the files in alphabetical order
#
for filename in filelist:
    print("\nFILE = "+filename)
    fullname = os.getcwd()+'/'+filename
    matisse_file=pyfits.open(fullname+'.fits', mode='readonly')
    header = matisse_file[0].header
    #
    try:
        if not 'MATISSE' in header['INSTRUME']:
            print (' Not for MATISSE in file : '+filename+'.fits')
            continue
    except Exception:
        print (' Bad header INTRUME in file : '+filename+'.fits')
        continue
    #
    try:
        datatype = header['HIERARCH ESO PRO TYPE']
        print (' ==> Datatype = %s'%datatype)
        if not datatype == 'REDUCED':
            print (' Not REDUCED data type in file : '+filename+'.fits')
            continue
    except Exception:
        print (' Bad header HIERARCH ESO PRO TYPE in file : '+filename+'.fits')
        continue
    #
    try:            
        datacatg = header['HIERARCH ESO PRO CATG']
        print (' ==> Datacatg = %s'%datacatg)
        if datacatg not in ('CALIB_RAW_INT','TARGET_RAW_INT'):
            print (' Not CALIB_RAW_INT or TARGET_RAW_INT for datacatg in file : '+filename+'.fits')
            continue
    except Exception:
        print (' Bad header HIERARCH ESO PRO CATG in file : '+filename+'.fits')
        continue
    #
    try:
        object = header['OBJECT']
        print (' ==> Object   = %s'%object)
    except Exception:
        print (' Bad header OBJECT in file : '+filename+'.fits')
        continue
    #
    try:
        detector = header['HIERARCH ESO DET CHIP NAME']
        print (' ==> Detector = %s'%detector)
    except Exception:
        print (' Bad header HIERARCH ESO DET CHIP NAME : '+filename+'.fits')
        continue
    #
    if detector not in ('HAWAII-2RG','AQUARIUS'):
        print (' Not a good detector in file : '+filename+'.fits')
        continue
    #
    band_L = (2.8,4.0)
    band_M = (4.5,5.0)
    band_N = (8.0,13.0)
    if detector == 'HAWAII-2RG':
        band = ("LM",(band_L,band_M))    # band, lamndaS_min, lambdaS_max
    else:
        band = ("N",(band_N,))           # band, lamnda_min, lambda_max
#    print (' ==> Band     = %s'%band[0])
#    print ( "  Band[1]  = ", band[1] )
    #
    if ( object in ("SOURCE") ) :
        lamp = True
    elif ( object in ("STD","SCIENCE") ) :
        lamp = False
    else : 
        print( " PROBLEM in Object header " )
        continue
    #
    band += (lamp,)                    # to check end of first part of plots

#    if os.path.isfile(fullname+"-OIFITS.pdf") and argoptions.overwrite_flag==False:
#        print ('  output PDF already exists, skip...')
#        continue

    oifitsdata = matisse_visual.Oifits(fullname+'.fits')
#    print ("==============================")
    matisse_visual.produce_oifits_report(oifitsdata,filename,band)
    matisse_file.close()
