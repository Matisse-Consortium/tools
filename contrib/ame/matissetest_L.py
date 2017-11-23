# -*- coding: utf-8 -*-
"""
Created on Mon Aug 01 14:17:54 2016

@author: ame
"""


from matissecubedata import *
import glob,os

plt.figure();

#final size of the interpolated + zero-padded images
band  = 1
dimfx=1024
dimfy=1024

# separation between 2 peaks in pixels
P0=3.72e-2 

#---------------------------#
# peak used for the recording
#---------------------------#
# without BCD     #
# IP57 -> peakNum=1
# IP13 -> peaknum=2
# IP35 -> peaknum=3
# IP37 -> peaknum=4
# IP15 -> peaknum=5
# IP17 -> peaknum=6

# with BCD        #
# IP75 -> peakNum=1
# IP31 -> peaknum=2
# IP17 -> peaknum=3
# IP15 -> peaknum=4
# IP37 -> peaknum=5
# IP35 -> peaknum=6


#Wavelength calibration 
y_in   = [1754,1498,1303]
lam_in = [2.87,3.23,3.5]
y0     = 1201
dimy   = 650
lam    = MatisseLamFit(y_in,lam_in,y0,dimy)


dir0   = '/home/fmillour/MATISSE/TESTS_NICE/OPD/'
dirfig = dir0+'/RapportOPDs/'

#peakNum=4
#dirdata=dir0+'data/BCOP1/'

#peakNum=5
#dirdata=dir0+'data/BCOP2/'

#peakNum=4
#dirdata=dir0+'data/BCOP3/'

peakNum=5
dirdata=dir0+'data/BCOP4/'


recompute=True

print("Working directory : {0}".format(dirdata))
dirs=os.listdir(dirdata)
dirs=[diri for diri in dirs if(os.path.isdir(dirdata+diri)) ]

DatafilesUp=glob.glob(dirdata+"*_Up.dat")
DatafilesDown=glob.glob(dirdata+"*_Down.dat")
print("{0} subdirectories found".format(len(dirs)))
print("{0} *_Up.dat files found".format(len(DatafilesUp)))
print("{0} *_Down.dat files found".format(len(DatafilesDown)))
for diri in dirs:
    diri=dirdata+diri+'/'
    print("Current Directory : {0} ".format(diri) )    
    mat=MatisseMotorCalib(diri,lam,dimfx,dimfy,P0,peakNum,band)

    if (os.path.isfile(mat.fileDown) and os.path.isfile(mat.fileUp) and not(recompute)):
        print("Datafiles already exists, skipping analyzis")
    else:
        print("Analazing data...")
        mat.computeAll()
        mat.save()

    filefig=mat.plot(dirfig)
    print("Plot save in {0}".format(filefig))