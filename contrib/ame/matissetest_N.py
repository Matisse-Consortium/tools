# -*- coding: utf-8 -*-
"""
Created on Mon Aug 01 14:17:54 2016

@author: ame
"""


from matissecubedata import *
import glob,os

plt.figure();

#final size of the interpolated + zero-padded images
band  = 2
dimfx = 1024
dimfy = 1024

# separation between 2 peaks in pixels
P0    = 4.29e-2 

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
y_in   = [230, 424, 471, 578, 695, 897, 950]
lam_in = [8.68, 9.86, 10.13, 10.8, 11.53, 12.7, 13.09]
y0     = 53
dimy   = 920
lam    = MatisseLamFit(y_in,lam_in,y0,dimy)



dir0   = '/home/fmillour/MATISSE/TESTS_NICE/OPD/'
dirfig = dir0+'/RapportOPDs/'

peakNum  = 1
dirdata  = dir0+'data/CPN1/'
darkFile = dirdata+'CPN1_CPN2_IP57_dark.fits'

#peakNum  = 3
#dirdata  = dir0+'data/CPN2/'
#darkFile = dirdata+'CPN1_CPN2_IP57_dark.fits'

peakNum  = 3
dirdata  = dir0+'data/CPN3/'
darkFile = dirdata+'CPN2_CPN3_IP35_dark.fits'

#peakNum  = 2
#dirdata  = dir0+'data/CPN4/'
#darkFile = dirdata+'CPN4_CPN3_IP13_dark.fits'


#recompute=True
recompute=False

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
    mat=MatisseMotorCalib(diri,lam,dimfx,dimfy,P0,peakNum,band, darkFile)

    if (os.path.isfile(mat.fileDown) and os.path.isfile(mat.fileUp) and not(recompute)):
        print("Datafiles already exists, skipping analysis")
    else:
        print("Analysing data...")
        mat.computeAll()
        mat.save()

    filefig=mat.plot(dirfig)
    print("Plot save in {0}".format(filefig))