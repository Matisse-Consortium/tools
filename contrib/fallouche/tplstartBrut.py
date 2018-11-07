import os,time,stat
from datetime import datetime
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import ndimage
from mpl_toolkits import mplot3d
import fnmatch
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from astropy.io import ascii
import operator
from collections import OrderedDict


direct=[]
directg=[]

pathF='/data/RawDataMatisse/2018-09-23/'
tplstart=sys.argv[1]
detect=sys.argv[2]
savedirectory=sys.argv[3]
if savedirectory[-1:]!='/':
    savedirectory=savedirectory+'/'
os.chdir(pathF)

os.system('dfits *.fits | fitsort tpl.start det.chip.name mjd-obs date iss.chop.st | grep ' + tplstart + '|grep '+ detect +'| sort -k 2  > outcome.txt')
files=np.loadtxt('outcome.txt',dtype=str,usecols=(0,3))


longueur=np.shape(files)[0]
nonchopper=[]
chopper=[]
comptnonchopper=[]
comptchopper=[]
for i in range(longueur):
    h=fits.open(files[i,0])
    chop=h[0].header['HIERARCH ESO ISS CHOP ST']
    tart=h['IMAGING_DATA'].data['TARTYP']
    l=len(tart)
    compteur=0
    for j in range(l):
        if tart[j]=='S':
            compteur=compteur+1
    
    if compteur<l and compteur > 0:
        chopper=np.append(chopper,files[i,0])
        comptchopper=np.append(comptchopper,i)
    if compteur==0:
        nonchopper=np.append(nonchopper,files[i,0])
        comptnonchopper=np.append(comptnonchopper,i)
    



np.savez(savedirectory+'pathAndfolderBrut.npz',pathF,nonchopper,chopper)

