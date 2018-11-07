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

pathF='/data/users/fal/RUN1C/'
#detect='HAWAI'
#tplstart='052136'
tplstart=sys.argv[1]
detect=sys.argv[2]
savedirectory=sys.argv[3]
if savedirectory != '/':
    savedirectory=savedirectory+'/'
if detect == 'AQUA':
    detecteur='IUS.rb'
if detect == 'HAWAI':
    detecteur='2RG.rb'

#os.chdir(pathF)

for r,d,f in os.walk(pathF):
    if r[22:26]=='2018':
        if r[33:37]=='Iter':
            if r[43:46]=='raw':
                if r[-6:]==detecteur:
                    if (r[68:76].replace('_','')==tplstart):
                        direct=np.append(direct,d)
                        directg=np.append(directg,r)




for f in os.walk(directg[0]):
    files=f[2]
    dossier=[]
    for j in files:
        if j[0:4]=='PHOT':
            dossier=np.append(dossier,j)

diconc={}
dicoc={}
date=[]
chopping=[]
longueur=len(dossier)
for i in range(longueur):
    h=fits.open(directg[0]+'/'+dossier[i])
    thedate=h[0].header['MJD-OBS']
    chop=h[0].header['HIERARCH ESO ISS CHOP ST']
    date=np.append(date,str(thedate))
    chopping=np.append(chopping,chop)
    if chop=='F':
        diconc[dossier[i]]=thedate
    if chop=='T':
        dicoc[dossier[i]]=thedate


orderc=OrderedDict(sorted(dicoc.items(), key=lambda x: x[1]))
danslordrec=orderc.keys()
ordernc=OrderedDict(sorted(diconc.items(), key=lambda x: x[1]))
danslordrenc=ordernc.keys()


print(danslordrenc)

chopper=danslordrenc[-1:]+danslordrec
del danslordrenc[-1:]
nonchopper=danslordrenc

pathF=directg[0]
print('pathF=',pathF)
print('chopper=',chopper)
print('nonchopper=',nonchopper)
print('tralala')


np.savez(savedirectory+'pathAndfolder.npz',pathF,nonchopper,chopper)
