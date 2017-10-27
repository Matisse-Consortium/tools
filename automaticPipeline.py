#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
  $Id:  $

  This file is part of the Matisse pipeline GUI series
  Copyright (C) 2017- Observatoire de la Côte d'Azur

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

  Created on Fri Oct 27

  $Author:  $
  $Date: $
  $Revision:  $
  $Name:  $
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import os
import shutil
import subprocess
import filecmp
from libAutoPipeline import matisseType,matisseAction,matisseRecipes,matisseCalib
from multiprocessing.pool import Pool
from functools import partial


def runEsorex(cmd):
    item=cmd.split()
    out=item[-1]+".log"
    err=item[-1]+".err"
    print cmd+" > "+out+" 2> "+err
    os.system(cmd+" > "+out+" 2> "+err)
    
repRaw='/data-matisse/TransFunc_LowFlux/OTHER'#'/home/pbe/RawAutomatic'
repResult='/home/pbe/ResultsAutomatic'
repArchive='/data-matisse/TransFunc_LowFlux/COSMETIC'#'/home/pbe/ArchiveAutomatic'    
listRaw = [os.path.join(repRaw, f) for f in os.listdir(repRaw) if os.path.isfile(os.path.join(repRaw, f)) and f[-5:] == '.fits']
listArchive = [os.path.join(repArchive, f) for f in os.listdir(repArchive) if os.path.isfile(os.path.join(repArchive, f)) and f[-5:] == '.fits']

# Determination of the number of Reduction Blocks
keyTplStart=[]

for filename in listRaw:
    hdu=fits.open(filename)
    temp=hdu[0].header['HIERARCH ESO TPL START']+'.'+hdu[0].header['HIERARCH ESO DET CHIP NAME']
    keyTplStart.append(temp)
    hdu.close()
keyTplStart=list(set(keyTplStart))                       


iterNumber=0
while True:
    iterNumber+=1
    print "Iteration ",iterNumber
    print " "
    if (iterNumber > 1):
        listIter=[]
        for iter in range(iterNumber-1):
            repIterPrev='/home/pbe/ResultsAutomatic/Iter'+str(iter+1)
            listRepIter= [os.path.join(repIterPrev, f) for f in os.listdir(repIterPrev) if os.path.isdir(os.path.join(repIterPrev, f))]
            for elt in listRepIter:
                listIter=listIter+[os.path.join(elt, f) for f in os.listdir(elt) if os.path.isfile(os.path.join(elt, f)) and f[-5:] == '.fits']

    listRedBlocks=[]
# Reduction Blocks List Initialization   
    for elt in keyTplStart:
           listRedBlocks.append({"action":" ","recipes":" ","param":" ","input":[],"calib":[],"status":0,"tplstart":" "})

# Fill the list of raw data in the Reduction Blocks List
    for filename in listRaw:
        hdu=fits.open(filename)
        stri=hdu[0].header['HIERARCH ESO TPL START']+'.'+hdu[0].header['HIERARCH ESO DET CHIP NAME']
        tag=matisseType(hdu[0].header)
        hdu.close()
        listRedBlocks[keyTplStart.index(stri)]["input"].append([filename,tag])

# Fill the list of actions,recipes,params in the Reduction Blocks List
    for elt in listRedBlocks:
        hdu=fits.open(elt["input"][0][0])
        keyTplStartCurrent=hdu[0].header['HIERARCH ESO TPL START']+'.'+hdu[0].header['HIERARCH ESO DET CHIP NAME']
        action=matisseAction(hdu[0].header,elt["input"][0][1])
        recipes,param=matisseRecipes(action)
        hdu.close()
        elt["action"]=action
        elt["recipes"]=recipes
        elt["param"]=param
        elt["tplstart"]=keyTplStartCurrent

# Fill the list of calib in the Reduction Blocks List from repArchive
    for elt in listRedBlocks:
        hdu=fits.open(elt["input"][0][0])
        calib,status=matisseCalib(hdu[0].header,elt["action"],listArchive,elt['calib'])
        hdu.close()
        elt["calib"]=calib
        elt["status"]=status
        
# Fill the list of calib in the Reduction Blocks List from repResult Iter i-1
    if (iterNumber > 1):
        for elt in listRedBlocks:
            hdu=fits.open(elt["input"][0][0])
            calib,status=matisseCalib(hdu[0].header,elt["action"],listIter,elt['calib'])
            hdu.close()
            elt["calib"]=calib
            elt["status"]=status
            
# Create the SOF files
    repIter=repResult+"/Iter"+str(iterNumber)
    if (os.path.isdir(repIter) == True):
        shutil.rmtree(repIter)
        os.mkdir(repIter)
    else:
        os.mkdir(repIter)

    listCmdEsorex=[]
    for elt in listRedBlocks:
        print elt["status"],elt["action"],elt["tplstart"]
#        print elt['input']
#        print elt['calib']
#        print '**********************************************'
        if (elt["status"]==1):
            sofname=repIter+"/"+elt["recipes"]+"."+elt["tplstart"]+".sof"       
            outputDir=repIter+"/"+elt["recipes"]+"."+elt["tplstart"]+".rb"
            os.mkdir(outputDir)
            fp=open(sofname,'w')
            for frame,tag in elt['input']:
                fp.write(frame+" "+tag+"\n")
            for frame,tag in elt['calib']:
                fp.write(frame+" "+tag+"\n")
            fp.close()
            cmd="esorex --output-dir="+outputDir+" "+elt['recipes']+" "+elt['param']+" "+sofname
            if (iterNumber > 1):
                sofnamePrev=repIterPrev+"/"+elt["recipes"]+"."+elt["tplstart"]+".sof"
                if (os.path.exists(sofnamePrev)):
                    if (filecmp.cmp(sofname,sofnamePrev)):
                        #outputDirPrev=repIterPrev+"/"+elt["tplstart"]+".rb"
                        #print "copy "+outputDirPrev+" in "+outputDir
                        #os.system("cp "+outputDirPrev+"/* "+outputDir+"/.")
                        shutil.rmtree(outputDir)
                    else:
                        listCmdEsorex.append(cmd)
                else:
                    listCmdEsorex.append(cmd)
            else:
                listCmdEsorex.append(cmd)
                
    if (listCmdEsorex == []):
        print "No more iteration to do"
        break
    else:
# Create a process pool with a maximum of 10 worker processes
        pool = Pool(processes=4)
# Map our function to a data set - number 1 through 20
        pool.map(runEsorex, listCmdEsorex)




