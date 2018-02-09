#!/usr/bin/env python
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
    print "    Recipes",item[2],"running..."
    val=item[1].split("=")
    os.system("cd "+val[1]+";"+cmd+" > "+out+" 2> "+err)



repRaw=""
repArchive=""
repResult=""
tplidsel=""
tplstartsel=""
nbCore=0
listArg=sys.argv
for elt in listArg:
    if ('--help' in elt):
        print "Usage: python automaticPipeline.py --dirRaw=RawDataPath [--dirCalib=CalibrationMapPath] [--dirResult=ResultsPath] [--nbCore=NumberOfCores]"
        sys.exit(0)

for elt in listArg:
    if ('--dirRaw' in elt):
        item=elt.split('=')
        repRaw=item[1]
    if ('--dirCalib' in elt):
        item=elt.split('=')
        repArchive=item[1]
    if ('--dirResult' in elt):
        item=elt.split('=')
        repResult=item[1]
    if ('--nbCore' in elt):
        item=elt.split('=')
        nbCore=int(item[1])
    if ('--tplID' in elt):
        item=elt.split('=')
        tplidsel=item[1]
    if ('--tplSTART' in elt):
        item=elt.split('=')
        tplstartsel=item[1]

print " "
print "-----------------------------------------------------------------------------------------"
if (repRaw == ""):
    print "ERROR : You have to specifiy a Raw Data Directory with --repRaw=DirectoryPath"
    sys.exit(0)
else:
    print '%-40s' % ("Raw Data Directory:",),repRaw
if (repArchive==""):
    repArchive="/data/CalibMap"
    print "Info: Calibration Directory not specified. We used the default directory"
print '%-40s' % ("Calibration Directory:",),repArchive
if (repResult==""):
        repResult=os.getcwd()
        print "Info : Results Directory not specified. We use current directory"
print '%-40s' % ("Results Directory:",),repResult
if (nbCore==0):
    nbCore=8
    print "Info : Number of Cores not specified. We use 8 cores"
print '%-40s' % ("Number of Cores:",),nbCore    
print "-----------------------------------------------------------------------------------------"


#repRaw='/home/pbe/RawAutomatic'#'/data-matisse/TransFunc_LowFlux/OTHER'#'/home/pbe/RawAutomatic'
#repResult='/home/pbe/ResultsAutomatic'
#repArchive='/home/pbe/ArchiveAutomatic'#'/data-matisse/TransFunc_LowFlux/COSMETIC'#'/home/pbe/ArchiveAutomatic'    
listRaw = [os.path.join(repRaw, f) for f in os.listdir(repRaw) if os.path.isfile(os.path.join(repRaw, f)) and f[-5:] == '.fits' and f[0:8]=="MATISSE_"]
if (repArchive != ""):
    listArchive = [os.path.join(repArchive, f) for f in os.listdir(repArchive) if os.path.isfile(os.path.join(repArchive, f)) and f[-5:] == '.fits']
else:
    listArchive =[]

# Sort lisRaw
listRawSorted=[]
for filename in listRaw:
    hdu=fits.open(filename)
    if ('HIERARCH ESO TPL START' in hdu[0].header and 'HIERARCH ESO DET CHIP NAME' in hdu[0].header) :
        tplid=hdu[0].header['HIERARCH ESO TPL ID']
        tplstart=hdu[0].header['HIERARCH ESO TPL START']
        if (tplidsel != "" and tplstartsel !=""):
            if (tplid == tplidsel and tplstart==tplstartsel):
                listRawSorted.append(filename)
        if (tplidsel != "" and tplstartsel ==""):
            if (tplid == tplidsel):
                listRawSorted.append(filename)
        if (tplidsel == "" and tplstartsel !=""):
            if (tplstart == tplstartsel):
                listRawSorted.append(filename)
        if (tplidsel == "" and tplstartsel ==""):
               listRawSorted.append(filename)
    hdu.close()
listRaw=listRawSorted

    
# Determination of the number of Reduction Blocks
keyTplStart=[]
listIterNumber=[]
for filename in listRaw:
    hdu=fits.open(filename)
    temp=hdu[0].header['HIERARCH ESO TPL START']+'.'+hdu[0].header['HIERARCH ESO DET CHIP NAME']
    keyTplStart.append(temp)
    hdu.close()
keyTplStart=list(set(keyTplStart))                       
for elt in keyTplStart:
    listIterNumber.append(0)

iterNumber=0
while True:
    iterNumber+=1
    print ""
    print "Iteration ",iterNumber
    print "-----------------------"
    if (iterNumber > 1):
        listIter=[]
        for iter in range(iterNumber-1):
            repIterPrev=repResult+'/Iter'+str(iter+1)
            listRepIter= [os.path.join(repIterPrev, f) for f in os.listdir(repIterPrev) if os.path.isdir(os.path.join(repIterPrev, f))]
            for elt in listRepIter:
                listIter=listIter+[os.path.join(elt, f) for f in os.listdir(elt) if os.path.isfile(os.path.join(elt, f)) and f[-5:] == '.fits']

    listRedBlocks=[]
# Reduction Blocks List Initialization  
    cpt=0
    for elt in keyTplStart:
        listRedBlocks.append({"action":" ","recipes":" ","param":" ","input":[],"calib":[],"status":0,"tplstart":" ","iter":listIterNumber[cpt]})
        cpt+=1
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
    cptStatusOne=0
    cptStatusZero=0
    cptToProcess=0
    cpt=0
    for elt in listRedBlocks:
        if (elt["status"]==1):
            cptStatusOne+=1
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
                        listIterNumber[cpt]=iterNumber
                        elt["iter"]=iterNumber
                        cptToProcess+=1
                        listCmdEsorex.append(cmd)
                else:
                    listIterNumber[cpt]=iterNumber
                    elt["iter"]=iterNumber
                    cptToProcess+=1
                    listCmdEsorex.append(cmd)
            else:
                listIterNumber[cpt]=iterNumber
                elt["iter"]=iterNumber
                cptToProcess+=1
                listCmdEsorex.append(cmd)
        else:
            cptStatusZero+=1
        cpt+=1
    print '%-40s' % ("Reduction Blocks to process:",),cptToProcess
                
    if (listCmdEsorex == []):
        print " "
        print "No more iteration to do"
        print "-----------------------"
        print " "
        print "Processing summary:"
        print " "
        for elt in listRedBlocks:
            if (elt["status"] == 1):
                msg="Processing done at iteration "+str(elt["iter"])
            else:
                if (elt["action"]=="NO-ACTION"):
                    msg="Data not taken into account by the Pipeline"
                else:
                    msg="Reduction Block not processed - Missing calibration"
            tplstart,detector = elt["tplstart"].split('.')
            print '%-24s' % (tplstart,),'%-14s' % (detector,),'%-30s' % (elt["action"],),msg

        break
    else:
# Create a process pool with a maximum of 10 worker processes
        pool = Pool(processes=nbCore)
# Map our function to a data set - number 1 through 20
        pool.map(runEsorex, listCmdEsorex)
    print '%-40s' % ("Reduction Blocks processed:",),cptStatusOne
    print '%-40s' % ("Reduction Blocks not processed:",),cptStatusZero
