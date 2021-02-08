#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:20:57 2020
@author: Ame
"""

import os
import pickle
from mat_autoPipeline import mat_autoPipeline
import argparse
import sys
#import datetime 
#import time
from tqdm import tqdm
import gc
from  mat_readReductionLog import mat_readReductionLog
#import psutil
    
_status=["NOT STARTED  ","PROCESSING...","PROCESSED  ","POSTPROCESSING...","COMPLETED    "]

class mat_dataReduction:

    def __init__(self):
        self.files=[]
        self.dirCalib=None
        self.dirResult=None
        self.paramL={}
        self.paramN={}
        self.postCmd=[]
        self.maxIter=4
        self.nbCore=1
        self.skipN=False
        self.skipL=False
        self.status=0
        self.filename=None
        self.logfile=None
        
    def addFilesFromCsh(self,filename):
        f=open(filename,"r")
        txt=f.read()
        f.close()
        
        cmd=txt.split(";")
        
        red=cmd[0]
        l=red.split('"')[1]
        files=l.replace("'","").replace("[","").replace("]","").replace(" ","").split(",")
        self.files.extend(files)
        
        return red
        
    def createScript(self,filename=None):
        txt="mat_autoPipeline.py \"{0}\" ".format(self.files)
        
        if self.dirResult!=None:
            txt+="-dirResult={0} ".format(self.dirResult)
        
        if self.dirCalib!=None:
            txt+="-dirCalib={0} ".format(self.dirCalib)
            
        if self.skipL==True:
            txt+="--skipL "
            
        if self.skipN==True:
            txt+="--skipN "       
            
        txt+="--maxIter={0} ".format(self.maxIter)
        txt+="--nbCore={0} ".format(self.nbCore)        
        
        if len(self.paramL)!=0:
            txt+="--paramL="
            for parami in self.paramL:
                txt+="/{0}={1}".format(parami,self.paramL[parami])
            txt+=" "

        if len(self.paramN)!=0:
            txt+="--paramN="
            for parami in self.paramN:
                txt+="/{0}={1}".format(parami,self.paramN[parami]) 
            txt+=" "
        
        if len(self.postCmd)!=0:
            txt+="; cd {0} ".format(self.dirResult)
            for cmdi in self.postCmd:
                txt+="; {0}".format(cmdi)
        
        if filename!=None:
            f=open(filename,"w")
            f.write(txt)
            f.close()
        return txt
        

    def save(self,filename=None):
        if filename==None:
            if self.filename!=None:
                filename= self.filename
            else:
                print("Error: filename not given")
                return
        f=open(filename,"wb")
        f.write(pickle.dumps(self.__dict__, protocol=2))
        f.close()
        self.filename=filename

    def load(self,filename):
        f=open(filename,"rb")
        dataPickle = f.read()
        f.close()
        self.__dict__ = pickle.loads(dataPickle)
        f.close()
        self.filename=filename
        
    def process(self):

        paramL="".join(["/{0}={1}".format(parami,self.paramL[parami]) for parami in self.paramL])
        paramN="".join(["/{0}={1}".format(parami,self.paramN[parami]) for parami in self.paramN])
  
        if self.status==0:
            self.status=1
            self.save()
            
            mat_autoPipeline(dirRaw="{0}".format(self.files),dirResult=self.dirResult,dirCalib=self.dirCalib,
                             nbCore=self.nbCore,paramL=paramL,paramN=paramN,
                             overwrite=1,maxIter=self.maxIter,skipL=self.skipL,skipN=self.skipN,resol="")
            gc.collect()
            self.status=2
            self.save()
        else:
            print("Error: Data reduction status is {0}.\nReset status to 0 to allow new data reduction".format(self.status))
      
    def postprocess(self):    
        if self.status==2:
            self.status=3
            self.save()                    
            if len(self.postCmd)!=0:
                txt="cd {0} ".format(self.dirResult)
                for cmdi in self.postCmd:
                    txt+="; {0}".format(cmdi)
                os.system(txt)          
            self.status=4
            self.save()
        elif self.status==3: 
            print("Error: Already performing post-processing")
        elif self.status==4: 
            print("Error: Post-processing already done.\n Set status to 2 if you want to perform it again")            
        elif self.status<2: 
            print("Error: Post-processing cannot be donne as data not reduced yet")    
      
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Launch Matisse batched data reduction from list of pkl files')
    
    parser.add_argument('action', metavar='action', type=str, help='Action : READ, STATUS, RUN')
    parser.add_argument('dir', metavar='dir', type=str, help='The directory where the pickle files are')
   
    parser.add_argument('--clean', default=0,  \
    help='Clean status', action='store_true')
    
    parser.add_argument('--process', default=1,  \
    help='do data processing step')  
    
    parser.add_argument('--postprocess', default=1,  \
    help='do post-processing step')    
        
    parser.add_argument('--nbCore', default=-1, type=int, \
    help='force number of Core used')    
        
        
    parser.add_argument('--logfile', default="",  \
    help='specify log file')    
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_ashowTransFunc.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        sys.exit(0)
        
        
    pklfiles=[os.path.abspath(os.path.join(args.dir,fi)) for fi in os.listdir(args.dir) if ".pkl" in fi]
    
    print("{0} pkl files founded in directory".format(len(pklfiles)))
   
    nred=0
    if args.logfile!="":
        n,ntot,nerr=mat_readReductionLog(args.logfile,ret=True)
        nred=len(n)
     
    process = psutil.Process(os.getpid())
    for i,fi in enumerate(pklfiles):
        #mem=process.memory_info().rss/1e6
        #print("Memory used : {0}Mo".format(mem))
        r=mat_dataReduction()    
        r.load(fi)
        if args.action=="RUN":
            print("Processing {0}".format(fi)) 
            if args.clean==True:
                r.status=0
            try:
                print("Creating dir {0}".format(r.dirResult))
                os.mkdir(r.dirResult)
            except:
                print("Dir {0} already exist".format(r.dirResult))
            if args.process==True:
                if args.nbCore!=-1:
                    r.nbCore=args.nbCore
                r.process()
            if args.postprocess==True:
                r.postprocess()
        elif args.action=="STATUS":
            txt=""  
            if args.logfile!="" and i<nred:
                if nerr[i]!=0:
                    txt="\t\t{0}/{1}\t => {2} error".format(n[i],ntot[i],nerr[i])
                else:
                    txt="\t\t{0}/{1}".format(n[i],ntot[i])
            print("{0} : {1}{2}".format(fi,_status[r.status],txt)) 
            
        elif args.action=="READ":
            print("**********************************************************")
            print("Reading : {0}".format(fi))
            print("{0} files to process".format(len(r.files)))
            print("dirResult={0}".format(r.dirResult))
            print("dirCalib={0}".format(r.dirCalib))
            print("paramL={0}".format(r.paramL))
            print("paramN={0}".format(r.paramN))
            print("skipL={0} skipN={1}".format(r.skipL,r.skipN))
            print("nbCore={0}".format(r.nbCore))
            print("postprocess={0}".format(r.postCmd))

            
            
        
        
 
   
    
    