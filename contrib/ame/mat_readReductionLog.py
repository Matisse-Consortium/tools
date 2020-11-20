#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 17:49:37 2020

@author: Ame
"""
import argparse
import sys
import os
import time



def removeAll(string,chars):
    for char in chars:
        string=string.replace(char,"")
    return string
    


def mat_readReductionLog(fname,showErrors=False,ret=False):
    
    """
    ctime=os.path.getctime(fname)
    mtime=os.path.getmtime(fname)
    print(time.ctime(ctime))
    print(time.ctime(mtime))
    dt=mtime-ctime
    """
    f=open(fname,mode="r")#,errors="ignore")
    data=f.read()
    
    
    
    f.close()  
    n=[]
    ntot=[]
    nerr=[]
    rec=[]
    err=[]
    red=[di for di in data.split("Reduction Blocks to process")[1:]]
    for redi in red:
        lines=redi.split("\n")
        ntoti=int(removeAll(lines[0]," :'\"),"))
        
        if ntoti!=0:  
            recipes=[]
            erri=[]
            for line in lines[1:]:
                if "Segmentation fault" in line:
                    erri.append(line)
                elif "esorex" in line:
                    recipes.append(line)
            
            n.append(len(recipes))
            ntot.append(ntoti)
            nerr.append(len(erri))
            err.append(erri)
            rec.append(recipes)
            
    nred=len(n)
    
    if showErrors==False and ret==False:        
        print("{0} MATISSE Reductions found in {1}".format(nred,fname))
        print("Current Status:")
        for i in range(nred):
            print("{0}/{1}\t: err={2}".format(n[i],ntot[i],nerr[i]))
           
        print("-------------------------------")
        print("Total of {0} recipes launched".format(sum(n)))
        #print("Total of {0} recipes launched" in {1}min ({2} min/recipes)".format(sum(n),dt,dt/sum(n)))
    elif showErrors==True:
        for i in range(nred):
                for errj in err[i]:
                    txt=errj.split(" ")[-1]
                    print(txt)    
    else:
        return (n,ntot,nerr)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read Matisse reduction log')
    
    parser.add_argument('filename', metavar='filename', type=str, \
        help='Name of the log file')   
    parser.add_argument('--showErrors', default=0, \
        help='Return list of failed reduction blocks', action='store_true')  
    
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_readReductionLog.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        sys.exit(0)
    
     
    mat_readReductionLog(args.filename,args.showErrors)           