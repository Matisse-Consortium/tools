#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on August 2018
@author: fmillour
"""

import glob
from   astropy.io import fits as fits
import wx
import os;
import argparse
import numpy as np
from mat_fileDialog import mat_FileDialog
from mat_show_oifits import open_oi,open_oi_dir,filter_oi_list
from   subprocess import call 


fileName = "/home/fmillour/reductionQueue.txt"
rawDir   = "/data/RawDataMatisse/"
destDir  = "/data/users/fmillour/DATA/"

###############################################################################

if __name__ == '__main__':
  print("Starting ...")

  #--------------------------------------------------------------------------
  parser = argparse.ArgumentParser(description='Queue processing.')
  #--------------------------------------------------------------------------
  parser.add_argument('--nbCore', default=1,type=int, \
                      help='Number Of Cores (default 1)')
  
  try:
    args = parser.parse_args()
  except:
    print("\n\033[93mRunning --help to be kind with you:\033[0m\n")
    parser.print_help()

  outName = os.path.splitext(fileName)[0]+"_Done"+os.path.splitext(fileName)[1];
    
  with open(fileName, "r") as f:
    lineList = f.readlines()
    
  for i,iline in enumerate(lineList):
    lineList[i] = iline.strip("\n")
      
  print(lineList);

  while("Look into the file"):
    with open(outName, "a+") as fo:
      with open(fileName, "r+") as f:
        lineList = f.readlines()
        f.seek(0)
        try:
          iline = lineList[0].strip("\n");
          print(iline)
                
          cmd = "export date="+iline+" ; export user=`whoami` ; echo 'export date=$date ; mat_autoPipeline.py "+rawDir+"$date --dirResult="+destDir+"$date --nbCore="+str(args.nbCore)+" --maxIter=2 ; echo \"Dear $user, the process for $date is done!\" | mail -s \"Process done $date\" $user@oca.eu' > run_$date.sh ; chmod +x run_$date.sh ; ./run_$date.sh > run_$date.log;"
        
          print(cmd);
          try:
            call(cmd, shell=True)
          
            for i in lineList[1:]:
              f.write(i)
            f.truncate()
            
            fo.write(iline+"\n")
          
          except:
            f.close()
            fo.close()
            print("autopipeline crashed!")
      
          cmd = "mat_tidyupOiFits.py "+destDir+iline+" > run_tidy_$date.log;"
          print(cmd);
          try:
            call(cmd, shell=True)
          except:
            f.close()
            fo.close()
            print("tidyup crashed!")

          cmd = "mat_autoCalib.py "+destDir+iline+"_OIFITS -o "+destDir+iline+"_OIFITS/../CALIBRATED > run_calib_$date.log;"
          print(cmd);
          try:
            call(cmd, shell=True)
          except:
            f.close()
            fo.close()
            print("autocalib crashed!")
          
        except:
          f.close()
          fo.close()
          print("No more line to process in "+fileName)
          break
        
  f.close()
  fo.close()
  print("Job done!")
